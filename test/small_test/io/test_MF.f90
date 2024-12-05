program test_io_mf

    use mpi
    use precision_mod       , only : dp
    use global_mod          , only : ierror, myrank, pi
    use grid_mod            , only : grid
    use multiphase_mod      , only : rho_0, rho_1, mu_0, mu_1, sigma, p_o
    use volume_of_fluid_mod , only : distance, vof, x_first
    use navier_stokes_mod
    use solver_mod
    use IO_mod

    implicit none

    ! Parameters
    real(dp), parameter :: a = 0.01_dp
    real(dp), parameter :: lambda = 1.0_dp
    real(dp), parameter :: rho_w = 1.0_dp
    real(dp), parameter :: wn = 2.0_dp*pi/lambda

    ! Variables
    integer          :: Nx, Ny, Nz, step, i, j
    real(dp)         :: Lx, Ly, Lz, time, dt, origin(3)
    type(grid)       :: comp_grid
    character(len=1) :: case

    ! Init MPI
    call mpi_init(ierror)
    call mpi_comm_rank(mpi_comm_world, myrank, ierror)

    ! Create the numerical grid
    Nx = 16
    Ny = 16*3
    Nz = 1
    Lx = lambda
    Ly = 3.0_dp*Lx
    Lz = Lx*float(Nz)/float(Nx)
    origin = [0.0_dp, -Ly/2.0_dp, 0.0_dp]
    call comp_grid%setup(Nx, Ny, Nz, Lx, Ly, Lz, origin, 1, 1)

    ! Setup the multiphase parameters
    rho_0 = rho_w
    rho_1 = rho_w
    mu_0 = 0.0182571749236_dp
    mu_1 = mu_0
    sigma = 1.0_dp
    
    ! Set the initial interface
    distance => wave

    ! Initialize the solver
    call init_solver(comp_grid)
    step = 0
    time = 0.0_dp
    call set_timestep(comp_grid, dt, 1.0_dp)

    call get_command_argument(1, case)
    select case(case)
    case('1')
        ! Perform 10 timestep
        do i = 1,10
            step = step + 1
            time = time + dt
            call advance_solution(comp_grid, step, dt)

            ! Save simulation state
            if (i == 5) call save_state(step)
        end do
        call save_state(step)
    case('2')
        ! Now load state 5 and perform again 4 timestep
        step = 5
        call load_state(step)

        ! Set appropriate vof advection order to make comparison
        if (x_first .eqv. .false.) then
            x_first = .true.
        else
            x_first = .false.
        endif

        do i = 6, 10
            step = step + 1
            time = time + dt
            call advance_solution(comp_grid, step, dt)
        end do
        call save_state(step*2)
    end select

    ! Clean memory
    call destroy_solver()
    call comp_grid%destroy()
    call mpi_finalize(ierror)

contains

    !===============================================================================================
    function wave(x, y) result (d)

        ! In/Out variables
        real(dp), intent(in) :: x, y

        ! Local variables
        real(dp) :: x1, y1, x2, y2, d

        x1 = x - comp_grid%delta/2.0_dp
        y1 = a*cos(wn*x1)

        x2 = x + comp_grid%delta/2.0_dp
        y2 = a*cos(wn*x2)
        
        d = -((x2 - x1)*(y1 - y) - (x1 - x)*(y2 - y1))/sqrt((x2 - x1)**2 + (y2 - y1)**2)

    end function
    !===============================================================================================

end program test_io_mf