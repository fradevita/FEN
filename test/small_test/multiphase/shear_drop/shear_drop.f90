program main

    use mpi
    use precision_mod       , only : dp
    use grid_mod            , only : grid, bc_type
    use volume_of_fluid_mod , only : distance, beta
    use multiphase_mod      , only : rho_0, rho_1, mu_0, mu_1, sigma
    use navier_stokes_mod   , only : set_timestep, v, dt_o
    use solver_mod          , only : init_solver, advance_solution, destroy_solver
    use solver_mod          , only : print_solver_status, save_fields
    use IO_mod

    implicit none

    ! Simulation parameters
    real(dp), parameter :: H = 1.0_dp    ! Lenght scale (half channel width)
    real(dp), parameter :: U = 1.0_dp    ! Wall velocity
    real(dp), parameter :: aOvH = 0.5_dp ! droplet size  over channel height 
    real(dp), parameter :: a = aOvH*H    ! droplet radius
    real(dp), parameter :: lm = 1.0_dp   ! viscosity ratio: mu_d/mu_f
    real(dp), parameter :: Re = 1.0_dp   ! Reynolds number = rho U 2a/mu
#if CASE==1
    real(dp), parameter :: Ca = 0.2_dp   ! Capillary number = U mu / sigma
#elif CASE==2
    real(dp), parameter :: Ca = 0.4_dp   ! Capillary number = U mu / sigma
#else
    real(dp), parameter :: Ca = 0.9_dp   ! Capillary number = U mu / sigma
#endif

    ! Variables  
    integer       :: ierror, Nx, Ny, Nz, step, out_id
    real(dp)      :: Lx, Ly, Lz, time, dt, Tmax, origin(3)
    type(grid)    :: comp_grid
    type(bc_type) :: bc(4)

    ! Initialize the MPI library
    call mpi_init(ierror)

    ! Create the grid
    Lx = 2.0_dp*H
    Ly = 2.0_dp*H
    Nx = 64
    Ny = 64
    Nz = 1
    Lz = Lx*float(Nz)/float(Nx)
    origin = [0.0_dp, 0.0_dp, 0.0_dp]
    bc(1)%s = 'Periodic'
    bc(2)%s = 'Periodic'
    bc(3)%s = 'Wall'
    bc(4)%s = 'Wall'
    call comp_grid%setup(Nx, Ny, Nz, Lx, Ly, Lz, origin, 8, 1, bc)

    ! Setup the initial VoF from a distance function
    beta = 1.0_dp      ! This value gives better results for this  test
    distance => circle

    ! Set the material properties
    rho_0 = 1.0_dp                     ! Neutrally buoyant
    rho_1 = rho_0                      ! Neutrally buoyant
    mu_0 = rho_0 * U * 2.0_dp * a / Re ! Set fluid viscosity with Reynolds number 
    mu_1 = lm * mu_0                   ! Set droplet viscosity with viscosity ratio
    sigma = U * mu_0 / Ca              ! Set surface tension with Capillary number

    ! Initialize the solver
    call init_solver(comp_grid)
    step = 0
    time = 0.0_dp
#if CASE==1
    Tmax = 1.0_dp
#elif CASE==2
    Tmax = 2.0_dp
#else
    Tmax = 3.0_dp
#endif 

    call set_timestep(comp_grid, dt, U)
    v%x%bc%top =  U
    v%x%bc%bottom = -U

    block
        integer :: i, j, k
        do k = comp_grid%lo(3),comp_grid%hi(3)
            do j = comp_grid%lo(2),comp_grid%hi(2)
                do i = comp_grid%lo(1),comp_grid%hi(1)
                    v%x%f(i,j,k) = -U + 2.0_dp*U*comp_grid%y(j)/Ly
                end do
            end do
        end do
    end block
    call v%update_ghost_nodes()

    call save_fields(0)

    !==== Start Time loop ===================================================================
    do while(time <= Tmax)

        step = step + 1
        time = time + dt

        ! Advance in time the solution
        call advance_solution(comp_grid, step, dt)

        ! Print solver status on log file
        call print_solver_status(stdout, step, time, dt)
        
        ! Compute Center of Mass quantites
        if (mod(step,100) == 0) call save_fields(step)
    end do

    ! Free the memory allocated by the solver
    call destroy_solver()

    ! Finalize the simulation
    call MPI_FINALIZE(ierror)

contains

    !========================================================================================
    function circle(x,y) result(d)

        ! In/Out variables
        real(dp), intent(in) :: x, y

        ! Local variables
        real(dp) :: x0, y0, d

        x0 = 1.0_dp
        y0 = 1.0_dp
        d = -(sqrt((x - x0)**2 + (y - y0)**2) - a)

    end function circle
    !========================================================================================

end program main
