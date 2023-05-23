program main

    use mpi
    use precision_mod       , only : dp
    use global_mod          , only : ierror
    use grid_mod
    use volume_of_fluid_mod , only : distance, vof
    use multiphase_mod      , only : rho_0, rho_1, mu_0, mu_1, sigma
    use navier_stokes_mod   , only : g, set_timestep, v
    use solver_mod
    use IO_mod
    
    implicit none

    integer          :: Nx, Ny, Nz, step, out_id
    real(dp)         :: Lx, Ly, Lz, time, dt, Tmax, origin(3)
    type(grid)       :: comp_grid
    type(bc_type)    :: bc(4)
    character(len=1) :: test_case
    
    ! Initialize the MPI library
    call mpi_init(ierror)

    call get_command_argument(1, test_case)

    ! Create the grid
    Lx = 1.0_dp
    Ly = 2.0_dp
    select case(test_case) 
    case('1')
        Nx = 64
        Ny = 128
    case('2')
        Nx = 128
        Ny = 256
    end select
    Nz = 1
    Lz = Lx*float(Nz)/float(Nx)
    origin = [0.0_dp, 0.0_dp, 0.0_dp]
    bc(1)%s = 'Periodic'
    bc(2)%s = 'Periodic'
    bc(3)%s = 'Wall'
    bc(4)%s = 'Wall'
    call comp_grid%setup(Nx, Ny, Nz, Lx, Ly, Lz, origin, 4, 1, bc)
    comp_grid%name = 'grid'

    ! Setup the initial VoF from a distance function
    distance => circle

    ! Set the material properties

    select case(test_case) 
    case('1')
        mu_0 = 10.0_dp
        mu_1 = 1.0_dp
        rho_0 = 1000.0_dp
        rho_1 = 100.0_dp
        sigma = 24.5_dp
    case('2')
        mu_0 = 10.0_dp
        mu_1 = 0.1_dp
        rho_0 = 1000.0_dp
        rho_1 = 1.0_dp
        sigma = 1.96_dp
    end select

    ! Set gravity
    g(2) = -0.98_dp
    
    ! Initialize the solver
    call init_solver(comp_grid)
    step = 0
    time = 0.0_dp
    Tmax = 3.0_dp
    call set_timestep(comp_grid, dt, 0.25_dp)

    ! Open output file
    open(newunit = out_id, file = 'output.txt') 

    !==== Start Time loop ===================================================================
    do while(time <= Tmax)

        step = step + 1
        time = time + dt

        ! Advance in time the solution
        call advance_solution(comp_grid, step, dt)

        ! Print solver status on log file
        call print_solver_status(stdout, step, time, dt)
        
        ! Compute Center of Mass quantites
        call point_quantities()
    end do

    ! Print final vof
    call vof%write('vof.dat')

    ! Free the memory allocated by the solver
    call destroy_solver()

    ! Finalize the simulation
    !call decomp_2d_finalize
    call MPI_FINALIZE(ierror)

contains

    !========================================================================================
    function circle(x,y) result(d)

        ! In/Out variables
        real(dp), intent(in) :: x, y

        ! Local variables
        real(dp) :: x0, y0, r, d

        x0 = 0.5_dp
        y0 = 0.5_dp
        r  = 0.25_dp
        d = -(sqrt((x - x0)**2 + (y - y0)**2) - r)

    end function circle
    !========================================================================================

    !========================================================================================
    subroutine point_quantities()

        integer  :: i, j, k
        real(dp) :: int_vof, xc, uc, x, delta

        int_vof = 0._dp
        xc = 0._dp
        uc = 0._dp
        delta = comp_grid%delta
        
        do k = comp_grid%lo(3),comp_grid%hi(3)
            do j = comp_grid%lo(2),comp_grid%hi(2)
                x = (j - 0.5)*delta
                do i = comp_grid%lo(1),comp_grid%hi(1)
                xc = xc + x*vof%f(i,j,k)*delta*delta
                uc = uc + 0.5_dp*v%y%f(i,j,k)*(vof%f(i,j+1,k) + vof%f(i,j,k))*delta*delta
                int_vof = int_vof + vof%f(i,j,k)*delta*delta
                end do
            end do
        end do
        
        if (comp_grid%nranks > 1) then
            call mpi_allreduce(mpi_in_place,xc,1,mpi_real8,mpi_sum,mpi_comm_world,ierror)
            call mpi_allreduce(mpi_in_place,uc,1,mpi_real8,mpi_sum,mpi_comm_world,ierror)
            call mpi_allreduce(mpi_in_place,int_vof,1,mpi_real8,mpi_sum,mpi_comm_world,ierror)
        endif

        if (comp_grid%rank == 0) then
            write(out_id,*) time, xc/int_vof, uc/int_vof, int_vof
            flush(out_id)
        endif

    end subroutine point_quantities
    !========================================================================================

end program main
