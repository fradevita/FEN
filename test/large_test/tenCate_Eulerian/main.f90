program main

    use mpi
    use precision_mod 
    use global_mod              , only : ierror, myrank, pi
    use grid_mod
    use ibm_mod
    use eulerian_sphere_mod
    use navier_stokes_mod       , only : g, set_timestep, viscosity, density, p
    use solver_mod
    use IO_mod                  , only : stdout
    use euclidean_mod           , only : ZEROV

    implicit none

    ! Parameters
    real(dp), parameter :: d = 15.0e-3_dp   ! sphere diameter       [m]
    real(dp), parameter :: Uinf = 0.128_dp  ! settling velocity   [m/s]
    real(dp), parameter :: rhop = 1120.0_dp ! particle density [kg/m^3]
    real(dp), parameter :: Lx = 10.0e-2_dp  ! Box x size
    real(dp), parameter :: Ly = 10.0e-2_dp  ! Box y size
    integer , parameter :: Nx = 128         ! Number of grid cells in x
    integer , parameter :: Ny = 128         ! Number of grid cells in y
    integer , parameter :: Nz = 208         ! Number of grid cells in z
    real(dp), parameter :: Lz = Nz*Lx/Nx    ! Box z size

    ! Variables
    integer              :: out_id, step
    real(dp)             :: rhof, muf, Tend, dt, time
    type(grid)  , target :: comp_grid
    type(bc_type)        :: bc(6)
    type(sphere), target :: S
    character(len= 1)    :: arg

    ! Setup MPI
    call mpi_init(ierror)
    call mpi_comm_rank(mpi_comm_world, myrank, ierror)

    ! Select fluid material properties based on the test case
    call get_command_argument(1, arg)
    select case(arg)
    case('1') ! Re = 1.5
        rhof = 970.0_dp
        muf  = 373.0e-3_dp
        Tend = 4.0_dp
        open(newunit = out_id, file = 'out1.csv')
        if (myrank == 0) write(out_id,'(A5)') 't,x,w'
    case('2') ! Re = 4.1
        rhof = 965.0_dp
        muf  = 212.0e-3_dp
        Tend = 2.5_dp
        open(newunit = out_id, file = 'out2.csv')
        if (myrank == 0) write(out_id,'(A5)') 't,x,w'
    case('3') ! Re = 11.6
        rhof = 962.0_dp
        muf  = 113.0e-3_dp
        Tend = 1.6_dp
        open(newunit = out_id, file = 'out3.csv')
        if (myrank == 0) write(out_id,'(A5)') 't,x,w'
    case('4') ! Re = 31.9
        rhof = 960.0_dp
        muf  = 58.0e-3_dp
        Tend = 1.2_dp
        open(newunit = out_id, file = 'out4.csv')
        if (myrank == 0) write(out_id,'(A5)') 't,x,w'
    end select

    !**** Define the computational grid ************************************************************
    bc(1)%s = 'Wall'
    bc(2)%s = 'Wall'
    bc(3)%s = 'Wall'
    bc(4)%s = 'Wall'
    bc(5)%s = 'Wall'
    bc(6)%s = 'Wall'
    call comp_grid%setup(Nx, Ny, Nz, Lx, Ly, Lz, [0.0_dp, 0.0_dp, 0.0_dp], 2, 4, bc)
    call comp_grid%print_json()
    !***********************************************************************************************

    !**** Define the solid sphere ******************************************************************
    S%density = rhof
    S%R = d/2.0_dp
    S%G => comp_grid
    call S%setup()
    allocate(solid_list(1))
    solid_list(1)%pS=> S
    S%center_of_mass%X(1:3) = [Lx/2.0_dp, Ly/2.0_dp, 120e-3_dp + d/2.0_dp]
    if (myrank == 0) write(out_id, '(*(E16.8,:,","))') 0.0_dp, S%center_of_mass%X(3), & 
                                                               S%center_of_mass%V(3)
    !***********************************************************************************************

    !**** Initilize the solver *********************************************************************
    viscosity = muf
    density = rhof
    call init_solver(comp_grid)
    step = 0
    time = 0.0_dp
    call set_timestep(comp_grid, dt, Uinf)
    dt = dt/2.0_dp

    ! Add explicilty weight to the sphere as external force
    g(3) = -9.80665_dp
    S%center_of_mass%Fe(3) = -9.80665_dp*S%volume()*(rhop - rhof)

    ! Initial pressure equal to hydrostatic pressure
    block
        integer :: i, k, j
        do k = comp_grid%lo(3),comp_grid%hi(3)
            do j = comp_grid%lo(2),comp_grid%hi(2)
                do i = comp_grid%lo(1),comp_grid%hi(1)
                    p%f(i,j,k) = -rhof*g(3)*(Lz - comp_grid%z(k))
                end do
            end do
        end do
        call p%update_ghost_nodes()
    end block
    call save_fields(step)
    !***********************************************************************************************

    !**** Start Time loop **************************************************************************
    time_loop: do while(time <= Tend)

        step = step + 1
        time = time + dt

        ! Advance in time the solution
        call advance_solution(comp_grid, step, dt)

        ! Advance solver status to log file
        call print_solver_status(stdout, step, time, dt)

        if (mod(step,100) == 0) then
            call save_fields(step)
        endif
        if (mod(step,1) == 0) then
            if (myrank == 0) then
                write(out_id, '(*(E16.8,:,","))') time, S%center_of_mass%X(3), S%center_of_mass%V(3)
                flush(out_id)
            endif
        endif

        if ( (S%center_of_mass%X(3) + d/2.0_dp) < 2*comp_grid%delta ) exit time_loop

    end do time_loop
    !***********************************************************************************************
    close(out_id)

    ! free memory
    call destroy_solver
    call comp_grid%destroy()

    call mpi_finalize(ierror)

end program main