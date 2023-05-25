program main

    use mpi
    use precision_mod      , only : dp
    use global_mod         , only : pi
    use grid_mod
    use vector_mod
    use volume_of_fluid_mod
    use IO_mod

    implicit none

    integer           :: ierror, Nx, Ny, Nz, step, Nstep
    real(dp)          :: Lx, Ly, Lz, dt, Tmax, time, int_phase_1, int_phase_2, origin(3)
    type(grid)        :: comp_grid
    type(vector)      :: v
    character(len=7)  :: ss
    character(len=16) :: filename

    ! Initialize the domain decomposition
    call mpi_init(ierror)

    ! The domain is a squared box of size pi
    Lx = pi
    Ly = pi
    
    ! Set the resolution
    Nx = 200
    Ny = 200
    Nz = 1
    Lz = Lx*float(Nz)/float(Nx)
    origin = [0.0_dp, 0.0_dp, 0.0_dp]

    ! Create the grid
    call comp_grid%setup(Nx, Ny, Nz, Lx, Ly, Lz, origin, 4, 1)

    ! Allocate all vof fields
    call allocate_vof_fields(comp_grid)
    
    ! Set the initial VoF from a distance function
    distance => circle
    call get_vof_from_distance
    
    ! Create the fixed velocity field
    call v%allocate(comp_grid, 1)
    call init_velocity()

    step = 0
    time = 0.0_dp
    dt = 0.00125_dp*pi
    Tmax = 10*pi
    Nstep = int(Tmax/dt)

    write(ss,'(I0.7)') step
    filename = 'data/vof_'//ss
    call vof%write(filename)

    !==== Start Time loop ===================================================================
    do while(time < Tmax)

        step = step + 1
        time = time + dt

        ! Advect the VoF interface
        call advect_vof(v, dt)

        if (step == Nstep/2) then
        v%x%f = -v%x%f
        v%y%f = -v%y%f
        call v%update_ghost_nodes()
        endif

        if (mod(step, 100) == 0) then
        write(ss,'(I0.7)') step
        filename = 'data/vof_'//ss
        call vof%write(filename)
        endif

        call check_vof_integral(int_phase_1, int_phase_2)

        if (comp_grid%rank == 0) write(stdout,10) 'step: ', step, 'time: ', time, 'dt: ', dt, &
        'phase 1 integral: ', int_phase_1, 'pahse 2 integral: ', int_phase_2
10 format(A6,I7,1x,A6,E13.6,1x,A4,E13.6,1x,A18,E13.6,1x,A18,E13.6)

    end do

    ! Syncornize processors
    call mpi_barrier(mpi_comm_world, ierror)
    
    ! Free the memory allocated by the solver
    call v%destroy()
    
    ! Finalize the simulation
    call comp_grid%destroy()
    call MPI_FINALIZE(ierror)

contains

    !========================================================================================
    function circle(x,y) result(d)

        ! In/Out variables
        real(dp), intent(in) :: x, y

        ! Local variables
        real(dp) :: x0, y0, r, d

        x0 = 0.5_dp*pi
        y0 = 0.2_dp*(pi + 1)
        r  = 0.2_dp*pi
        d = sqrt((x - x0)**2 + (y - y0)**2) - r

    end function circle
    !========================================================================================

    !========================================================================================
    subroutine init_velocity()
        
        integer  :: i, j, k
        real(dp) :: x, y

        do k = comp_grid%lo(3),comp_grid%hi(3)
            do j = comp_grid%lo(2),comp_grid%hi(2)
                do i = comp_grid%lo(1),comp_grid%hi(1)
                x = float(i)*comp_grid%delta
                y = (float(j) - 0.5_dp)*comp_grid%delta
                v%x%f(i,j,k) = sin(x)*cos(y)
                x = (float(i) - 0.5_dp)*comp_grid%delta
                y = float(j)*comp_grid%delta
                v%y%f(i,j,k) = -cos(x)*sin(y)
                end do
            end do
        end do

        ! Update halos and bc
        call v%update_ghost_nodes()

    end subroutine init_velocity
    !========================================================================================

end program main
