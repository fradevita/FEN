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
    character(len=3)  :: ss
    character(len=11) :: filename
    character(len=1)  :: test_case

    ! Initialize the domain decomposition
    call mpi_init(ierror)
    
    ! The domain is a unit squared box
    Lx = 1.0_dp
    Ly = 1.0_dp

    ! Select case
    call get_command_argument(1, test_case)

    ! Set the resolution
    select case(test_case)
    case('1')
        Nx = 50
        Ny = 50
    case('2')
        Nx = 100
        Ny = 100
    case('3')
        Nx = 200
        Ny = 200
    end select

    Nz = 1
    Lz = Lx*float(Nz)/float(Nx)
    origin = [0.0_dp, 0.0_dp, 0.0_dp]

    ! Create the grid
    call comp_grid%setup(Nx, Ny, Nz, Lx, Ly, Lz, origin, 1, 1)

    ! Allocate all vof fields
    call allocate_vof_fields(comp_grid)
    
    ! Set the initial VoF from a distance function
    distance => Zalesak
    call get_vof_from_distance
    
    ! Create the fixed velocity field
    call v%allocate(comp_grid, 1)
    call init_velocity()

    step = 0
    time = 0.0_dp
    dt = 0.00125_dp*pi
    Tmax = 2*pi
    Nstep = int(Tmax/dt)

    ! Write output
    write(ss,'(I0.3)') Nx
    filename = 'ini_'//ss//'.raw'
    call vof%write(filename)

    !==== Start Time loop =========================================================================
    do while(time < Tmax)

        step = step + 1
        time = time + dt

        ! Advect the VoF interface
        call advect_vof(v, dt)

        call check_vof_integral(int_phase_1, int_phase_2)

        if (comp_grid%rank == 0) write(stdout,10) 'step: ', step, 'time: ', time, 'dt: ', dt, &
        'phase 1 integral: ', int_phase_1, 'pahse 2 integral: ', int_phase_2
10 format(A6,I7,1x,A6,E13.6,1x,A4,E13.6,1x,A18,E13.6,1x,A18,E13.6)

    end do

    ! Syncornize processors
    call mpi_barrier(mpi_comm_world, ierror)

    ! Write output
    write(ss,'(I0.3)') Nx
    filename = 'fin_'//ss//'.raw'
    call vof%write(filename)

    ! Free the memory allocated by the solver
    call v%destroy()
    
    ! Finalize the simulation
    call comp_grid%destroy()
    call MPI_FINALIZE(ierror)

contains

    !==============================================================================================
    function Zalesak(x,y) result(d)

        ! In/Out variables
        real(dp), intent(in) :: x, y
        real(dp)             :: d

        ! Local variables
        real(dp) :: xc, yc, zc, R, dl(2)
        real(dp) :: x1, x2, y1, y2, dx, dy

        xc = 0.5_dp
        yc = 0.75_dp
        R = 0.15_dp
        x1 = 0.5_dp - 0.025_dp
        x2 = 0.5_dp + 0.025_dp
        y1 = 0.6_dp
        y2 = 0.85_dp

        dl(1) = -(sqrt( (x - xc)**2 + (y - yc)**2) - R)
        dx = min(dabs(x - x1),dabs(x - x2))
        dy = min(dabs(y - y1),dabs(y - y2))
        if ( dabs(x - 0.5_dp) <= 0.025_dp .and. dabs(y - 0.725) <= 0.125_dp) then
            dl(2) = -min(dx,dy)
        else
            if ( dabs(x - 0.5_dp) <= 0.025_dp) then
                dl(2) = dy
            elseif ( dabs(y - 0.725_dp) <= 0.125_dp) then
                dl(2) = dx
            else
                dl(2) = sqrt(dx**2 + dy**2)
            endif
        endif
        d = -minval(dl)

    end function Zalesak
    !==============================================================================================

    !==============================================================================================
    subroutine init_velocity()
        
        integer  :: i, j, k
        real(dp) :: x, y

        do k = comp_grid%lo(3),comp_grid%hi(3)
            do j = comp_grid%lo(2),comp_grid%hi(2)
                do i = comp_grid%lo(1),comp_grid%hi(1)
                x = float(i)*comp_grid%delta
                y = (float(j) - 0.5_dp)*comp_grid%delta
                v%x%f(i,j,k) = 0.5_dp - y
                x = (float(i) - 0.5_dp)*comp_grid%delta
                y = float(j)*comp_grid%delta
                v%y%f(i,j,k) = x - 0.5_dp
                end do
            end do
        end do

        ! Update halos and bc
        call v%update_ghost_nodes()

    end subroutine init_velocity
    !==============================================================================================

end program main
