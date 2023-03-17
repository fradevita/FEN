program main

  use mpi
  use precision      , only : dp
  use constants      , only : pi
  use class_Grid     , only : base_grid, bc_type
  use class_Vector
  use volume_of_fluid, only : allocate_vof_fields, get_h_from_vof, distance
  use volume_of_fluid, only : vof, get_vof_from_distance, advect_vof, check_vof_integral
  use io
  use json      
  use decomp_2d    

  implicit none

  integer           :: ierror, Nx, Ny, Nz, step, Nstep
  real(dp)          :: Lx, Ly, Lz, dt, Tmax, time, int_phase_1, int_phase_2, origin(3)
  type(vector)      :: v
  type(bc_type)     :: bc(4)
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
  bc(1)%s = 'Periodic'
  bc(2)%s = 'Periodic'
  bc(3)%s = 'Periodic'
  bc(4)%s = 'Periodic'

  ! Create the grid
  call base_grid%setup(Nx, Ny, Nz, Lx, Ly, Lz, origin, 4, 1, bc)

  ! Allocate all vof fields
  call allocate_vof_fields
  
  ! Set the initial VoF from a distance function
  distance => circle
  call get_vof_from_distance
  
  ! Create the fixed velocity field
  call v%allocate(1)
  call init_velocity()

  step = 0
  time = 0.0_dp
  dt = 0.00125_dp*pi
  Tmax = 10*pi
  Nstep = int(Tmax/dt)

  write(ss,'(I0.7)') step
  filename = 'data/vof_'//ss
  call vof%write(filename)
  call print_setup_json(dt)

  !==== Start Time loop ===================================================================
  do while(time < Tmax)

    step = step + 1
    time = time + dt

    ! Advect the VoF interface
    call advect_vof(v, dt)

    if (step == Nstep/2) then
      v%x%f = -v%x%f
      v%y%f = -v%y%f
      call v%apply_bc()
    endif

    if (mod(step, 100) == 0) then
      write(ss,'(I0.7)') step
      filename = 'data/vof_'//ss
      call vof%write(filename)
    endif

    call check_vof_integral(int_phase_1, int_phase_2)

    if (base_grid%rank == 0) write(stdout,10) 'step: ', step, 'time: ', time, 'dt: ', dt, &
      'phase 1 integral: ', int_phase_1, 'pahse 2 integral: ', int_phase_2
10 format(A6,I7,1x,A6,E13.6,1x,A4,E13.6,1x,A18,E13.6,1x,A18,E13.6)

  end do

  ! Syncornize processors
  call mpi_barrier(mpi_comm_world, ierror)
  
  ! Free the memory allocated by the solver
  call v%destroy()
  
  ! Finalize the simulation
  call decomp_2d_finalize
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

    do k = base_grid%lo(3),base_grid%hi(3)
      do j = base_grid%lo(2),base_grid%hi(2)
        do i = base_grid%lo(1),base_grid%hi(1)
          x = float(i)*base_grid%delta
          y = (float(j) - 0.5_dp)*base_grid%delta
          v%x%f(i,j,k) = sin(x)*cos(y)
          x = (float(i) - 0.5_dp)*base_grid%delta
          y = float(j)*base_grid%delta
          v%y%f(i,j,k) = -cos(x)*sin(y)
        end do
      end do
    end do

    ! Update halos and bc
    call v%apply_bc()

  end subroutine init_velocity
  !========================================================================================

end program main
