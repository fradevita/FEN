program test_fields

  use mpi
  use precision , only : dp
  use constants , only : pi
  use class_Grid, only : base_grid, bc_type
  use class_Scalar
  use class_Vector
  use tensors   , only : tensor

  implicit none

  ! Parameters
  real(dp), parameter :: tol = 1.0e-15_dp

  ! Variables
  integer       :: ierror, Nx, Ny, Nz, prow, pcol, i, j, k
  real(dp)      :: Lx, Ly, Lz, e, x, y, z, origin(3)
  type(bc_type) :: bc(6)
  type(scalar)  :: s
  type(vector)  :: v
  type(tensor)  :: T

  ! Init the MPI
  call MPI_init(ierror)

  ! Set the grid
  Nx = 64
  Ny = 64
  Nz = 64
  Lx = 1.0_dp
  Ly = 1.0_dp
  Lz = 1.0_dp
  origin = [0.0_dp, 0.0_dp, 0.0_dp]
  prow = 1
  pcol = 8
  bc(1)%s = 'Periodic'
  bc(2)%s = 'Periodic'
  bc(3)%s = 'Periodic'
  bc(4)%s = 'Periodic'
  bc(5)%s = 'Periodic'
  bc(6)%s = 'Periodic'
  call base_grid%setup(Nx, Ny, Nz, Lx, Ly, Lz, origin, prow, pcol, bc)

  ! Allocate memory for the scalar field s
  call s%allocate(1)

  ! Initialize the scalar field s
  do k = base_grid%lo(3),base_grid%hi(3)
     z = (float(k) - 0.5_dp)*base_grid%delta
     do j = base_grid%lo(2),base_grid%hi(2)
        y = (float(j) - 0.5_dp)*base_grid%delta
        do i = base_grid%lo(1),base_grid%hi(1)
           x = (float(i) - 0.5_dp)*base_grid%delta
           s%f(i,j,k) = sample_function(x, y, z)
        end do
     end do
  end do

  ! Apply periodic BC
  call s%apply_bc()

  ! Check that the ghost nodes have been correctly filled
  call check_ghost_nodes(s)

  ! Destroy the scalar s
  call s%destroy()

  ! Allocate memory for the vector field v
  call v%allocate(1)

  ! Initialize the vector field v
  do k = base_grid%lo(3),base_grid%hi(3)
     z = (float(k) - 0.5_dp)*base_grid%delta 
     do j = base_grid%lo(2),base_grid%hi(2)
        y = (float(j) - 0.5_dp)*base_grid%delta 
        do i = base_grid%lo(1),base_grid%hi(1)
           x = (float(i) - 0.5_dp)*base_grid%delta
           v%x%f(i,j,k) = sample_function(x,y,z)
           v%y%f(i,j,k) = sample_function(x,y,z)
           v%z%f(i,j,k) = sample_function(x,y,z)
        end do
     end do
  end do

  ! Apply periodic BC
  call v%apply_bc()

  ! Check that the ghost nodes have been correctly filled
  call check_ghost_nodes(v%x)
  call check_ghost_nodes(v%y)
  call check_ghost_nodes(v%z)
  
  ! Free the memory allocate for the vector v
  call v%destroy()

  ! Allocate memory for the tensor field T
  call T%allocate(1)

  ! Initialize the tensor field T
  do k = base_grid%lo(3),base_grid%hi(3)
     z = (float(k) - 0.5_dp)*base_grid%delta 
     do j = base_grid%lo(2),base_grid%hi(2)
        y = (float(j) - 0.5_dp)*base_grid%delta 
        do i = base_grid%lo(1),base_grid%hi(1)
           x = (float(i) - 0.5_dp)*base_grid%delta
           T%x%x%f(i,j,k) = sample_function(x,y,z)
           T%x%y%f(i,j,k) = sample_function(x,y,z)
           T%x%z%f(i,j,k) = sample_function(x,y,z)
           T%y%x%f(i,j,k) = sample_function(x,y,z)
           T%y%y%f(i,j,k) = sample_function(x,y,z)
           T%y%z%f(i,j,k) = sample_function(x,y,z)
           T%z%x%f(i,j,k) = sample_function(x,y,z)
           T%z%y%f(i,j,k) = sample_function(x,y,z)
           T%z%z%f(i,j,k) = sample_function(x,y,z)
        end do
     end do
  end do

  ! Apply periodic BC
  call T%apply_bc()

  ! Check that the ghost nodes have been correctly filled
  !call check_ghost_nodes(T%x%x)
  !call check_ghost_nodes(T%x%y)
  !call check_ghost_nodes(T%x%z)
  !call check_ghost_nodes(T%y%x)
  !call check_ghost_nodes(T%y%y)
  !call check_ghost_nodes(T%y%z)
  !call check_ghost_nodes(T%z%x)
  !call check_ghost_nodes(T%z%y)
  !call check_ghost_nodes(T%z%z)
  
  ! Free the memory allocate for the vector v
  call T%destroy()

  ! Destroy the grid
  call base_grid%destroy()

  ! Finalize the simulation
  call MPI_FINALIZE(ierror)

contains

  !========================================================================================
  pure real(dp) function sample_function(x, y, z)

    real(dp), intent(in) :: x, y, z

    sample_function = sin(2.0_dp*pi*x)*cos(2.0_dp*pi*y)*sin(2.0_dp*pi*z)

  end function sample_function
  !========================================================================================

  !========================================================================================
  subroutine check_ghost_nodes(s)

    use io
    
    type(scalar), intent(in) :: s

    ! Check that the ghost nodes have been correctly filled
    ! ** Left boundary **
    i = base_grid%hi(1) - 1
    x = (float(i) - 0.5_dp)*base_grid%delta
    do k = base_grid%lo(3),base_grid%hi(3)
       z = (float(k) - 0.5_dp)*base_grid%delta
       do j = base_grid%lo(2),base_grid%hi(2)
          y = (float(j) - 0.5_dp)*base_grid%delta
          e = abs(s%f(i,j,k) - sample_function(x,y,z))
          if (e > tol) then
             call print_error_message('Error greather than tollerance on left boundary')
          endif
       end do
    end do

    ! ** Right boundary **
    i = base_grid%hi(1) + 1
    x = (float(i) - 0.5_dp)*base_grid%delta
    do k = base_grid%lo(3),base_grid%hi(3)
       z = (float(k) - 0.5_dp)*base_grid%delta
       do j = base_grid%lo(2),base_grid%hi(2)
          y = (float(j) - 0.5_dp)*base_grid%delta
          e = abs(s%f(i,j,k) - sample_function(x,y,z))
          if (e > tol) then
             print *, i, j, k, e
             call print_error_message('Error greather than tollerance on right boundary')
          endif
       end do
    end do

    ! ** Top boundary **
    j = base_grid%hi(2) + 1
    y = (float(j) - 0.5_dp)*base_grid%delta
    do k = base_grid%lo(3),base_grid%hi(3)
       z = (float(k) - 0.5_dp)*base_grid%delta
       do i = base_grid%lo(1),base_grid%hi(1)
          x = (float(i) - 0.5_dp)*base_grid%delta
          e = abs(s%f(i,j,k) - sample_function(x,y,z))
          if (e > tol) then
             call print_error_message('Error greather than tollerance on top boundary')
          endif
       end do
    end do

    ! ** Bottom boundary **
    j = base_grid%lo(2) - 1
    y = (float(j) - 0.5_dp)*base_grid%delta
    do k = base_grid%lo(3),base_grid%hi(3)
       z = (float(k) - 0.5_dp)*base_grid%delta 
       do i = base_grid%lo(1),base_grid%hi(1)
          x = (float(i) - 0.5_dp)*base_grid%delta
          e = abs(s%f(i,j,k) - sample_function(x,y,z))
          if (e > tol) then
             call print_error_message('Error greather than tollerance on bottom boundary')
          endif
       end do
    end do

    ! ** Front boundary **
    k = base_grid%hi(3) + 1
    z = (float(k) - 0.5_dp)*base_grid%delta
    do j = base_grid%lo(2),base_grid%hi(2)
       y = (float(j) - 0.5_dp)*base_grid%delta
       do i = base_grid%lo(1),base_grid%hi(1)
          x = (float(i) - 0.5_dp)*base_grid%delta
          e = abs(s%f(i,j,k) - sample_function(x,y,z))
          if (e > tol) then
             call print_error_message('Error greather than tollerance on front boundary')
          endif
       end do
    end do

    ! ** Back boundary **
    k = base_grid%lo(3) - 1
    z = (float(k) - 0.5_dp)*base_grid%delta
    do j = base_grid%lo(2),base_grid%hi(2)
       y = (float(j) - 0.5_dp)*base_grid%delta
       do i = base_grid%lo(1),base_grid%hi(1)
          x = (float(i) - 0.5_dp)*base_grid%delta
          e = abs(s%f(i,j,k) - sample_function(x,y,z))
          if (e > tol) then
             call print_error_message('Error greather than tollerance on back boundary')
          endif
       end do
    end do
    
  end subroutine check_ghost_nodes
  !========================================================================================
  
end program test_fields
