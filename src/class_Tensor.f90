module class_Tensor

  ! This module contains the tensor type definition and procedures.
  ! Thetensor type is a type with three vectors.

  use precision, only : dp
  use class_Vector

  implicit none

  type tensor
     ! Three vector fields
     type(vector) :: x, y, z

   contains
     procedure :: allocate
     procedure :: apply_bc
     procedure :: symmetric
     procedure :: destroy     
  end type tensor

contains

  !========================================================================================
  subroutine allocate(T, l)

    ! In/Out variables
    class(tensor), intent(inout)           :: T
    integer      , intent(in   ), optional :: l

    ! By default the field has zero ghost nodes, if l is present set the number of
    ! ghost nodes per side equal to l
    if (present(l)) then
       T%x%x%gl = l
       T%x%y%gl = l
       T%x%z%gl = l
       T%y%x%gl = l
       T%y%y%gl = l
       T%y%z%gl = l
       T%z%x%gl = l
       T%z%y%gl = l
       T%z%z%gl = l
    endif

    call T%x%allocate(T%x%x%gl)
    call T%y%allocate(T%x%x%gl)
#if DIM == 3
    call T%z%allocate(T%x%x%gl)
#endif

  end subroutine allocate
  !========================================================================================

  !========================================================================================
  subroutine apply_bc(T)

    use mpi
    use class_Grid, only : base_grid
    use halo      , only : field_halo_update

    ! In/Out variables
    class(tensor), intent(inout) :: T

    ! Local variables
    integer :: ierror, l

    ! First update halos
    call field_halo_update(T%x%x%f, T%x%x%gl)
    call field_halo_update(T%x%y%f, T%x%x%gl)
    call field_halo_update(T%y%x%f, T%x%x%gl)
    call field_halo_update(T%y%y%f, T%x%x%gl)
#if DIM==3
    call field_halo_update(T%x%z%f, T%x%x%gl)
    call field_halo_update(T%y%z%f, T%x%x%gl)
    call field_halo_update(T%z%x%f, T%x%x%gl)
    call field_halo_update(T%z%y%f, T%x%x%gl)
    call field_halo_update(T%z%z%f, T%x%x%gl)
#endif

    if (base_grid%periodic_bc(1)) then
       do l = 1,T%x%x%gl
          T%x%x%f(base_grid%lo(1)-l,:,:) = T%x%x%f(base_grid%hi(1)-l+1,:,:)
          T%x%x%f(base_grid%hi(1)+l,:,:) = T%x%x%f(base_grid%lo(1)+l-1,:,:)
          T%x%y%f(base_grid%lo(1)-l,:,:) = T%x%y%f(base_grid%hi(1)-l+1,:,:)
          T%x%y%f(base_grid%hi(1)+l,:,:) = T%x%y%f(base_grid%lo(1)+l-1,:,:)
          T%y%x%f(base_grid%lo(1)-l,:,:) = T%y%x%f(base_grid%hi(1)-l+1,:,:)
          T%y%x%f(base_grid%hi(1)+l,:,:) = T%y%x%f(base_grid%lo(1)+l-1,:,:)
          T%y%y%f(base_grid%lo(1)-l,:,:) = T%y%y%f(base_grid%hi(1)-l+1,:,:)
          T%y%y%f(base_grid%hi(1)+l,:,:) = T%y%y%f(base_grid%lo(1)+l-1,:,:)
#if DIM==3
          T%x%z%f(base_grid%lo(1)-l,:,:) = T%x%z%f(base_grid%hi(1)-l+1,:,:)
          T%x%z%f(base_grid%hi(1)+l,:,:) = T%x%z%f(base_grid%lo(1)+l-1,:,:)
          T%y%z%f(base_grid%lo(1)-l,:,:) = T%y%z%f(base_grid%hi(1)-l+1,:,:)
          T%y%z%f(base_grid%hi(1)+l,:,:) = T%y%z%f(base_grid%lo(1)+l-1,:,:)
          T%z%x%f(base_grid%lo(1)-l,:,:) = T%z%x%f(base_grid%hi(1)-l+1,:,:)
          T%z%x%f(base_grid%hi(1)+l,:,:) = T%z%x%f(base_grid%lo(1)+l-1,:,:)
          T%z%y%f(base_grid%lo(1)-l,:,:) = T%z%y%f(base_grid%hi(1)-l+1,:,:)
          T%z%y%f(base_grid%hi(1)+l,:,:) = T%z%y%f(base_grid%lo(1)+l-1,:,:)
          T%z%z%f(base_grid%lo(1)-l,:,:) = T%z%z%f(base_grid%hi(1)-l+1,:,:)
          T%z%z%f(base_grid%hi(1)+l,:,:) = T%z%z%f(base_grid%lo(1)+l-1,:,:)
#endif
       end do
    endif
    
    ! This is necessary to avoit deadlock
   call mpi_barrier(mpi_comm_world, ierror)
    
  end subroutine apply_bc
  !========================================================================================

  !========================================================================================
  subroutine symmetric(T)

    ! Replace T with its symmetric tensor

    ! In/Out variables
    class(tensor), intent(inout) :: T

    ! Local variables
    type(tensor) :: T_temp

    ! Allocate memory
    call T_temp%allocate(1)

    ! Save T
    T_temp = T

    ! Replace with its symmetric tensors
    T%x%y%f = 0.5_dp*(T_temp%x%y%f + T_temp%y%x%f)
    T%y%x%f = T%x%y%f
#if DIM==3
    T%x%z%f = 0.5_dp*(T_temp%x%z%f + T_temp%z%x%f)
    T%z%x%f = T%x%z%f
    T%y%z%f = 0.5_dp*(T_temp%y%z%f + T_temp%z%y%f)
    T%z%y%f = T%y%z%f
#endif

    call T_temp%destroy
    call T%apply_bc()
    
  end subroutine symmetric
  !========================================================================================
  
  !========================================================================================
  subroutine destroy(T)

    class(tensor), intent(inout) :: T

    ! Free all memory allocated by vector field v
    call T%x%destroy()
    call T%y%destroy()
#if DIM == 3
    call T%z%destroy()
#endif

  end subroutine destroy
  !========================================================================================

end module class_Tensor
