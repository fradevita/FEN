module halo

  implicit none

  private
  public :: field_halo_update

contains

  !========================================================================================
  subroutine field_halo_update(f, l)

    use precision , only : dp
    use class_Grid, only : base_grid
    use decomp_2d , only : update_halo

    ! In/Out variables
    integer   , intent(in   ) :: l
    real(dp)  , intent(inout) :: f(base_grid%lo(1)-l:base_grid%hi(1)+l, &
                                   base_grid%lo(2)-l:base_grid%hi(2)+l, &
                                   base_grid%lo(3)-l:base_grid%hi(3)+l)

    ! Local variables
    integer :: lo(3), hi(3)
    real(dp), dimension(:,:,:), allocatable :: fh

    lo = base_grid%lo
    hi = base_grid%hi
    
    ! Call decomp_2d function to update halos
    call update_halo(f(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), fh, level = l, opt_global = .true.)

    ! Copy into f
    f(lo(1):hi(1),lo(2)-l:hi(2)+l,lo(3)-l:hi(3)+l) = fh(:,lo(2)-l:hi(2)+l,lo(3)-l:hi(3)+l)

    ! Free memroy
    deallocate(fh)

  end subroutine field_halo_update
  !========================================================================================

end module halo
