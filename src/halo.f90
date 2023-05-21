!< Module to define procedures for halo update
module halo_mod

    implicit none

    private
    public :: update_halos

contains

    !==============================================================================================
    subroutine update_halos(f, G, l)

        use precision_mod, only : dp
        use grid_mod     , only : grid
        use decomp_2d    , only : update_halo

        ! In/Out variables
        integer   , intent(in   ) :: l                        !< ghost nodes level
        type(grid), intent(in   ) :: G                        !< Grid where f is defined
        real(dp)  , intent(inout) :: f(G%lo(1)-l:G%hi(1)+l, & !< Input array
                                       G%lo(2)-l:G%hi(2)+l, &
                                       G%lo(3)-l:G%hi(3)+l)

        ! Local variables
        integer :: lo(3), hi(3)
        real(dp), dimension(:,:,:), allocatable :: fh

        lo = G%lo
        hi = G%hi
        
        ! Call decomp_2d function to update halos
        call update_halo(f(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), fh, level = l, opt_global = .true.)

        ! Copy into f
        f(lo(1):hi(1),lo(2)-l:hi(2)+l,lo(3)-l:hi(3)+l) = fh(:,lo(2)-l:hi(2)+l,lo(3)-l:hi(3)+l)

        ! Free memroy
        deallocate(fh)

    end subroutine update_halos
    !==============================================================================================

end module
