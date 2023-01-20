module non_newtonian

   ! Module containing procedures for the computation of viscosity for non-newtonian simulations.

   use precision, only : dp

   implicit none

   real(dp) :: mu_max, mu_min, n, kc

   interface
      subroutine viscosity_model(mu, D)
         use class_Scalar
         use class_Tensor
         type(tensor), intent(in   ) :: D
         type(scalar), intent(inout) :: mu
      end subroutine viscosity_model
   end interface
   procedure(viscosity_model), pointer :: compute_viscosity => power_law

contains

   !========================================================================================
   subroutine power_law(mu, D)

      use class_Grid  , only : base_grid
      use class_Scalar
      use class_Tensor, only : tensor
      use constants   , only : small

      ! Compute the local viscosity using the power law

      ! In/Out variables
      type(tensor), intent(in   ) :: D
      type(scalar), intent(inout) :: mu

      ! Local variables
      integer :: i, j, k, lo(3), hi(3)
      real(dp) :: gamma_dot

      lo = base_grid%lo
      hi = base_grid%hi

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               gamma_dot = (4.0_dp*D%x%x%f(i,j,k)**2 + 4.0_dp*D%y%y%f(i,j,k)**2 + &
                            8.0_dp*D%x%y%f(i,j,k)**2)*0.5_dp

               if (gamma_dot < 1.0e-7_dp) then
                  mu%f(i,j,k) = mu_min
               else
                  mu%f(i,j,k) = kc*gamma_dot**((n - 1.0_dp)/2.0_dp)
               endif
            end do
         end do
      end do

      call mu%apply_bc()

   end subroutine power_law
   !========================================================================================

end module non_newtonian
