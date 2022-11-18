!> This module contains definitions and procedures for the Immersed Boundary Method
module ibm

   use precision           , only : dp
   use class_Vector        , only : vector
   use class_Eulerian_Solid, only : eulerian_Solid_pointer

   implicit none

   type(vector)                              :: Fe                     !< Eulerian forcing vector field
   type(Eulerian_Solid_pointer), allocatable :: Eulerian_Solid_list(:) !< List of pointer to solid objects

contains

   !========================================================================================
   subroutine init_ibm

        ! Allocate memory for the forcing field that will be added to the RHS of the
        ! momentum equation.
        call Fe%allocate()

   end subroutine init_ibm
   !========================================================================================

    !========================================================================================
   subroutine compute_ibm_forcing(v, RHS, dt, F)

        use eulerian_ibm        , only : compute_eulerian_ibm_forcing => compute_ibm_forcing

        ! In/Out variables
        type(vector), intent(in   ) :: v   !< Velocity field to be forced
        type(vector), intent(in   ) :: RHS !< RHS of the momentum equation
        real(dp)    , intent(in   ) :: dt  !< timestep
        type(vector), intent(inout) :: F   !< forcing field

        F%x%f = 0.0_dp
        F%y%f = 0.0_dp
#if DIM==3
        F%x%f = 0.0_dp
#endif

        ! Evaluate the forcing due to eulerian solids
        call compute_eulerian_ibm_forcing(v, RHS, eulerian_Solid_list, dt, F)

   end subroutine compute_ibm_forcing
   !========================================================================================

   !========================================================================================
   subroutine destroy_ibm
    
        deallocate(Eulerian_Solid_list)

   end subroutine destroy_ibm
   !========================================================================================

end module ibm
