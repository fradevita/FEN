!> This module contains definitions and procedures for the Immersed Boundary Method
module ibm_mod

    use precision_mod       , only : dp
    use vector_mod          , only : vector
    use solid_mod           , only : solid_pointer
    use eulerian_solid_mod  , only : eulerian_solid, eulerian_Solid_pointer
    use lagrangian_solid_mod, only : lagrangian_solid, lagrangian_solid_pointer

    implicit none

    type(vector)                                :: Fe                       !< Eulerian forcing vector field
    class(solid_pointer)          , allocatable :: solid_list(:)            !< List of pointer to solid objects
    class(eulerian_solid_pointer) , allocatable :: eulerian_solid_list(:)   !< List of pointer to eulerian solid objects
    type(lagrangian_Solid_pointer), allocatable :: lagrangian_solid_list(:) !< List of pointer to lagrangian solid objects

contains

    !==============================================================================================
    subroutine init_ibm(comp_grid)

        !use eulerian_ibm
        use grid_mod
        use eulerian_solid_mod  , only : eulerian_solid
        use eulerian_ibm_mod    , only : init_eulerian_ibm
        use lagrangian_solid_mod, only : lagrangian_solid
        use lagrangian_ibm_mod  , only : init_lagrangian_ibm => init_ibm
        
        ! In/Out variables
        type(grid), intent(in) :: comp_grid

        ! Local variables
        integer :: n, nls, nes

        ! Allocate memory for the forcing field that will be added to the RHS of the
        ! momentum equation.
        call Fe%allocate(comp_grid)

        ! Cycle over solid list to find number of lagrangian and eulerian solid
        nes = 0
        do n = 1,size(solid_list)
            select type(var => solid_list(n)%pS)
            class is (eulerian_solid)
                nes = nes + 1
            class is (lagrangian_solid)
                nls = nls + 1
            end select
        end do
        if (nes > 0) allocate(eulerian_solid_list(nes))
        nes = 1
        do n = 1,size(solid_list)
            select type(var => solid_list(n)%pS)
            class is (eulerian_solid)
                eulerian_solid_list(nes)%pS => var
                nes = nes + 1
            end select
        end do
        if (nls > 0) allocate(lagrangian_solid_list(nes))
        nls = 1
        do n = 1,size(solid_list)
            select type(var => solid_list(n)%pS)
            class is (lagrangian_solid)
                lagrangian_solid_list(nls)%pS => var
                nls = nls + 1
            end select
        end do

        ! If the list of eulerian solid has been created, initialize the eulerian ibm variables 
        if (allocated(eulerian_solid_list)) call init_eulerian_ibm(eulerian_solid_list, comp_grid)

        ! If the list of lagrangian solid has been created, initialize the lagrangian ibm variables 
        if (allocated(lagrangian_solid_list)) call init_lagrangian_ibm(comp_grid)

    end subroutine init_ibm
    !==============================================================================================

    !==============================================================================================
    subroutine apply_ibm_forcing(v, dt)

        use eulerian_ibm_mod  , only : eulerian_forcing_velocity   => forcing_velocity
        use lagrangian_ibm_mod, only : lagrangian_forcing_velocity => forcing_velocity

        ! In/Out variables
        type(vector), intent(inout) :: v   !< Velocity field to be forced
        real(dp)    , intent(in   ) :: dt  !< timestep

        ! Evaluate the forcing due to eulerian solids
        if (allocated(Eulerian_Solid_list)) then
            call eulerian_forcing_velocity(v, eulerian_Solid_list, Fe, dt)
        endif

        if (allocated(Lagrangian_Solid_list)) then
            call lagrangian_forcing_velocity(v, lagrangian_Solid_list, dt)
        endif

    end subroutine apply_ibm_forcing
    !==============================================================================================

    !==============================================================================================
    subroutine destroy_ibm

        use lagrangian_ibm_mod, only : destroy_lagrangian_ibm => destroy_ibm

        !deallocate(Eulerian_Solid_list)
        call destroy_lagrangian_ibm()
        deallocate(lagrangian_solid_list)

    end subroutine destroy_ibm
    !==============================================================================================

end module
