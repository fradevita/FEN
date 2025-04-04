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
        use eulerian_solid_mod     , only : eulerian_solid
        use eulerian_ibm_mod       , only : init_eulerian_ibm
        use lagrangian_solid_mod   , only : lagrangian_solid
        use lagrangian_solid_2D_mod, only : lagrangian_solid_2D
        use lagrangian_ibm_mod     , only : init_lagrangian_ibm => init_ibm
        
        ! In/Out variables
        type(grid), intent(in) :: comp_grid

        ! Local variables
        integer :: n, nls, nes

        ! Allocate memory for the forcing field that will be added to the RHS of the
        ! momentum equation.
        call Fe%allocate(comp_grid)

        ! Cycle over solid list to find number of lagrangian and eulerian solid
        nes = 0
        nls = 0
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
        if (nls > 0) allocate(lagrangian_solid_list(nls))
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

    !===============================================================================================
    subroutine print_solid(step)

        use grid_mod
        use scalar_mod
        use global_mod

        ! In/Out variables
        integer, intent(in) :: step

        ! Local variables
        integer               :: i, j, k, b, nb
        real(dp)              :: x, y, z
        real(dp), allocatable :: distance(:)
        type(grid), pointer   :: G
        type(scalar)          :: temp
        character(len=7 )     :: sn
        character(len=18)     :: filename

        if (allocated(eulerian_solid_list)) then
            nb = size(eulerian_solid_list)
            allocate(distance(nb))
            G => eulerian_solid_list(1)%pS%G
            call temp%allocate(G)

            do k = G%lo(3),G%hi(3)
                do j = G%lo(2),G%hi(2)
                    do i = G%lo(1),G%hi(1)

                        x = G%x(i) + stagger(1, 0)*G%delta
                        y = G%y(j) + stagger(2, 0)*G%delta
                        z = G%z(k) + stagger(3, 0)*G%delta

                        ! Compute distance from every solid body
                        do b = 1,nb
                            distance(b) = eulerian_solid_list(b)%pS%distance([x,y,z])
                        end do
                        temp%f(i,j,k) = minval(distance)
                        
                    end do
                end do
            end do
            write(sn,'(I0.7)') step
            filename = 'data/S_'//sn//'.raw'
            call temp%write(filename)
            deallocate(distance)
            call temp%destroy()
        endif

    end subroutine
    !===============================================================================================

    !==============================================================================================
    subroutine destroy_ibm

        use lagrangian_ibm_mod, only : destroy_lagrangian_ibm => destroy_ibm

        !deallocate(Eulerian_Solid_list)
        call destroy_lagrangian_ibm()
        deallocate(lagrangian_solid_list)

    end subroutine destroy_ibm
    !==============================================================================================

end module
