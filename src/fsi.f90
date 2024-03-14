module fsi_mod

    implicit none

contains

    !===============================================================================================
    subroutine weak_coupling_solver(comp_grid, step, dt)

        use precision_mod     , only : dp
        use grid_mod
        use ibm_mod           , only : lagrangian_solid_list, eulerian_solid_list, Fe
        use lagrangian_ibm_mod, only : lagrangian_compute_hydrodynamic_loads => compute_hydrodynamic_loads
        use eulerian_ibm_mod  , only : tag_cells, compute_hydrodynamic_loads
        use navier_stokes_mod , only : navier_stokes_solver, v, p, mu, rho, density, g
#ifdef MF
        use mpi
        use global_mod, only : ierror
#endif
        use solid_mod
        use eulerian_solid_mod

        ! In/Out variables
        type(grid), intent(in   ) :: comp_grid
        integer   , intent(in   ) :: step
        real(dp)  , intent(inout) :: dt

        ! Local variables
        integer                        :: b
        class(eulerian_solid), pointer :: Es
#ifdef MF
        integer                        :: ie(3)
        real(dp)                       :: rhof
#endif

        if (allocated(Eulerian_Solid_list)) then
            
            ! First solve solid object motion
            do b = 1,size(eulerian_solid_list)

                ! Evalute hydrodynamic loads on solid object b
                call compute_hydrodynamic_loads(eulerian_solid_list(b)%pS, v, p, mu, rho, g, Fe, density)

                ! Advance solid object b
                call eulerian_solid_list(b)%pS%advance(dt)

                ! Check if the solid object is outside of the domain and in case traslate it
                call Eulerian_Solid_list(b)%pS%check_periodicity()
            end do

            ! Update eularian ibm variables
            call tag_cells(Eulerian_Solid_list)
        
        endif

        ! If solving solid body dynamics with Lagrangian IBM
        if (allocated(lagrangian_solid_list)) then
        
            ! For every solid body:
            do b = 1,size(lagrangian_solid_list)

                ! For deformable solids, velocity and acceleration on forcing elements are 
                ! evaluated by interpolation from mass points. For rigid solids, velocity and 
                ! acceleration of forcing elements are updated withing the rigid_body subroutine.  
                if (lagrangian_solid_list(b)%pS%is_deformable) then
                    call lagrangian_solid_list(b)%pS%update_forcing_elements()
                endif

                ! Compute hydrodynamic loads with fluid fields at step n
                call lagrangian_compute_hydrodynamic_loads(lagrangian_solid_list(b)%pS, v, p, mu, &
                                                                rho, g)
   
                ! Solve the structure dynamics with the computed hydrodynamic loads
#ifdef MF
                ! For multiphase flows find the density of the surrounding fluid
                ie = comp_grid%closest_grid_node(lagrangian_solid_list(b)%pS%center_of_mass%X(1:3), 0)
                if (ie(2) >= comp_grid%lo(2) .and. ie(2) <= comp_grid%hi(2)) then
                    rhof = rho%f(ie(1),ie(2),1)
                else
                    rhof = 0.0_dp
                endif
                call mpi_allreduce(mpi_in_place,rhof,1,mpi_real8,mpi_sum,mpi_comm_world,ierror)
                call lagrangian_solid_list(b)%pS%advance(dt, g, rhof)
#else
                call lagrangian_solid_list(b)%pS%advance(dt, g, density)
#endif
                call check_periodicity(lagrangian_solid_list(b)%pS, comp_grid)            
            end do
        endif

        ! then solve fluid motion
        call navier_stokes_solver(comp_grid, step, dt)

    end subroutine weak_coupling_solver
    !===============================================================================================

    !===============================================================================================
    subroutine check_periodicity(obj, comp_grid)
        
        use lagrangian_solid_mod , only : lagrangian_solid
        use grid_mod             , only : grid
        use precision_mod        , only : dp

        ! In/Out variables
        class(lagrangian_solid), intent(inout) :: obj
        type(grid)             , intent(in   ) :: comp_grid

        ! Local variables
        integer  :: n
        real(dp) :: DL
        logical  :: is_outside

        ! Set the flag to true
        is_outside = .true.

        ! Check if the body is entirely outside
        if (comp_grid%periodic_bc(1)) then
            do n = 1,obj%number_of_mass_points
                if (obj%mass_points(n)%X(1) > (comp_grid%Lx - comp_grid%origin(1)) .or. &
                    obj%mass_points(n)%X(1) < comp_grid%origin(1)) then
                else
                    ! If one point is inside the domain set to false the flag for the body
                    is_outside = .false.
                end if
            end do
        endif

        ! If all the mass points are outside of the domain translate it
        if (is_outside .and. (comp_grid%periodic_bc(1) .eqv. .true.)) then

            ! Check if need to translate to the left or to the right
            DL = merge(-comp_grid%Lx, comp_grid%Lx, obj%mass_points(1)%X(1) > &
                        (comp_grid%Lx - comp_grid%origin(1)))
            ! Traslate mass points
            do n = 1,obj%number_of_mass_points
                obj%mass_points(n)%X(1) = obj%mass_points(n)%X(1) + DL
            end do
            ! Traslate center of mass
            obj%center_of_mass%X(1) = obj%center_of_mass%X(1) + DL
            ! Update solid
            call obj%update()
        end if

    end subroutine check_periodicity
    !===============================================================================================

end module
