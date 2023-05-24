module lagrangian_ibm_mod

    ! This module contains all procedures for the Lagrangian Immersed Boundary Method
    ! using the Mooving-Least-Square method. The main procedures for the IBM are
    ! the forcing of the velocity field and the computation of the hydrodynamic loads.

    use precision_mod       , only : dp
    use grid_mod
    use vector_mod
    use lagrangian_solid_mod, only : solid

    implicit none

    ! Forces on the Eulerian grid
    type(vector) :: F

    ! Number of cycle of the forcing.
    integer :: Nfstep = 1

    ! distance of the probe from the body
    real(dp) :: h_probe

contains

    !========================================================================================
    subroutine init_ibm(comp_grid)

        ! This subroutine initialize the variables of the module

        type(grid), intent(in) :: comp_grid

        call F%allocate(comp_grid, 1)

        ! By default set the distance of the probe to 1.2 dx
        h_probe = 1.2_dp*comp_grid%delta

    end subroutine init_ibm
    !========================================================================================

    !========================================================================================
    subroutine forcing_velocity(v, solid_array, dt)

        use mpi
        use global_mod          , only : Ndim
        use mls_mod             , only : Ne, interpolate, phi
        use lagrangian_solid_mod, only : edge, solid_pointer

        ! In/Out variables
        real(dp)    , intent(in   ) :: dt
        type(vector), intent(inout) :: v
        type(solid_pointer), intent(in) :: solid_array(:)

        ! Local variables
        integer    :: b, fstep, l, q, sj, si, ii, jj, ie(3)
        real(dp)   :: lm(Ndim), Vl(Ndim+1), Fl, c
        type(edge) :: l_edge

        ! Cycle over the number of solid bodies
        solid_body_cycle: do b = 1,size(solid_array)

            ! The forcing is applyied iteratively Nfstep times
            forcing_cycle: do fstep = 1,Nfstep

                ! Set to zero the force on the eulerian grid
                F%x%f = 0.0_dp
                F%y%f = 0.0_dp
#if DIM==3
                F%z%f = 0.0_dp
#endif

                ! Cycle over all the edges of the solid body b
                edges_cycle: do l = 1,solid_array(b)%pS%number_of_edges

                    ! Select the edge
                    l_edge = solid_array(b)%pS%edges(l)

                    ! Select the lagrangian marker
                    lm = l_edge%C%X

                    ! If the point is outside the domain translate it
                    lm = traslate(lm, v%G)

                    ! **** X ****
                    ! Find the closest Eulerain u node to lm
                    ie = v%G%closest_grid_node([lm(1), lm(2), 0.0_dp], 1)

                    ! Interpolate u on the lagrangian marker lm
                    Vl = interpolate(v%x, lm, ie, 1)
                    
                    ! Compute Lagrangian force on the Lagrangian marker,
                    ! eq (17) of de Tullio & Pascazio JCP 2016.
                    Fl = (l_edge%C%V(1) - Vl(1))/dt

                    ! Compute the scaling factor, eq. 20 of de Tullio & Pascazio JCP 2016.
                    c = l_edge%l/v%G%delta

                    ! Transfer the Lagrangian force to Eulerian force,
                    ! eq (18) of de Tullio & Pascazio JCP 2016.
                    ! Select the rank containing the Eulerian point
                    if (ie(2) >= v%G%lo(2) .and. ie(2) <= v%G%hi(2)) then
                        q = 1
                        do sj = -1,1
                            jj = ie(2) + sj
                            do si = -1,1
                                ii = ie(1) + si
                                F%x%f(ii,jj,1) = F%x%f(ii,jj,1) + c*phi(1,q)*Fl
                                q = q + 1
                            end do
                        end do
                    endif

                    ! **** Y ****
                    ! Find the closest Eulerain v node to lm
                    ie = v%G%closest_grid_node([lm(1), lm(2), 0.0_dp], 2)

                    ! Interpolate v on the lagrangian marker lm
                    Vl = interpolate(v%y, lm, ie, 2)

                    ! Compute Lagrangian force on the Lagrangian marker,
                    ! eq (17) of de Tullio & Pascazio JCP 2016.
                    Fl = (l_edge%C%V(2) - Vl(1))/dt

                    ! Transfer Lagrangian force to Eulerian force,
                    ! eq (18) of de Tullio & Pascazio JCP 2016.
                    if (ie(2) >= v%G%lo(2) .and. ie(2) <= v%G%hi(2)) then
                        q = 1
                        do sj = -1,1
                            jj = ie(2) + sj
                            do si = -1,1
                                ii = ie(1) + si
                                F%y%f(ii,jj,1) = F%y%f(ii,jj,1) + c*phi(1,q)*Fl
                                q = q + 1
                            end do
                        end do
                    end if

                end do edges_cycle

                ! Apply eulerian force to the velocity field due to solid body b
                call F%update_ghost_nodes()
                v%x%f = v%x%f + F%x%f*dt
                v%y%f = v%y%f + F%y%f*dt
                call v%update_ghost_nodes()

            end do forcing_cycle

        end do solid_body_cycle

    end subroutine forcing_velocity
    !========================================================================================

    !========================================================================================
    subroutine compute_hydrodynamic_loads(obj, v, p, mu, rho, g)

        ! Compute hydrodynamic loads on solid obj given by the velocity field v, pressure
        ! field p, of the flow with viscosity mu, density rho and body force g.

        use global_mod          , only : Ndim, tdof
        use scalar_mod
        use lagrangian_solid_mod, only : solid, edge

        ! In/Out variables
        real(dp)    , intent(in   )         :: g(2)
        type(scalar), intent(in   )         :: p, mu, rho
        type(vector), intent(in   )         :: v
        type(solid ), intent(inout), target :: obj

        ! Local variables
        integer                             :: l, n
        real(dp)                            :: Xcm(Ndim), sum_l
        real(dp), dimension(:), allocatable :: tau11, tau12, tau22, pl
        real(dp), dimension(:), allocatable :: tau11m, tau12m, tau22m, plm
        type(edge), pointer                 :: l_edge

        ! Allocate local shear stress and pressure, the size is the number of
        ! lagrangian markers
        allocate(tau11(obj%number_of_edges))
        allocate(tau12(obj%number_of_edges))
        allocate(tau22(obj%number_of_edges))
        allocate(   pl(obj%number_of_edges))

        ! Evaluate stresses on the surface using the probe on the positive direction of the norm
        call compute_stresses(obj, v, p, mu, rho, g, tau11, tau12, tau22, pl, 1)

        ! If the structure is open it is necessary to compute forces also on the negative side of the
        ! normal
        if (obj%is_open) then
            ! Allocate memory
            allocate(tau11m(obj%number_of_edges))
            allocate(tau12m(obj%number_of_edges))
            allocate(tau22m(obj%number_of_edges))
            allocate(   plm(obj%number_of_edges))

            call compute_stresses(obj, v, p, mu, rho, g, tau11m, tau12m, tau22m, plm, -1)

            ! Add tau and p on both sides
            tau11 = tau11 - tau11m
            tau12 = tau12 - tau12m
            tau22 = tau22 - tau22m
            pl = pl - plm
        end if

        ! Get local forces from stresses
        Xcm = obj%center_of_mass%X(1:tdof)
        do l = 1,obj%number_of_edges
            l_edge => obj%edges(l)
            
            ! Shear stress force
            l_edge%C%Fv(1) = (tau11(l)*l_edge%n(1) + tau12(l)*l_edge%n(2))*l_edge%l
            l_edge%C%Fv(2) = (tau22(l)*l_edge%n(2) + tau12(l)*l_edge%n(1))*l_edge%l
            
            ! Pressure force
            l_edge%C%Fp(1) = -pl(l)*l_edge%n(1)*l_edge%l
            l_edge%C%Fp(2) = -pl(l)*l_edge%n(2)*l_edge%l

            ! Total force
            l_edge%C%Fh(1) = l_edge%C%Fv(1) + l_edge%C%Fp(1)
            l_edge%C%Fh(2) = l_edge%C%Fv(2) + l_edge%C%Fp(2)
            
            ! Torque around Xcm
            l_edge%C%Fh(3) = -l_edge%C%Fh(1)*(l_edge%C%X(2) - Xcm(2)) + &
                              l_edge%C%Fh(2)*(l_edge%C%X(1) - Xcm(1))
        end do

        ! Distribute the forces from the Lagrangian markers to the mass points. The force
        ! on the mass point is given by the weigthed average between the froces on the
        ! lagragian markers using as weight the size of the edges.
        do n = 1,obj%number_of_mass_points
            obj%mass_points(n)%Fh = 0.0_dp
            sum_l = 0.0_dp
            do l = 1,obj%mass_points(n)%number_of_edges
                obj%mass_points(n)%Fh = obj%mass_points(n)%Fh + &
                obj%edges(obj%mass_points(n)%edges_index(l))%C%Fh(1:2)* &
                obj%edges(obj%mass_points(n)%edges_index(l))%l
                sum_l = sum_l + obj%edges(obj%mass_points(n)%edges_index(l))%l
            end do
            obj%mass_points(n)%Fh = obj%mass_points(n)%Fh/sum_l
        end do

        ! Free memory for local shear stress and pressure
        deallocate(tau11, tau12, tau22, pl)
        if (obj%is_open) deallocate(tau11m, tau12m, tau22m, plm)

    end subroutine compute_hydrodynamic_loads
    !==============================================================================================

    !==============================================================================================
    subroutine compute_stresses(obj, v, p, mu, rho, g, tau11, tau12, tau22, pl, probe_sign)

        ! Evaluate stresses on probe location. This additional function is defined to avoid
        ! duplication between positive and negative norm direction 

        use mpi
        use global_mod          , only : Ndim, ierror
        use scalar_mod
        use lagrangian_solid_mod, only : edge
        use mls_mod             , only : interpolate

        ! In/Out variables
        type(solid) , intent(in   ) :: obj        
        real(dp)    , intent(in   ) :: g(2)
        type(scalar), intent(in   ) :: p, mu, rho
        type(vector), intent(in   ) :: v
        real(dp)    , intent(inout) :: tau11(:)
        real(dp)    , intent(inout) :: tau12(:)
        real(dp)    , intent(inout) :: tau22(:)
        real(dp)    , intent(inout) :: pl(:)
        integer     , intent(in   ) :: probe_sign

        ! Local variables
        integer    :: l, ie(3)
        real(dp)   :: X_probe(Ndim), Ul(Ndim+1), Vl(Ndim+1), ppl(Ndim+1)
        type(edge) :: l_edge
        
        ! Cycle over all the edges (lagrangian markers) of the solid
        do l = 1,obj%number_of_edges

            ! Select the local edge
            l_edge = obj%edges(l)

            ! Set local stress and pressure to zero
            tau11(l) = 0.0_dp
            tau12(l) = 0.0_dp
            tau22(l) = 0.0_dp
            pl(l) = 0.0_dp

            ! Generate the probe in the given direction of the normal
            X_probe = l_edge%C%X + l_edge%n*h_probe*probe_sign

            ! Check periodicity and in case translate it
            X_probe = traslate(X_probe, v%G)

            ! **** U ****
            ! Find the closest Eulerian u node to the probe
            ie = v%G%closest_grid_node([X_probe(1), X_probe(2), 0.0_dp], 1)

            ! Interpolate u and its derivatives on the probe
            Ul = interpolate(v%x, X_probe, ie, 1)

            ! **** V ****
            ! Find the closest Eulerian v node to the probe
            ie = v%G%closest_grid_node([X_probe(1), X_probe(2), 0.0_dp], 2)

            ! Interpolate v and its derivatives on the probe
            Vl = interpolate(v%y, X_probe, ie, 2)

            ! **** P ****
            ! Find the closest Eulerian p node
            ie = v%G%closest_grid_node([X_probe(1), X_probe(2), 0.0_dp], 0)

            ! Interpolate p and its derivatives on the probe
            ppl = interpolate(p, X_probe, ie, 0)

            ! Compute pressure on the Lagrangian marker
            ! eq. 46 of de Tullio and Pascazio JCP 2016.
            if (ie(2) >= v%G%lo(2) .and. ie(2) <= v%G%hi(2)) then
                pl(l) = ppl(1) + h_probe*((l_edge%C%A(1) - g(1))*l_edge%n(1)  + &
                                          (l_edge%C%A(2) - g(2))*l_edge%n(2)) * &
                                           rho%f(ie(1),ie(2),1)*probe_sign

                ! Viscous stresses
                tau11(l) = 2.0_dp*mu%f(ie(1),ie(2),1)*Ul(2)
                tau12(l) = mu%f(ie(1),ie(2),1)*(Ul(3) + Vl(2))
                tau22(l) = 2.0_dp*mu%f(ie(1),ie(2),1)*Vl(3)

            else

                tau11(l) = 0.0_dp
                tau12(l) = 0.0_dp
                tau22(l) = 0.0_dp
                pl(l) = 0.0_dp

            endif

        end do

        ! MPI communication for the integration
        call mpi_allreduce(mpi_in_place,tau11,obj%number_of_edges,mpi_real8,mpi_sum,mpi_comm_world,ierror)
        call mpi_allreduce(mpi_in_place,tau12,obj%number_of_edges,mpi_real8,mpi_sum,mpi_comm_world,ierror)
        call mpi_allreduce(mpi_in_place,tau22,obj%number_of_edges,mpi_real8,mpi_sum,mpi_comm_world,ierror)
        call mpi_allreduce(mpi_in_place,   pl,obj%number_of_edges,mpi_real8,mpi_sum,mpi_comm_world,ierror)

    end subroutine
    !==============================================================================================

    !==============================================================================================
    subroutine advance_structure(obj, step, dt, g, density, comp_grid)

        ! This subroutine perform one timestep of the structural solver
        use global_mod

        ! In/Out varialbes
        integer    , intent(in   ) :: step
        real(dp)   , intent(in   ) :: dt, g(:), density
        type(solid), intent(inout) :: obj        
        type(grid) , intent(in   ) :: comp_grid

        ! Local variables
        integer  :: n, substep
        real(dp) :: bf(dofs), Vnp1(dofs), Xnp1(dofs)

        ! If the solid is deformable
        if (obj%is_deformable) then

            ! For deformable structures solve the velocity verlet algorithm 
            ! for every mass point.

            ! Setup the external forces vector as sum of gravity forces and hydrodynamic forces
            do n = 1,obj%number_of_mass_points
                obj%mass_points(n)%Fe = obj%mass_points(n)%Fh + &
                                        obj%mass_points(n)%m*g(1:tdof)
            end do

            if (obj%is_open) then 
                ! Must account directly for buoyancy since the volume is zero
                do n = 1,obj%number_of_mass_points
                obj%mass_points(n)%Fe = obj%mass_points(n)%Fe - &
                                        density*obj%Vol*g(1:tdof)/obj%number_of_mass_points
                end do
            endif

            ! Solve deformation due to internal forces nsubstesp per timestep
            do substep = 1,obj%nsubsteps
                ! Velocity Verlet
                call obj%velocity_verlet(dt/real(obj%nsubsteps, dp))
            end do

        else

            ! For rigid solid body solve Newton's second law for the center of mass

            ! First integrate the local forces to get center of mass forces
            call obj%integrate_hydrodynamic_forces()

            ! Add body forces
            bf = 0.0_dp
            if (obj%is_open) then
                ! For open solid buoyancy must be added explicitly
                bf = obj%Vol*(obj%rho - density)*g
            else
                bf = obj%M*g
            endif

            ! Find new values
            Vnp1 = obj%center_of_mass%V + dt*(obj%center_of_mass%Fh + bf)/obj%M
            Xnp1 = obj%center_of_mass%X + dt*0.5_dp*(Vnp1 + obj%center_of_mass%V)

            ! Save traslation and rotation
            obj%tra = Xnp1(1:tdof) - obj%center_of_mass%X(1:tdof)
            obj%rot = Xnp1(tdof+1:tdof+1:rdof) - obj%center_of_mass%X(tdof+1:tdof+1:rdof)

            ! Update solution
            obj%center_of_mass%V = Vnp1
            obj%center_of_mass%X = Xnp1
            obj%center_of_mass%A = (obj%center_of_mass%Fh + bf)/obj%M

            ! Moves all lagrangian points accordingly
            call obj%rigid_body_motion()
        endif

        ! Check if the body is outside of the domain and in case traslate it
        call check_periodicity(obj, comp_grid)

    end subroutine advance_structure
    !========================================================================================

    !========================================================================================
    function traslate(X, comp_grid) result(X1)

        use global_mod , only : Ndim

        ! In/out variables
        type(grid), intent(in) :: comp_grid
        real(dp)  , intent(in) :: X(Ndim)
        real(dp)               :: X1(Ndim)

        X1 = X

        if (comp_grid%periodic_bc(1) .eqv. .true.) then
            if (X(1) >= comp_grid%Lx) then
                X1(1) = X(1) - comp_grid%Lx
            elseif (X(1) <= comp_grid%origin(1)) then
                X1(1) = X(1) + comp_grid%Lx
            endif
        endif

    end function traslate
    !========================================================================================

    !========================================================================================
    subroutine check_periodicity(obj, comp_grid)

        ! In/Out variables
        type(solid), intent(inout) :: obj
        type(grid) , intent(in   ) :: comp_grid

        ! Local variables
        integer  :: n
        real(dp) :: DL

        ! Set the flag to true
        obj%is_out = .true.

        ! Check if the body is entirely outside
        if (comp_grid%periodic_bc(1) .eqv. .true.) then
            do n = 1,obj%number_of_mass_points
                if (obj%mass_points(n)%X(1) > (comp_grid%Lx - comp_grid%origin(1)) .or. &
                    obj%mass_points(n)%X(1) < comp_grid%origin(1)) then
                else
                    ! If one point is inside the domain set to false the flag for the body
                    obj%is_out = .false.
                end if
            end do
        endif

        ! If all the mass points are outside of the domain translate it
        if (obj%is_out .and. comp_grid%periodic_bc(1) .eqv. .true.) then

            ! Check if need to translate to the left or to the right
            DL = merge(-comp_grid%Lx, comp_grid%Lx, obj%mass_points(1)%X(1) > &
                (comp_grid%Lx - comp_grid%origin(1)))
            do n = 1,obj%number_of_mass_points
                obj%mass_points(n)%X(1) = obj%mass_points(n)%X(1) + DL
            end do
            obj%center_of_mass%X(1) = obj%center_of_mass%X(1) + DL
            obj%is_out = .false.

            call obj%update_lagrangian_markers()
        end if

    end subroutine check_periodicity
    !========================================================================================

    !========================================================================================
    subroutine destroy_ibm

        ! Free the allocated memory
        call F%destroy()

    end subroutine destroy_ibm
    !========================================================================================

end module