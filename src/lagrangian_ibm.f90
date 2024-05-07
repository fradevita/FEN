module lagrangian_ibm_mod

    ! This module contains all procedures for the Lagrangian Immersed Boundary 
    ! Method using the Mooving-Least-Square method. The main procedures for the 
    ! IBM are the forcing of the velocity field and the computation of the 
    ! hydrodynamic loads.

    use precision_mod       , only : dp
    use grid_mod
    use vector_mod
    use lagrangian_solid_mod, only : lagrangian_solid

    implicit none

    ! Forces on the Eulerian grid
    type(vector) :: F

    ! Number of cycle of the forcing (by default 1).
    integer :: Nfstep = 1

    ! distance of the probe from the body
    real(dp) :: h_probe

contains

    !===============================================================================================
    subroutine init_ibm(comp_grid)

        ! Objective: initialize the variables of the module

        type(grid), intent(in) :: comp_grid !< Eulerian grid

        ! Allocate memory for eulerian forces.
        call F%allocate(comp_grid, 1)

        ! By default set the distance of the probe to 1.2 dx
        h_probe = 1.0_dp*comp_grid%delta

    end subroutine init_ibm
    !===============================================================================================

    !===============================================================================================
    subroutine forcing_velocity(v, solid_array, dt)

        ! Objective: compute the eulerian forcing to impose no-slip and no-penetration
        !            boundary condition on solid boundary.

        use global_mod             , only : Ndim
        use mls_mod                , only : Ne, interpolate, phi
        use lagrangian_solid_mod   , only : lagrangian_solid_pointer, forcing_element
        use lagrangian_solid_1D_mod, only : lagrangian_solid_1D
        use lagrangian_solid_2D_mod, only : lagrangian_solid_2D

        ! In/Out variables
        real(dp)                      , intent(in   ) :: dt
        type(vector)                  , intent(inout) :: v
        type(lagrangian_solid_pointer), intent(in   ) :: solid_array(:)

        ! Local variables
        integer                        :: b, fstep, n, q, si, sj, ii, jj, kk, ie(3)
        real(dp)                       :: Xl(3), Vl(Ndim+1), Fl, c, factor
        type(forcing_element), pointer :: forcElem
#if DIM==3
        integer                       :: sk
#endif

        ! Cycle over the number of solid bodies
        solid_body_cycle: do b = 1,size(solid_array)

            select type(var => solid_array(b)%pS)
            class is(lagrangian_solid_1D)
                factor = v%G%delta
            class is(lagrangian_solid_2D)
                factor = 1.0_dp
            end select

            ! The forcing is applyied iteratively Nfstep times
            forcing_cycle: do fstep = 1,Nfstep

                ! Set to zero the force on the Eulerian grid
                call F%setToValue(0.0_dp)

                ! Cycle over forcing points
                marker_cycle: do n = 1,size(solid_array(b)%pS%forcing_elements) 

                    ! Select the local forcing element
                    forcElem => solid_array(b)%pS%forcing_elements(n)
#if DIM==3
                    Xl = forcElem%C%X
#else
                    Xl(1:2) = forcElem%C%X
                    Xl(3) = 0.0_dp
#endif
                    ! Compute the local scaling factor, eq. 20 of de Tullio & Pascazio 
                    ! JCP 2016.
                    c = factor*forcElem%A/v%G%delta**2

                    ! If the point is outside the domain traslate it
                    call traslate(Xl, v%G)

                    !**** U ************************************************************************
                    ! Find the closest Eulerain x-face to Xl
                    ie = v%G%closest_grid_node(Xl, 1)

                    ! Interpolate u on the lagrangian marker Xl
                    Vl = interpolate(v%x, Xl, ie, 1)
                    
                    ! Compute Lagrangian force on the Lagrangian marker,
                    ! eq (17) of de Tullio & Pascazio JCP 2016.
                    Fl = (forcElem%C%V(1) - Vl(1))/dt

                    ! Transfer the Lagrangian force to Eulerian force,
                    ! eq (18) of de Tullio & Pascazio JCP 2016.
                    ! Select the rank containing the Eulerian point
                    if (ie(2) >= v%G%lo(2) .and. ie(2) <= v%G%hi(2) .and. &
                        ie(3) >= v%G%lo(3) .and. ie(3) <= v%G%hi(3)) then
                        q = 1
                        kk = 1
#if DIM==3
                        do sk = -1,1
                            kk = ie(3) + sk
#endif
                            do sj = -1,1
                                jj = ie(2) + sj
                                do si = -1,1
                                    ii = ie(1) + si
                                    F%x%f(ii,jj,kk) = F%x%f(ii,jj,kk) + c*phi(1,q)*Fl
                                    q = q + 1
                                end do
                            end do
#if DIM==3
                        end do
#endif
                    endif
                    !*******************************************************************************

                    !**** V ************************************************************************
                    ! Find the closest Eulerain y-face to Xl
                    ie = v%G%closest_grid_node(Xl, 2)

                    ! Interpolate v on the lagrangian marker lm
                    Vl = interpolate(v%y, Xl, ie, 2)

                    ! Compute Lagrangian force on the Lagrangian marker,
                    ! eq (17) of de Tullio & Pascazio JCP 2016.
                    Fl = (forcElem%C%V(2) - Vl(1))/dt

                    ! Transfer Lagrangian force to Eulerian force,
                    ! eq (18) of de Tullio & Pascazio JCP 2016.
                    if (ie(2) >= v%G%lo(2) .and. ie(2) <= v%G%hi(2) .and. &
                        ie(3) >= v%G%lo(3) .and. ie(3) <= v%G%hi(3)) then
                        q = 1
                        kk = 1
#if DIM==3
                        do sk = -1,1
                            kk = ie(3) + sk
#endif
                            do sj = -1,1
                                jj = ie(2) + sj
                                do si = -1,1
                                    ii = ie(1) + si
                                    F%y%f(ii,jj,kk) = F%y%f(ii,jj,kk) + c*phi(1,q)*Fl
                                    q = q + 1
                                end do
                            end do
#if DIM==3
                        end do
#endif
                    end if

#if DIM==3
                    !**** W ************************************************************************
                    ! Find the closest Eulerain z-face to Xl
                    ie = v%G%closest_grid_node(Xl, 3)

                    ! Interpolate v on the lagrangian marker lm
                    Vl = interpolate(v%z, Xl, ie, 3)

                    ! Compute Lagrangian force on the Lagrangian marker,
                    ! eq (17) of de Tullio & Pascazio JCP 2016.
                    Fl = (forcElem%C%V(3) - Vl(1))/dt

                    ! Transfer Lagrangian force to Eulerian force,
                    ! eq (18) of de Tullio & Pascazio JCP 2016.
                    if (ie(2) >= v%G%lo(2) .and. ie(2) <= v%G%hi(2) .and. &
                        ie(3) >= v%G%lo(3) .and. ie(3) <= v%G%hi(3)) then
                        q = 1
                        do sk = -1,1
                            kk = ie(3) + sk
                            do sj = -1,1
                                jj = ie(2) + sj
                                do si = -1,1
                                    ii = ie(1) + si
                                    F%z%f(ii,jj,kk) = F%z%f(ii,jj,kk) + c*phi(1,q)*Fl
                                    q = q + 1
                                end do
                            end do
                        end do
                    end if
#endif

                end do marker_cycle

                ! Apply eulerian force to the velocity field due to solid body b
                v%x%f = v%x%f + F%x%f*dt
                v%y%f = v%y%f + F%y%f*dt
#if DIM==3
                v%z%f = v%z%f + F%z%f*dt
#endif
                call v%update_ghost_nodes()

            end do forcing_cycle

        end do solid_body_cycle

    end subroutine forcing_velocity
    !===============================================================================================

    !===============================================================================================
    subroutine compute_hydrodynamic_loads(obj, v, p, mu, rho, g)

        ! Compute hydrodynamic loads on solid obj given by the velocity field v, 
        ! pressure field p, of the flow with viscosity mu, density rho and 
        ! body force g.

        use global_mod          , only : Ndim, tdof, rdof, dofs
        use euclidean_mod       , only : crossProduct
        use scalar_mod
        use lagrangian_solid_mod, only : forcing_element

        ! In/Out variables
        real(dp)               , intent(in   )         :: g(:)
        type(scalar)           , intent(in   )         :: p, mu, rho
        type(vector)           , intent(in   )         :: v
        class(lagrangian_solid), intent(inout), target :: obj

        ! Local variables
        integer                             :: l, Nfe
        real(dp)                            :: Xcm(Ndim), d(tdof), Fh_x_d(3)
        type(forcing_element)  , pointer    :: forcElem
        real(dp), dimension(:), allocatable :: tau11, tau12, tau22, pl
        real(dp), dimension(:), allocatable :: tau11m, tau12m, tau22m, plm
        real(dp), dimension(:), allocatable :: tau13, tau23, tau33
        real(dp), dimension(:), allocatable :: tau13m, tau23m, tau33m

        ! Number of forcing Elements
        Nfe = size(obj%forcing_elements)

        ! Allocate local shear stress and pressure, the size is the number of
        ! forcing elements
        allocate(tau11(Nfe))
        allocate(tau12(Nfe))
        allocate(tau22(Nfe))
        allocate(   pl(Nfe))
#if DIM==3
        allocate(tau13(Nfe))
        allocate(tau23(Nfe))
        allocate(tau33(Nfe))
#endif

        ! Evaluate stresses on the surface using the probe on the positive direction of the norm
        call compute_stresses(obj, v, p, mu, rho, g, pl, tau11, tau12, tau22, &
                                tau13, tau23, tau33, 1)

        ! If the structure is open it is necessary to compute forces also on the negative side of
        ! the normal
        if (obj%is_open) then
            ! Allocate memory
            allocate(tau11m(Nfe))
            allocate(tau12m(Nfe))
            allocate(tau22m(Nfe))
            allocate(   plm(Nfe))
#if DIM==3
            allocate(tau13m(Nfe))
            allocate(tau23m(Nfe))
            allocate(tau33m(Nfe))
#endif

            call compute_stresses(obj, v, p, mu, rho, g, plm, tau11m, tau12m, tau22m, &
                                    tau13m, tau23m, tau33m, -1)

            ! Add tau and p on both sides
            pl = pl - plm
            tau11 = tau11 - tau11m
            tau12 = tau12 - tau12m
            tau22 = tau22 - tau22m
#if DIM==3
            tau13 = tau13 - tau13m
            tau23 = tau23 - tau23m
            tau33 = tau33 - tau33m
#endif
        end if

        ! Get local forces from stresses
        Xcm = obj%center_of_mass%X(1:tdof)
        do l = 1,Nfe
            forcElem => obj%forcing_elements(l)
            
            ! Shear stress forces
            forcElem%C%Fv(1) = (tau11(l)*forcElem%n(1) + tau12(l)*forcElem%n(2))*forcElem%A
            forcElem%C%Fv(2) = (tau12(l)*forcElem%n(1) + tau22(l)*forcElem%n(2))*forcElem%A

            ! Pressure force
            forcElem%C%Fp(1) = -pl(l)*forcElem%n(1)*forcElem%A
            forcElem%C%Fp(2) = -pl(l)*forcElem%n(2)*forcElem%A
#if DIM==3
            forcElem%C%Fv(1) = forcElem%C%Fv(1) + tau13(l)*forcElem%n(3)*forcElem%A
            forcElem%C%Fv(2) = forcElem%C%Fv(2) + tau23(l)*forcElem%n(3)*forcElem%A
            forcElem%C%Fv(3) = (tau13(l)*forcElem%n(1) + tau23(l)*forcElem%n(2) + &
                                tau33(l)*forcElem%n(3))*forcElem%A
            forcElem%C%Fp(3) = -pl(l)*forcElem%n(3)*forcElem%A
#endif
            ! Total hydrodynamic forces
            forcElem%C%Fh(1:tdof) = forcElem%C%Fv(1:tdof) + forcElem%C%Fp(1:tdof)

            ! Torque around Xcm
            d = forcElem%C%X(1:tdof) - Xcm
#if DIM==3
            Fh_x_d = crossProduct(d, forcElem%C%Fh)
            forcElem%C%Fh(tdof+1:dofs) = Fh_x_d
#else
            Fh_x_d = crossProduct([d(1), d(2), 0.0_dp],[forcElem%C%Fh(1), forcElem%C%Fh(2), 0.0_dp])
            forcElem%C%Fh(3) = Fh_x_d(3)
#endif
        end do

        ! Transfer forces from lagrangian marker to mass point
        if (obj%is_deformable) call obj%interpolate_from_forcing_element_to_mass_point()

        ! Free memory for local shear stress and pressure
        deallocate(tau11, tau12, tau22, pl)
        if (obj%is_open) deallocate(tau11m, tau12m, tau22m, plm)
#if DIM==3
        deallocate(tau13, tau23, tau33)
        if (obj%is_open) deallocate(tau13m, tau23m, tau33m)
#endif

    end subroutine compute_hydrodynamic_loads
    !===============================================================================================

    !===============================================================================================
    subroutine compute_stresses(obj, v, p, mu, rho, g, pl, tau11, tau12, tau22, &
                                    tau13, tau23, tau33, probe_sign)

        ! Objective: evaluate stresses on probe location.

        use mpi
        use global_mod          , only : Ndim, ierror, myrank
        use scalar_mod
        use lagrangian_solid_mod
        use mls_mod             , only : interpolate

        ! In/Out variables
        class(lagrangian_solid), intent(in   ), target :: obj        
        real(dp)               , intent(in   )         :: g(:)
        type(scalar)           , intent(in   )         :: p, mu, rho
        type(vector)           , intent(in   )         :: v
        real(dp)               , intent(inout)         :: tau11(:)
        real(dp)               , intent(inout)         :: tau12(:)
        real(dp)               , intent(inout)         :: tau22(:)
        real(dp)               , intent(inout)         :: tau13(:)
        real(dp)               , intent(inout)         :: tau23(:)
        real(dp)               , intent(inout)         :: tau33(:)
        real(dp)               , intent(inout)         :: pl(:)
        integer                , intent(in   )         :: probe_sign

        ! Local variables
        integer                       :: l, ie(3), Nfe
        real(dp)                      :: X_probe(3), Ul(Ndim+1), Vl(Ndim+1), ppl(Ndim+1)
        real(dp)                      :: dpdn, visc
#if DIM==3
        real(dp)                      :: Wl(Ndim+1)
#endif
        type(forcing_element), pointer :: forcElem
      
        ! Number of forcin elements
        Nfe = size(obj%forcing_elements)

        ! Cycle over all the forcing elements of the solid
        do l = 1,Nfe

            ! Select the local edge
            forcElem => obj%forcing_elements(l)

            ! Set local stress and pressure to zero
               pl(l) = 0.0_dp
            tau11(l) = 0.0_dp
            tau12(l) = 0.0_dp
            tau22(l) = 0.0_dp
#if DIM==3  
            tau13(l) = 0.0_dp
            tau23(l) = 0.0_dp
            tau33(l) = 0.0_dp
#endif

            ! Generate the probe in the given direction of the normal
            X_probe(1:Ndim) = forcElem%C%X + forcElem%n*h_probe*probe_sign 

            ! Check periodicity and in case translate it
            call traslate(X_probe, v%G)

            ! **** U *******************************************************************************
            ! Find the closest Eulerian u node to the probe
            ie = v%G%closest_grid_node(X_probe, 1)

            ! Interpolate u and its derivatives on the probe
            Ul = interpolate(v%x, X_probe, ie, 1)

            ! **** V *******************************************************************************
            ! Find the closest Eulerian v node to the probe
            ie = v%G%closest_grid_node(X_probe, 2)

            ! Interpolate v and its derivatives on the probe
            Vl = interpolate(v%y, X_probe, ie, 2)

#if DIM==3
            ! **** W *******************************************************************************
            ! Find the closest Eulerian v node to the probe
            ie = v%G%closest_grid_node(X_probe, 3)

            ! Interpolate v and its derivatives on the probe
            Wl = interpolate(v%z, X_probe, ie, 3)
#endif

            ! **** P *******************************************************************************
            ! Find the closest Eulerian p node
            ie = v%G%closest_grid_node(X_probe, 0)

            ! Interpolate p and its derivatives on the probe
            ppl = interpolate(p, X_probe, ie, 0)
            
            ! Compute pressure on the Lagrangian marker
            ! eq. 46 of de Tullio and Pascazio JCP 2016.
            if ((ie(2) >= v%G%lo(2) .and. ie(2) <= v%G%hi(2)) .and. &
                (ie(3) >= v%G%lo(3) .and. ie(3) <= v%G%hi(3))) then

                ! Normal pressure gradient on the lagrangian marker
                dpdn = -rho%f(ie(1),ie(2),ie(3))*((forcElem%C%A(1) - g(1))*forcElem%n(1)  + &
                                                  (forcElem%C%A(2) - g(2))*forcElem%n(2))*probe_sign
#if DIM==3
                dpdn = dpdn -rho%f(ie(1),ie(2),ie(3))*(forcElem%C%A(3) - g(3))*&
                                forcElem%n(3)*probe_sign
#endif
                pl(l) = ppl(1) - h_probe*dpdn

                ! Viscous stresses: tau_ij = 2 mu D_ij
                visc = mu%f(ie(1),ie(2),ie(3))
                tau11(l) = 2.0_dp*visc*Ul(2)
                tau12(l) = visc*(Ul(3) + Vl(2))
                tau22(l) = 2.0_dp*visc*Vl(3)
#if DIM==3
                tau13(l) = visc*(Ul(4) + Wl(2))
                tau23(l) = visc*(Vl(4) + Wl(3))
                tau33(l) = 2.0_dp*visc*Wl(4)
#endif
            else
                   pl(l) = 0.0_dp
                tau11(l) = 0.0_dp
                tau12(l) = 0.0_dp
                tau22(l) = 0.0_dp
#if DIM==3
                tau13(l) = 0.0_dp
                tau23(l) = 0.0_dp
                tau33(l) = 0.0_dp
#endif
            endif

        end do

        ! MPI communication for the integration
        call mpi_allreduce(mpi_in_place, tau11, Nfe, mpi_real8, mpi_sum, mpi_comm_world, ierror)
        call mpi_allreduce(mpi_in_place, tau12, Nfe, mpi_real8, mpi_sum, mpi_comm_world, ierror)
        call mpi_allreduce(mpi_in_place, tau22, Nfe, mpi_real8, mpi_sum, mpi_comm_world, ierror)
        call mpi_allreduce(mpi_in_place,    pl, Nfe, mpi_real8, mpi_sum, mpi_comm_world, ierror)
#if DIM==3
        call mpi_allreduce(mpi_in_place, tau13, Nfe, mpi_real8, mpi_sum, mpi_comm_world, ierror)
        call mpi_allreduce(mpi_in_place, tau23, Nfe, mpi_real8, mpi_sum, mpi_comm_world, ierror)
        call mpi_allreduce(mpi_in_place, tau33, Nfe, mpi_real8, mpi_sum, mpi_comm_world, ierror)
#endif

    end subroutine
    !===============================================================================================

    !===============================================================================================
    subroutine traslate(X, comp_grid)

        use global_mod , only : Ndim

        ! In/out variables
        type(grid), intent(in   ) :: comp_grid
        real(dp)  , intent(inout) :: X(3)

        if (comp_grid%periodic_bc(1)) then
            if (X(1) > comp_grid%Lx) then
                X(1) = X(1) - comp_grid%Lx
            elseif (X(1) < comp_grid%origin(1)) then
                X(1) = X(1) + comp_grid%Lx
            endif
        endif

    end subroutine traslate
    !===============================================================================================

    !===============================================================================================
    subroutine destroy_ibm

        ! Free the allocated memory
        call F%destroy()

    end subroutine destroy_ibm
    !===============================================================================================

end module
