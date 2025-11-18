module lagrangian_ibm_mod

    ! This module contains all procedures for the Lagrangian Immersed Boundary 
    ! Method using the Mooving-Least-Square method. The main procedures for the 
    ! IBM are the forcing of the velocity field and the computation of the 
    ! hydrodynamic loads.

    use precision_mod       , only : dp
    use grid_mod            , only : grid
    use vector_mod          , only : vector
    use lagrangian_solid_mod, only : solid => lagrangian_solid, forcing_element

    implicit none

    ! Forces on the Eulerian grid
    type(vector) :: F

    ! Number of cycle for the forcing (by default 1).
    integer :: Nfstep = 1

    ! two probes location for stress evaluation
    logical  :: mls_derivatives = .false.
    integer  :: nprobes = 2
    real(dp) :: hp1, hp2

    interface
        subroutine stress_model(obj, v, p, mu, rho, g)
            use precision_mod, only : dp
            use scalar_mod   , only : scalar
            use vector_mod   , only : vector
            import solid
            ! In/Out variables
            real(dp)    , intent(in   )         :: g(:)
            type(scalar), intent(in   )         :: p, mu, rho
            type(vector), intent(in   )         :: v
            class(solid), intent(inout), target :: obj
        end subroutine stress_model
    end interface

    procedure(stress_model), pointer :: compute_hydrodynamic_loads => hydrodynamic_loads_from_probes

contains

    !===============================================================================================
    subroutine init_ibm(comp_grid)

        ! Objective: initialize the variables of the module

        type(grid), intent(in) :: comp_grid !< Eulerian grid

        ! Allocate memory for eulerian forces.
        call F%allocate(comp_grid, 1)

        ! By default set the distance for the two probes
        hp1 = 2.0_dp*comp_grid%delta
        hp2 = 3.0_dp*comp_grid%delta

    end subroutine init_ibm
    !===============================================================================================

    !===============================================================================================
    subroutine forcing_velocity(v, solid_array, dt)

        ! Objective: compute the eulerian forcing to impose no-slip and 
        !            no-penetration boundary condition on solid boundary.

        use global_mod             , only : Ndim
        use mls_mod                , only : Ne, m, get_phi, alpha_s
        use lagrangian_solid_mod   , only : solid_pointer => lagrangian_solid_pointer
        use lagrangian_solid_1D_mod, only : solid_1D => lagrangian_solid_1D
        use lagrangian_solid_2D_mod, only : solid_2D => lagrangian_solid_2D

        ! In/Out variables
        real(dp)           , intent(in   ) :: dt
        type(vector)       , intent(inout) :: v
        type(solid_pointer), intent(in   ) :: solid_array(:)

        ! Local variables
        integer                        :: b, fstep, n, ie(3), q, si, sj, ii, jj, kk
        real(dp)                       :: delta, factor, Xl(3), c, fk(Ne), xs(Ndim,Ne)
        real(dp)                       :: phi(m,Ne), Ul, Fl, ds(Ndim,Ne)
        type(forcing_element), pointer :: forcElem
#if DIM==3
        integer                        :: sk
#endif

        delta = v%x%G%delta

        ! Cycle over the number of solid bodies
        solid_body_cycle: do b = 1,size(solid_array)

            select type(var => solid_array(b)%pS)
            class is(solid_1D)
                factor = delta
            class is(solid_2D)
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
                    ! Compute the local scaling factor, eq. 20 of de Tullio & Pascazio JCP 2016.
                    c = factor*forcElem%A/delta**2

                    ! If the point is outside the domain traslate it
                    call traslate(Xl, v%G)

                    !**** U ************************************************************************
                    ! Find the closest Eulerain x-face to Xl
                    ie = v%G%closest_grid_node(Xl, 1)

#ifdef MPI
                    ! select the proper rank
                    if ( (ie(2) >= v%x%G%lo(2) .and. ie(2) <= v%x%G%hi(2)) .and. &
                         (ie(3) >= v%x%G%lo(3) .and. ie(3) <= v%x%G%hi(3))) then 
#endif

                        ! Build the array of positions and f values in the support domain
                        q = 1
                        xs = 0.0_dp
                        fk = 0.0_dp
                        kk = 1
#if DIM==3
                        do sk = -1,1
                            kk = ie(3) + sk
#endif
                            do sj = -1,1
                                jj = ie(2) + sj
                                do si = -1,1
                                    ii = ie(1) + si
                                    xs(1,q) = v%x%G%x(ie(1)) + 0.5_dp*delta + si*delta
                                    xs(2,q) = v%x%G%y(ie(2)) + sj*delta
                                    fk(q) = v%x%f(ii,jj,kk)
#if DIM==3
                                    xs(3,q) = v%x%G%z(ie(3)) + sk*delta
#endif
                                    q = q + 1
                                end do
                            end do
#if DIM==3
                        end do
#endif
                        !call build_support_domain(v%x, ie, 1, xs, fk)

                        ! Evaluate shape function
                        ds = alpha_s*delta
                        phi = get_Phi(Xl, xs, ds, .false.)

                        ! Interpolate velocity on the Lagrangian marker, 
                        ! eq (14) of de Tullio & Pascazio JCP 2016.
                        Ul = 0.0_dp
                        do q = 1, Ne
                            Ul = Ul + phi(1,q)*fk(q)
                        end do

                        ! Compute Lagrangian force on the Lagrangian marker,
                        ! eq (17) of de Tullio & Pascazio JCP 2016.
                        Fl = (forcElem%C%V(1) - Ul)/dt

                        ! Transfer the Lagrangian force to Eulerian force,
                        ! eq (18) of de Tullio & Pascazio JCP 2016.
                        ! Select the rank containing the Eulerian point
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
#ifdef MPI
                    endif ! exit rank
#endif
                    !*******************************************************************************

                    !**** V ************************************************************************
                    ! Find the closest Eulerain y-face to Xl
                    ie = v%G%closest_grid_node(Xl, 2)
#ifdef MPI
                    ! select the proper rank
                    if ( (ie(2) >= v%x%G%lo(2) .and. ie(2) <= v%x%G%hi(2)) .and. &
                         (ie(3) >= v%x%G%lo(3) .and. ie(3) <= v%x%G%hi(3))) then 
#endif
                        ! Build the array of positions and f values in the support domain
                        q = 1
                        xs = 0.0_dp
                        fk = 0.0_dp
                        kk = 1
#if DIM==3
                        do sk = -1,1
                            kk = ie(3) + sk
#endif
                            do sj = -1,1
                                jj = ie(2) + sj
                                do si = -1,1
                                    ii = ie(1) + si
                                    xs(1,q) = v%G%x(ie(1)) + si*delta
                                    xs(2,q) = v%G%y(ie(2)) + 0.5_dp*delta + sj*delta
                                    fk(q) = v%y%f(ii,jj,kk)
#if DIM==3
                                    xs(3,q) = v%G%z(ie(3)) + sk*delta
#endif
                                    q = q + 1
                                end do
                            end do
#if DIM==3
                        end do
#endif
                        !call build_support_domain(v%y, ie, 2, xs, fk)

                        ! Evaluate shape function
                        ds = alpha_s*delta
                        phi = get_Phi(Xl, xs, ds, .false.)

                        ! Interpolate velocity on the Lagrangian marker, 
                        ! eq (14) of de Tullio & Pascazio JCP 2016.
                        Ul = 0.0_dp
                        do q = 1, Ne
                            Ul = Ul + phi(1,q)*fk(q)
                        end do

                        ! Compute Lagrangian force on the Lagrangian marker,
                        ! eq (17) of de Tullio & Pascazio JCP 2016.
                        Fl = (forcElem%C%V(2) - Ul)/dt

                        ! Transfer the Lagrangian force to Eulerian force,
                        ! eq (18) of de Tullio & Pascazio JCP 2016.
                        ! Select the rank containing the Eulerian point
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
#ifdef MPI
                    endif ! exit rank
#endif
                    !*******************************************************************************

#if DIM==3
                    !**** W ************************************************************************
                    ! Find the closest Eulerain z-face to Xl
                    ie = v%G%closest_grid_node(Xl, 3)

#ifdef MPI
                    ! select the proper rank
                    if ( (ie(2) >= v%x%G%lo(2) .and. ie(2) <= v%x%G%hi(2)) .and. &
                         (ie(3) >= v%x%G%lo(3) .and. ie(3) <= v%x%G%hi(3))) then 
#endif
                        ! Build the array of positions and f values in the support domain
                        q = 1
                        xs = 0.0_dp
                        fk = 0.0_dp
                        kk = 1

                        do sk = -1,1
                            kk = ie(3) + sk
                            do sj = -1,1
                                jj = ie(2) + sj
                                do si = -1,1
                                    ii = ie(1) + si
                                    xs(1,q) = v%x%G%x(ie(1)) + si*delta
                                    xs(2,q) = v%x%G%y(ie(2)) + sj*delta
                                    fk(q) = v%z%f(ii,jj,kk)
                                    xs(3,q) = v%x%G%z(ie(3)) + 0.5_dp*delta + sk*delta
                                    q = q + 1
                                end do
                            end do
                        end do

                        ! Evaluate shape function
                        ds = alpha_s*delta
                        phi = get_Phi(Xl, xs, ds, .false.)

                        ! Interpolate velocity on the Lagrangian marker, 
                        ! eq (14) of de Tullio & Pascazio JCP 2016.
                        Ul = 0.0_dp
                        do q = 1, Ne
                            Ul = Ul + phi(1,q)*fk(q)
                        end do

                        ! Compute Lagrangian force on the Lagrangian marker,
                        ! eq (17) of de Tullio & Pascazio JCP 2016.
                        Fl = (forcElem%C%V(3) - Ul)/dt

                        ! Transfer the Lagrangian force to Eulerian force,
                        ! eq (18) of de Tullio & Pascazio JCP 2016.
                        ! Select the rank containing the Eulerian point
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
#ifdef MPI
                    endif ! exit rank
#endif
                   
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
    subroutine hydrodynamic_loads_from_probes(obj, v, p, mu, rho, g)

        ! Compute hydrodynamic loads on solid object obj given by the velocity 
        ! field v, pressure field p, of the fluid with viscosity mu, density rho 
        ! and body force g.

        use global_mod          , only : Ndim, tdof, rdof, dofs
        use euclidean_mod       , only : crossProduct
        use scalar_mod          , only : scalar
        use tensor_mod          , only : tensor
        use lagrangian_solid_mod, only : forcing_element

        ! In/Out variables
        real(dp)    , intent(in   )         :: g(:)
        type(scalar), intent(in   )         :: p, mu, rho
        type(vector), intent(in   )         :: v
        class(solid), intent(inout), target :: obj

        ! Local variables
        integer                               :: l, Nfe
        real(dp)                              :: Xcm(Ndim), r(tdof), Fh_x_r(3)
        type(tensor)                          :: tau
        type(forcing_element)   , pointer     :: forcElem
        real(dp), dimension(:)  , allocatable :: pl_np, pl_nm
        real(dp), dimension(:,:), allocatable :: taul_np, taul_nm

        ! Evaluate global stress tensor
        if (mls_derivatives) then
            ! Interpolation is made on velocity
        else
            ! Interpolation is made on shear stress, must evaluate tau
            call tau%allocate(v%G, 1)

            ! Here tau is evaluated on cell vertex
            block
                integer  :: i, j, k, im, ip, jm, jp
                real(dp) :: idelta, dudx, dudy, dvdx, dvdy
#if DIM==3
                integer  :: kp, km
                real(dp) :: dudz, dvdz, dwdx, dwdy, dwdz
#endif
                idelta = 1.0_dp/p%G%delta
                do k = v%G%lo(3),v%G%hi(3)
#if DIM==3
                    kp = k + 1
                    km = k - 1
#endif
                    do j = v%G%lo(2),v%G%hi(2)
                        jp = j + 1
                        jm = j - 1
                        do i = v%G%lo(1),v%G%hi(1)
                            ip = i + 1
                            im = i - 1
                            dudx = (v%x%f( i, j, k) - v%x%f(im, j, k))*idelta
                            dudy = (v%x%f( i,jp, k) + v%x%f(im,jp, k) - &
                                    v%x%f( i,jm, k) - v%x%f(im,jm, k))*idelta*0.25_dp
                            dvdx = (v%y%f(ip, j, k) + v%y%f(ip,jm, k) - &
                                    v%y%f(im, j, k) - v%y%f(im,jm, k))*idelta*0.25_dp
                            dvdy = (v%y%f(i,j,k) - v%y%f(i,jm,k))*idelta

                            tau%x%x%f(i,j,k) = mu%f(i,j,k)*(dudx + dudx)
                            tau%x%y%f(i,j,k) = mu%f(i,j,k)*(dudy + dvdx)
                            tau%y%x%f(i,j,k) = tau%x%y%f(i,j,k)
                            tau%y%y%f(i,j,k) = mu%f(i,j,k)*(dvdy + dvdy)
#if DIM==3
                            dudz = (v%x%f( i, j,kp) + v%x%f(im, j,kp) - &
                                    v%x%f( i, j,km) - v%x%f(im, j,km))*idelta*0.25_dp
                            dvdz = (v%y%f( i, j,kp) + v%y%f( i,jm,kp) - &
                                    v%y%f( i, j,km) - v%y%f( i,jm,km))*idelta*0.25_dp
                            dwdx = (v%z%f(ip, j, k) + v%z%f(ip, j,km) - &
                                    v%z%f(im, j, k) - v%z%f(im, j,km))*idelta*0.25_dp
                            dwdy = (v%z%f( i,jp, k) + v%z%f( i,jp,km) - &
                                    v%z%f( i,jm, k) - v%z%f(i, jm,km))*idelta*0.25_dp
                            dwdz = (v%z%f( i, j, k) - v%z%f( i, j,km))*idelta

                            tau%x%z%f(i,j,k) = mu%f(i,j,k)*(dudz + dwdx)
                            tau%y%z%f(i,j,k) = mu%f(i,j,k)*(dvdz + dwdy)
                            tau%z%x%f(i,j,k) = tau%x%z%f(i,j,k)
                            tau%z%y%f(i,j,k) = tau%y%z%f(i,j,k)
                            tau%z%z%f(i,j,k) = mu%f(i,j,k)*(dwdz + dwdz)
#endif
                        end do
                    end do
                end do
            end block
            call tau%update_ghost_nodes()
        endif

        ! Number of forcing Elements
        Nfe = size(obj%forcing_elements)

        ! Allocate local (l) stress, the size is the number of forcing elements
        allocate(pl_np(Nfe))
        allocate(taul_np(6,Nfe)) ! order is: tauxx, tauxy, tauyy, tauzx, tauzy, tauzz

        ! Evaluate stresses on the surface using the probe on the positive direction of the norm
        if (mls_derivatives) then
            call compute_stresses_DMLS(obj, v, p, rho, mu, g, pl_np, taul_np, 1)
        else
            call compute_stresses(obj, tau, p, rho, g, pl_np, taul_np, 1)
        endif

        ! If the structure is open it is necessary to compute forces also on the negative side of
        ! the normal
        if (obj%is_open) then
            ! Allocate memory
            allocate(pl_nm(Nfe))
            allocate(taul_nm(6, Nfe))

            call compute_stresses(obj, tau, p, rho, g, pl_nm, taul_nm, -1)

            ! Add tau and p on both sides
            pl_np = pl_np - pl_nm
            taul_np = taul_np - taul_nm
        end if

        obj%p_s = pl_np
        obj%tau_s = taul_np

        ! Get local forces from stresses
        Xcm = obj%center_of_mass%X(1:tdof)
        do l = 1,Nfe
            forcElem => obj%forcing_elements(l)

            ! Shear stress forces
            forcElem%C%Fv(1) = (taul_np(1,l)*forcElem%n(1) + taul_np(2,l)*forcElem%n(2))*forcElem%A
            forcElem%C%Fv(2) = (taul_np(2,l)*forcElem%n(1) + taul_np(3,l)*forcElem%n(2))*forcElem%A

            ! Pressure force
            forcElem%C%Fp(1) = -pl_np(l)*forcElem%n(1)*forcElem%A
            forcElem%C%Fp(2) = -pl_np(l)*forcElem%n(2)*forcElem%A
#if DIM==3
            forcElem%C%Fv(1) = forcElem%C%Fv(1) + taul_np(4,l)*forcElem%n(3)*forcElem%A
            forcElem%C%Fv(2) = forcElem%C%Fv(2) + taul_np(5,l)*forcElem%n(3)*forcElem%A
            forcElem%C%Fv(3) = (taul_np(4,l)*forcElem%n(1) + taul_np(5,l)*forcElem%n(2) + &
                                taul_np(6,l)*forcElem%n(3))*forcElem%A
            forcElem%C%Fp(3) = -pl_np(l)*forcElem%n(3)*forcElem%A
#endif
            ! Total hydrodynamic forces
            forcElem%C%Fh(1:tdof) = forcElem%C%Fv(1:tdof) + forcElem%C%Fp(1:tdof)

            ! Torque around Xcm
            r = forcElem%C%X(1:tdof) - Xcm
#if DIM==3
            Fh_x_r = crossProduct(r, forcElem%C%Fh)
            forcElem%C%Fh(tdof+1:dofs) = Fh_x_r
#else
            Fh_x_r = crossProduct([r(1), r(2), 0.0_dp],[forcElem%C%Fh(1), forcElem%C%Fh(2), 0.0_dp])
            forcElem%C%Fh(3) = Fh_x_r(3)
#endif
        end do

        ! If solving for deformable structure, transfer forces from lagrangian marker to mass point
        if (obj%is_deformable) call obj%interpolate_from_forcing_element_to_mass_point()

        ! Free memory for local shear stress and pressure
        deallocate(taul_np, pl_np)
        if (obj%is_open) deallocate(taul_nm, pl_nm)

    end subroutine hydrodynamic_loads_from_probes
    !===============================================================================================

    !===============================================================================================
    subroutine compute_stresses(obj, tau, p, rho, g, pl, taul, probe_sign)

        ! Objective: evaluate stresses on solid surface.

        use mpi
        use global_mod, only : Ndim, ierror, myrank
        use scalar_mod, only : scalar
        use tensor_mod, only : tensor
        use mls_mod   , only : Ne, m, get_phi, alpha_s

        ! In/Out variables
        class(solid), intent(in   ), target :: obj        
        real(dp)    , intent(in   )         :: g(:)
        type(scalar), intent(in   )         :: p, rho
        type(tensor), intent(in   )         :: tau
        real(dp)    , intent(inout)         :: taul(:,:)
        real(dp)    , intent(inout)         :: pl(:)
        integer     , intent(in   )         :: probe_sign

        ! Local variables
        integer                        :: Nfe, l, ie(3), q, ii, jj, kk, si, sj
        real(dp)                       :: delta, a, dpdn, Xp1(3), xs(Ndim,Ne), pk(Ne), tauk(6,Ne)
        real(dp)                       :: phi(m,Ne), pp1, taup1(6), Xp2(3), pp2, taup2(6), ds(Ndim,Ne)
#if DIM==3
        integer                        :: sk
#endif
        type(forcing_element), pointer :: forcElem

        delta = p%G%delta

        ! Number of forcing elements
        Nfe = size(obj%forcing_elements)

        ! Two probe location factor
        a = hp1/hp2

        ! Set stresses to zero
        pl = 0.0_dp
        taul = 0.0_dp

        ! Cycle over all the forcing elements of the solid
        do l = 1,Nfe

            ! Select the local forcing elemeng
            forcElem => obj%forcing_elements(l)

            ! **** First probe *********************************************************************
            pp1 = 0.0_dp
            taup1 = 0.0_dp
            dpdn = 0.0_dp 

            ! Generate the probe in the given normal direction
            Xp1(1:Ndim) = forcElem%C%X + forcElem%n*hp1*probe_sign
            
            ! Check periodicity and in case translate it
            !call traslate(X_probe, v%G)

            ! Find the closest Eulerian cell center to the probe
            ie = p%G%closest_grid_node(Xp1, 0)
#ifdef MPI
            ! select the proper rank
            if ( (ie(2) >= p%G%lo(2) .and. ie(2) <= p%G%hi(2)) .and. &
                 (ie(3) >= p%G%lo(3) .and. ie(3) <= p%G%hi(3))) then
#endif
                dpdn = -rho%f(ie(1),ie(2),ie(3))*((forcElem%C%A(1) - g(1))*forcElem%n(1)  + &
                                                  (forcElem%C%A(2) - g(2))*forcElem%n(2))*probe_sign
#if DIM==3
                dpdn = dpdn - &
                        rho%f(ie(1),ie(2),ie(3))*((forcElem%C%A(3) - g(3))*forcElem%n(3))*probe_sign
#endif

                ! Build the array of positions and field values in the support domain
                q = 1
                xs = 0.0_dp
                pk = 0.0_dp
                tauk = 0.0_dp
                kk = 1
#if DIM==3
                do sk = -1,1
                    kk = ie(3) + sk
#endif
                    do sj = -1,1
                        jj = ie(2) + sj
                        do si = -1,1
                            ii = ie(1) + si
                            xs(1,q)   = p%G%x(ie(1)) + si*delta
                            xs(2,q)   = p%G%y(ie(2)) + sj*delta
                            pk(q)     = p%f(ii,jj,kk)
                            tauk(1,q) = tau%x%x%f(ii,jj,kk)
                            tauk(2,q) = tau%x%y%f(ii,jj,kk)
                            tauk(3,q) = tau%y%y%f(ii,jj,kk)
#if DIM==3
                            xs(3,q)   = p%G%z(ie(3)) + sk*delta
                            tauk(4,q) = tau%z%x%f(ii,jj,kk)
                            tauk(5,q) = tau%z%y%f(ii,jj,kk)
                            tauk(6,q) = tau%z%z%f(ii,jj,kk)
#endif
                            q = q + 1
                        end do
                    end do
#if DIM==3
                end do
#endif

                ! Evaluate shape function for the probe
                ds = alpha_s*delta
                phi = get_Phi(Xp1, xs, ds, .false.)

                ! Interpolate pressure and viscous stresses on probe 1
                do q = 1, Ne
                    pp1      = pp1 + phi(1,q)*pk(q)
                    taup1(1) = taup1(1) + phi(1,q)*tauk(1,q)
                    taup1(2) = taup1(2) + phi(1,q)*tauk(2,q)
                    taup1(3) = taup1(3) + phi(1,q)*tauk(3,q)
#if DIM==3
                    taup1(4) = taup1(4) + phi(1,q)*tauk(4,q)
                    taup1(5) = taup1(5) + phi(1,q)*tauk(5,q)
                    taup1(6) = taup1(6) + phi(1,q)*tauk(6,q)
#endif
                end do
#ifdef MPI
            endif ! exit rank
#endif

            if (nprobes > 1) then
                pp2 = 0.0_dp
                taup2 = 0.0_dp
                ! **** Second probe ****************************************************************
                Xp2(1:Ndim) = forcElem%C%X + forcElem%n*hp2*probe_sign 
                ie = p%G%closest_grid_node(Xp2, 0)
#ifdef MPI
                ! select the proper rank
                if ( (ie(2) >= p%G%lo(2) .and. ie(2) <= p%G%hi(2)) .and. &
                     (ie(3) >= p%G%lo(3) .and. ie(3) <= p%G%hi(3))) then 
#endif
                    ! Build the array of positions and f values in the support domain
                    q = 1
                    xs = 0.0_dp
                    pk = 0.0_dp
                    tauk = 0.0_dp
                    kk = 1
#if DIM==3
                    do sk = -1,1
                        kk = ie(3) + sk
#endif
                        do sj = -1,1
                            jj = ie(2) + sj
                            do si = -1,1
                                ii = ie(1) + si
                                xs(1,q)   = p%G%x(ie(1)) + si*delta
                                xs(2,q)   = p%G%y(ie(2)) + sj*delta
                                pk(q)     = p%f(ii,jj,kk)
                                tauk(1,q) = tau%x%x%f(ii,jj,kk)
                                tauk(2,q) = tau%x%y%f(ii,jj,kk)
                                tauk(3,q) = tau%y%y%f(ii,jj,kk)
#if DIM==3
                                xs(3,q)   = p%G%z(ie(3)) + sk*delta
                                tauk(4,q) = tau%z%x%f(ii,jj,kk)
                                tauk(5,q) = tau%z%y%f(ii,jj,kk)
                                tauk(6,q) = tau%z%z%f(ii,jj,kk)
#endif
                                q = q + 1
                            end do
                        end do
#if DIM==3
                    end do
#endif
    
                    ! Evaluate shape function for probe 1
                    ds = alpha_s*delta
                    phi = get_Phi(Xp2, xs, ds, .false.)
    
                    ! Interpolate pressure and viscous stresses on probe 1
                    do q = 1, Ne
                        pp2      = pp2 + phi(1,q)*pk(q)
                        taup2(1) = taup2(1) + phi(1,q)*tauk(1,q)
                        taup2(2) = taup2(2) + phi(1,q)*tauk(2,q)
                        taup2(3) = taup2(3) + phi(1,q)*tauk(3,q)
#if DIM==3
                        taup2(4) = taup2(4) + phi(1,q)*tauk(4,q)
                        taup2(5) = taup2(5) + phi(1,q)*tauk(5,q)
                        taup2(6) = taup2(6) + phi(1,q)*tauk(6,q)
#endif
                    end do
#ifdef MPI
                endif ! exit rank
#endif

                ! **** Surface values***************************************************************
                ! Extrapolate shear stress from probes on surface
                taul(1:3,l) = (taup1(1:3) - a*taup2(1:3))/(1.0_dp - a)
#if DIM==3
                taul(4:6,l) = (taup1(4:6) - a*taup2(4:6))/(1.0_dp - a)
#endif

                pl(l) = pp1/hp1**2 - pp2/hp2**2 - dpdn*(1.0_dp/hp1 - 1.0_dp/hp2)
                pl(l) = pl(l)/(1.0_dp/hp1**2 - 1.0_dp/hp2**2)
            else
                taul(1:3,l) = taup1(1:3)
                pl(l) = pp1 -hp1*dpdn
            endif
        end do

#ifdef MPI
        call mpi_allreduce(mpi_in_place,   pl,  Nfe, mpi_real8, mpi_sum, mpi_comm_world, ierror)
        call mpi_allreduce(mpi_in_place, taul,  6*Nfe, mpi_real8, mpi_sum, mpi_comm_world, ierror)
#endif
    
    end subroutine
    !===============================================================================================

    !===============================================================================================
    subroutine compute_stresses_DMLS(obj, v, p, rho, mu, g, pl, taul, probe_sign)

        ! Objective: evaluate stresses on solid surface using derivatives of shape function.

        use mpi
        use global_mod, only : Ndim, ierror, myrank
        use scalar_mod, only : scalar
        use tensor_mod, only : tensor
        use mls_mod

        ! In/Out variables
        class(solid), intent(in   ), target :: obj        
        real(dp)    , intent(in   )         :: g(:)
        type(scalar), intent(in   )         :: p, rho, mu
        type(vector), intent(in   )         :: v
        real(dp)    , intent(inout)         :: taul(:,:)
        real(dp)    , intent(inout)         :: pl(:)
        integer     , intent(in   )         :: probe_sign

        ! Local variables
        integer                        :: Nfe, l, ie(3), q, si, sj, sk, ii, jj, kk
        real(dp)                       :: delta, a, pp1, Up1(Ndim+1), Vp1(Ndim+1), Xp1(3), dpdn
        real(dp)                       :: xs(Ndim,Ne), pk(Ne), ds(Ndim,Ne), phi(m,Ne), uk(Ne)
        real(dp)                       :: pp2, Up2(Ndim+1), Vp2(Ndim+1), Xp2(3)
        real(dp)                       :: dUdx, dUdy, dVdx, dVdy 
#if DIM==3
        real(dp)                       :: Wp1(Ndim+1), Wp2(Ndim+1), dUdz, dVdz, dWdx, dWdy, dWdz
#endif
        type(forcing_element), pointer :: forcElem

        delta = p%G%delta

        ! Number of forcing elements
        Nfe = size(obj%forcing_elements)

        ! Two probe location factor
        a = hp1/hp2

        ! Set stresses to zero
        pl = 0.0_dp
        taul = 0.0_dp

        ! Cycle over all the forcing elements of the solid
        do l = 1,Nfe

            ! Select the local forcing elemeng
            forcElem => obj%forcing_elements(l)

            ! **** First probe *********************************************************************
            pp1 = 0.0_dp
            Up1 = 0.0_dp
            Vp1 = 0.0_dp
            dpdn = 0.0_dp
#if DIM==3
            Wp1 = 0.0_dp
#endif
            ! Generate the probes in the given normal direction
            Xp1(1:Ndim) = forcElem%C%X + forcElem%n*hp1*probe_sign 
            
            ! Check periodicity and in case translate it
            !call traslate(X_probe, v%G)

            ! Find the closest Eulerian cell center to probe 1 and interpolate pressure
            ie = p%G%closest_grid_node(Xp1, 0)
#ifdef MPI
            ! select the proper rank
            if ( (ie(2) >= p%G%lo(2) .and. ie(2) <= p%G%hi(2)) .and. &
                 (ie(3) >= p%G%lo(3) .and. ie(3) <= p%G%hi(3))) then
#endif
                dpdn = -rho%f(ie(1),ie(2),ie(3))*((forcElem%C%A(1) - g(1))*forcElem%n(1)  + &
                                                  (forcElem%C%A(2) - g(2))*forcElem%n(2))*probe_sign
#if DIM==3
                dpdn = dpdn - &
                        rho%f(ie(1),ie(2),ie(3))*((forcElem%C%A(3) - g(3))*forcElem%n(3))*probe_sign
#endif
                ! Build the array of positions and field values in the support domain
                q = 1
                xs = 0.0_dp
                pk = 0.0_dp
                kk = 1
#if DIM==3
                do sk = -1,1
                    kk = ie(3) + sk
#endif
                    do sj = -1,1
                        jj = ie(2) + sj
                        do si = -1,1
                            ii = ie(1) + si
                            xs(1,q)   = p%G%x(ie(1)) + si*delta
                            xs(2,q)   = p%G%y(ie(2)) + sj*delta
                            pk(q)     = p%f(ii,jj,kk)
#if DIM==3
                            xs(3,q)   = p%G%z(ie(3)) + sk*delta
#endif
                            q = q + 1
                        end do
                    end do
#if DIM==3
                end do
#endif

                !call build_support_domain(p, )

                ! Evaluate shape function for the probe
                ds = alpha_s*delta
                phi = get_Phi(Xp1, xs, ds, .false.)

                ! Interpolate pressure and viscous stresses on probe 1
                do q = 1, Ne
                    pp1 = pp1 + phi(1,q)*pk(q)
                end do
#ifdef MPI
            endif ! exit rank
#endif
            ! Find the closest Eulerian x-face to probe 1 and interpolate u
            ie = p%G%closest_grid_node(Xp1, 1)
#ifdef MPI
            ! select the proper rank
            if ( (ie(2) >= p%G%lo(2) .and. ie(2) <= p%G%hi(2)) .and. &
                 (ie(3) >= p%G%lo(3) .and. ie(3) <= p%G%hi(3))) then
#endif
                ! Build the array of positions and field values in the support domain
                q = 1
                xs = 0.0_dp
                uk = 0.0_dp
                kk = 1
#if DIM==3
                do sk = -1,1
                    kk = ie(3) + sk
#endif
                    do sj = -1,1
                        jj = ie(2) + sj
                        do si = -1,1
                            ii = ie(1) + si
                            xs(1,q)   = p%G%x(ie(1)) + 0.5_dp*delta + si*delta
                            xs(2,q)   = p%G%y(ie(2)) + sj*delta
                            uk(q)     = v%x%f(ii,jj,kk)
#if DIM==3
                            xs(3,q)   = p%G%z(ie(3)) + sk*delta
#endif
                            q = q + 1
                        end do
                    end do
#if DIM==3
                end do
#endif
                ! Evaluate shape function for the probe
                ds = alpha_s*delta
                phi = get_Phi(Xp1, xs, ds, .true.)

                ! Interpolate pressure and viscous stresses on probe 1
                do q = 1, Ne
                    Up1(:) = Up1(:) + phi(:,q)*uk(q)
                end do
            endif ! exit rank
            
            ! Find the closest Eulerian y-face to probe 1 and interpolate v
            ie = p%G%closest_grid_node(Xp1, 2)
#ifdef MPI
            ! select the proper rank
            if ( (ie(2) >= p%G%lo(2) .and. ie(2) <= p%G%hi(2)) .and. &
                 (ie(3) >= p%G%lo(3) .and. ie(3) <= p%G%hi(3))) then
#endif
                ! Build the array of positions and field values in the support domain
                q = 1
                xs = 0.0_dp
                uk = 0.0_dp
                kk = 1
#if DIM==3
                do sk = -1,1
                    kk = ie(3) + sk
#endif
                    do sj = -1,1
                        jj = ie(2) + sj
                        do si = -1,1
                            ii = ie(1) + si
                            xs(1,q)   = p%G%x(ie(1)) + si*delta
                            xs(2,q)   = p%G%y(ie(2)) + 0.5_dp*delta + sj*delta
                            uk(q)     = v%y%f(ii,jj,kk)
#if DIM==3
                            xs(3,q)   = p%G%z(ie(3)) + sk*delta
#endif
                            q = q + 1
                        end do
                    end do
#if DIM==3
                end do
#endif
                ! Evaluate shape function for the probe
                ds = alpha_s*delta
                phi = get_Phi(Xp1, xs, ds, .true.)

                ! Interpolate pressure and viscous stresses on probe 1
                do q = 1, Ne
                    Vp1(:) = Vp1(:) + phi(:,q)*uk(q)
                end do
            endif ! exit rank

#if DIM==3
            ! Find the closest Eulerian z-face to probe 1 and interpolate w
            ie = p%G%closest_grid_node(Xp1, 3)
#ifdef MPI
            ! select the proper rank
            if ( (ie(2) >= p%G%lo(2) .and. ie(2) <= p%G%hi(2)) .and. &
                 (ie(3) >= p%G%lo(3) .and. ie(3) <= p%G%hi(3))) then
#endif
                ! Build the array of positions and field values in the support domain
                q = 1
                xs = 0.0_dp
                uk = 0.0_dp
                kk = 1

                do sk = -1,1
                    kk = ie(3) + sk
                    do sj = -1,1
                        jj = ie(2) + sj
                        do si = -1,1
                            ii = ie(1) + si
                            xs(1,q)   = p%G%x(ie(1)) + si*delta
                            xs(2,q)   = p%G%y(ie(2)) + sj*delta
                            uk(q)     = v%z%f(ii,jj,kk)
                            xs(3,q)   = p%G%z(ie(3)) + 0.5_dp*delta + sk*delta
                            q = q + 1
                        end do
                    end do
                end do

                ! Evaluate shape function for the probe
                ds = alpha_s*delta
                phi = get_Phi(Xp1, xs, ds, .true.)

                ! Interpolate pressure and viscous stresses on probe 1
                do q = 1, Ne
                    Wp1(:) = Wp1(:) + phi(:,q)*uk(q)
                end do
#ifdef MPI
            endif ! exit rank
#endif
#endif

            if (nprobes > 1) then
                ! **** Second probe ****************************************************************
                pp2 = 0.0_dp
                Up2 = 0.0_dp
                Vp2 = 0.0_dp
#if DIM==3
                Wp2 = 0.0_dp
#endif

                Xp2(1:Ndim) = forcElem%C%X + forcElem%n*hp2*probe_sign 

                ! Check periodicity and in case translate it
                !call traslate(Xp2, v%G)

                ! Find the closest Eulerian cell center to probe 1 and interpolate pressure
                ie = p%G%closest_grid_node(Xp2, 0)
#ifdef MPI
                    ! select the proper rank
                    if ( (ie(2) >= p%G%lo(2) .and. ie(2) <= p%G%hi(2)) .and. &
                         (ie(3) >= p%G%lo(3) .and. ie(3) <= p%G%hi(3))) then
#endif
                        ! Build the array of positions and field values in the support domain
                        q = 1
                        xs = 0.0_dp
                        pk = 0.0_dp
                        kk = 1
#if DIM==3
                        do sk = -1,1
                            kk = ie(3) + sk
#endif
                            do sj = -1,1
                                jj = ie(2) + sj
                                do si = -1,1
                                    ii = ie(1) + si
                                    xs(1,q)   = p%G%x(ie(1)) + si*delta
                                    xs(2,q)   = p%G%y(ie(2)) + sj*delta
                                    pk(q)     = p%f(ii,jj,kk)
#if DIM==3
                                    xs(3,q)   = p%G%z(ie(3)) + sk*delta
#endif
                                    q = q + 1
                                end do
                            end do
#if DIM==3
                        end do
#endif
                        ! Evaluate shape function for the probe
                        ds = alpha_s*delta
                        phi = get_Phi(Xp2, xs, ds, .false.)
                
                        ! Interpolate pressure and viscous stresses on probe 1
                        do q = 1, Ne
                            pp2 = pp2 + phi(1,q)*pk(q)
                        end do
#ifdef MPI
                    endif ! exit rank
#endif
                    ! Find the closest Eulerian x-face to probe 1 and interpolate u
                    ie = p%G%closest_grid_node(Xp2, 1)
#ifdef MPI
                    ! select the proper rank
                    if ( (ie(2) >= p%G%lo(2) .and. ie(2) <= p%G%hi(2)) .and. &
                         (ie(3) >= p%G%lo(3) .and. ie(3) <= p%G%hi(3))) then
#endif
                        ! Build the array of positions and field values in the support domain
                        q = 1
                        xs = 0.0_dp
                        uk = 0.0_dp
                        kk = 1
#if DIM==3
                        do sk = -1,1
                            kk = ie(3) + sk
#endif
                            do sj = -1,1
                                jj = ie(2) + sj
                                do si = -1,1
                                    ii = ie(1) + si
                                    xs(1,q)   = p%G%x(ie(1)) + 0.5_dp*delta + si*delta
                                    xs(2,q)   = p%G%y(ie(2)) + sj*delta
                                    uk(q)     = v%x%f(ii,jj,kk)
#if DIM==3
                                    xs(3,q)   = p%G%z(ie(3)) + sk*delta
#endif
                                    q = q + 1
                                end do
                            end do
#if DIM==3
                        end do
#endif
                        ! Evaluate shape function for the probe
                        ds = alpha_s*delta
                        phi = get_Phi(Xp2, xs, ds, .true.)
        
                        ! Interpolate pressure and viscous stresses on probe 1
                        do q = 1, Ne
                            Up2(:) = Up2(:) + phi(:,q)*uk(q)
                        end do
#ifdef MPI
                    endif ! exit rank
#endif                     
                    ! Find the closest Eulerian y-face to probe 1 and interpolate v
                    ie = p%G%closest_grid_node(Xp2, 2)
#ifdef MPI
                    ! select the proper rank
                    if ( (ie(2) >= p%G%lo(2) .and. ie(2) <= p%G%hi(2)) .and. &
                         (ie(3) >= p%G%lo(3) .and. ie(3) <= p%G%hi(3))) then
#endif
                        ! Build the array of positions and field values in the support domain
                        q = 1
                        xs = 0.0_dp
                        uk = 0.0_dp
                        kk = 1
#if DIM==3
                        do sk = -1,1
                            kk = ie(3) + sk
#endif
                            do sj = -1,1
                                jj = ie(2) + sj
                                do si = -1,1
                                    ii = ie(1) + si
                                    xs(1,q)   = p%G%x(ie(1)) + si*delta
                                    xs(2,q)   = p%G%y(ie(2)) + 0.5_dp*delta + sj*delta
                                    uk(q)     = v%y%f(ii,jj,kk)
#if DIM==3
                                    xs(3,q)   = p%G%z(ie(3)) + sk*delta
#endif
                                    q = q + 1
                                end do
                            end do
#if DIM==3
                        end do
#endif
                        ! Evaluate shape function for the probe
                        ds = alpha_s*delta
                        phi = get_Phi(Xp2, xs, ds, .true.)
                
                        ! Interpolate pressure and viscous stresses on probe 1
                        do q = 1, Ne
                            Vp2(:) = Vp2(:) + phi(:,q)*uk(q)
                        end do
#ifdef MPI
                    endif ! exit rank
#endif

#if DIM==3
                    ! Find the closest Eulerian z-face to probe 1 and interpolate w
                    ie = p%G%closest_grid_node(Xp2, 3)
#ifdef MPI
                    ! select the proper rank
                    if ( (ie(2) >= p%G%lo(2) .and. ie(2) <= p%G%hi(2)) .and. &
                         (ie(3) >= p%G%lo(3) .and. ie(3) <= p%G%hi(3))) then
#endif
                       ! Build the array of positions and field values in the support domain
                        q = 1
                        xs = 0.0_dp
                        uk = 0.0_dp
                        kk = 1
                
                        do sk = -1,1
                            kk = ie(3) + sk
                            do sj = -1,1
                                jj = ie(2) + sj
                                do si = -1,1
                                    ii = ie(1) + si
                                    xs(1,q)   = p%G%x(ie(1)) + si*delta
                                    xs(2,q)   = p%G%y(ie(2)) + sj*delta
                                    uk(q)     = v%z%f(ii,jj,kk)
                                    xs(3,q)   = p%G%z(ie(3)) + sk*delta
                                    q = q + 1
                                end do
                            end do
                        end do

                        ! Evaluate shape function for the probe
                        ds = alpha_s*delta
                        phi = get_Phi(Xp2, xs, ds, .true.)
        
                        ! Interpolate pressure and viscous stresses on probe 1
                        do q = 1, Ne
                            Wp2(:) = Wp2(:) + phi(:,q)*uk(q)
                        end do
#ifdef MPI
                    endif ! exit rank
#endif
#endif
                
                ! **** Surface values***************************************************************
                ! Extrapolate shear rate from probes on surface
                dUdx = (Up1(2) - a*Up2(2))/(1.0_dp - a)
                dUdy = (Up1(3) - a*Up2(3))/(1.0_dp - a)
                dVdx = (Vp1(2) - a*Vp2(2))/(1.0_dp - a)
                dVdy = (Vp1(3) - a*Vp2(3))/(1.0_dp - a)
#if DIM==3
                dUdz = (Up1(4) - a*Up2(4))/(1.0_dp - a)
                dVdz = (Vp1(4) - a*Vp2(4))/(1.0_dp - a)
                dWdx = (Wp1(2) - a*Wp2(2))/(1.0_dp - a)
                dWdy = (Wp1(3) - a*Wp2(3))/(1.0_dp - a)
                dWdz = (Wp1(4) - a*Wp2(4))/(1.0_dp - a)
#endif
                ! Extrapolate pressure
                pl(l) = pp1/hp1**2 - pp2/hp2**2 - dpdn*(1.0_dp/hp1 - 1.0_dp/hp2)
                pl(l) = pl(l)/(1.0_dp/hp1**2 - 1.0_dp/hp2**2)
            else
                dUdx = Up1(2)
                dUdy = Up1(3)
                dVdx = Vp1(2)
                dVdy = Vp1(3)
#if DIM==3
                dUdz = Up1(4)
                dVdz = Vp1(4)
                dWdx = Wp1(2)
                dWdy = Wp1(3)
                dWdz = Wp1(4)
#endif
                pl(l) = pp1 - hp1*dpdn
            endif

            ! Compute shear stress on surface
            ie = v%G%closest_grid_node(forcElem%C%X, 0)
            if ((ie(2) >= v%G%lo(2) .and. ie(2) <= v%G%hi(2)) .and. &
                (ie(3) >= v%G%lo(3) .and. ie(3) <= v%G%hi(3))) then

                taul(1,l) = mu%f(ie(1), ie(2), ie(3))*(dUdx + dUdx)
                taul(2,l) = mu%f(ie(1), ie(2), ie(3))*(dUdy + dVdx)
                taul(3,l) = mu%f(ie(1), ie(2), ie(3))*(dVdy + dVdy)
#if DIM==3
                taul(4,l) = mu%f(ie(1), ie(2), ie(3))*(dUdz + dWdx)
                taul(5,l) = mu%f(ie(1), ie(2), ie(3))*(dVdz + dWdy)
                taul(6,l) = mu%f(ie(1), ie(2), ie(3))*(dWdz + dWdz)
#endif
            endif
        end do

        ! MPI communication for the integration
        call mpi_allreduce(mpi_in_place, taul(1,:), Nfe, mpi_real8, mpi_sum, mpi_comm_world, ierror)
        call mpi_allreduce(mpi_in_place, taul(2,:), Nfe, mpi_real8, mpi_sum, mpi_comm_world, ierror)
        call mpi_allreduce(mpi_in_place, taul(3,:), Nfe, mpi_real8, mpi_sum, mpi_comm_world, ierror)
#if DIM==3
        ! MPI communication for the integration
        call mpi_allreduce(mpi_in_place, taul(4,:), Nfe, mpi_real8, mpi_sum, mpi_comm_world, ierror)
        call mpi_allreduce(mpi_in_place, taul(5,:), Nfe, mpi_real8, mpi_sum, mpi_comm_world, ierror)
        call mpi_allreduce(mpi_in_place, taul(6,:), Nfe, mpi_real8, mpi_sum, mpi_comm_world, ierror)
#endif
    end subroutine
    !===============================================================================================

    !===============================================================================================
    subroutine build_support_domain(s, ie, loc, xk, sk)

        use global_mod, only : Ndim, stagger
        use scalar_mod, only : scalar
        use mls_mod   , only : Ne

        ! In/Out variables
        type(scalar), intent(in   ) :: s           ! input scalar field
        integer     , intent(in   ) :: ie(3)       ! logical coordinates of closest point
        integer     , intent(in   ) :: loc         ! location of the scalar in the grid
        real(dp)    , intent(inout) :: xk(Ndim,Ne) ! suport domain nodes
        real(dp)    , intent(inout) :: sk(Ne)      ! support domain values

        ! Local variables
        integer :: i, j, k, q

        q = 1
        k = 1
#if DIM==3
        do k = ie(3)-1,ie(3)+1
#endif
            do j = ie(2)-1,ie(2)+1
                do i = ie(1)-1,ie(1)+1
                    xk(1,q) = s%G%x(i) + stagger(1,loc)*s%G%delta
                    xk(2,q) = s%G%y(j) + stagger(2,loc)*s%G%delta
#if DIM==3
                    xk(3,q) = s%G%z(k) + stagger(3,loc)*s%G%delta
#endif
                    sk(  q) = s%f(i,j,k)
                    q       = q + 1
                end do
            end do
#if DIM==3
        end do
#endif

    end subroutine build_support_domain
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
