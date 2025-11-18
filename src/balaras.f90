module balaras_mod

    ! This module contains procedures for stress model of Wang et al JCP 2019.

    use precision_mod       , only : dp
    use scalar_mod          , only : scalar
    use vector_mod          , only : vector
    use lagrangian_solid_mod, only : lagrangian_solid, forcing_element

    implicit none

    private
    public :: compute_hydrodynamic_loadsB

contains

    !===============================================================================================
    subroutine compute_hydrodynamic_loadsB(obj, v, p, mu, rho, g)

        ! Compute hydrodynamic loads on solid object obj given by the velocity 
        ! field v, pressure field p, of the flow with viscosity mu, density rho 
        ! and body force g.

        use global_mod   , only : Ndim, tdof, rdof, dofs
        use euclidean_mod, only : crossProduct

        ! In/Out variables
        real(dp)               , intent(in   )         :: g(:)
        type(scalar)           , intent(in   )         :: p, mu, rho
        type(vector)           , intent(in   )         :: v
        class(lagrangian_solid), intent(inout), target :: obj

        ! Local variables
        integer                               :: l, Nfe
        real(dp)                              :: Xcm(Ndim), r(tdof), Fh_x_r(3)
        type(forcing_element)   , pointer     :: forcElem
        real(dp), dimension(:)  , allocatable :: pl_np, pl_nm
        real(dp), dimension(:,:), allocatable :: fvl_np, fvl_nm

        ! Number of forcing Elements
        Nfe = size(obj%forcing_elements)

        ! Allocate local stress, the size is the number of forcing elements
        allocate(pl_np(Nfe))
        allocate(fvl_np(3,Nfe))

        ! Evaluate stresses on the surface using the probe on the positive direction of the norm
        call compute_stresses(obj, v, p, rho, mu, g, pl_np, fvl_np, 1)

        ! If the structure is open it is necessary to compute forces also on the negative side of
        ! the normal
        if (obj%is_open) then
            ! Allocate memory
            allocate(pl_nm(Nfe))
            allocate(fvl_nm(3, Nfe))

            call compute_stresses(obj, v, p, rho, mu, g, pl_nm, fvl_nm, -1)

            ! Add tau and p on both sides
            pl_np = pl_np - pl_nm
            fvl_np = fvl_np - fvl_nm
        end if

        ! store pressure for output purpose
        obj%p_s = pl_np

        ! Get local forces from stresses
        Xcm = obj%center_of_mass%X(1:tdof)
        do l = 1,Nfe
            obj%tau_s(1:3,l) = fvl_np(:,l)      ! save viscous forces in tau_s
            forcElem => obj%forcing_elements(l)

            ! Shear stress forces
            forcElem%C%Fv(1) = fvl_np(1,l)*forcElem%A
            forcElem%C%Fv(2) = fvl_np(2,l)*forcElem%A

            ! Pressure force
            forcElem%C%Fp(1) = -pl_np(l)*forcElem%n(1)*forcElem%A
            forcElem%C%Fp(2) = -pl_np(l)*forcElem%n(2)*forcElem%A
#if DIM==3
            forcElem%C%Fv(3) = fvl_np(3,l)*forcElem%A
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
        deallocate(fvl_np, pl_np)
        if (obj%is_open) deallocate(fvl_nm, pl_nm)

    end subroutine compute_hydrodynamic_loadsB
    !===============================================================================================

    !===============================================================================================
    subroutine compute_stresses(obj, v, p, rho, mu, g, pl, fvl, probe_sign)

        ! Objective: evaluate stresses on solid surface.

        use mpi
        use global_mod   , only : Ndim, ierror, myrank
        use euclidean_mod, only : dotProduct, distance, ZEROV
        use mls_mod

        ! In/Out variables
        class(lagrangian_solid), intent(in   ), target :: obj        
        real(dp)               , intent(in   )         :: g(:)
        type(scalar)           , intent(in   )         :: p, rho, mu
        type(vector)           , intent(in   )         :: v
        real(dp)               , intent(inout)         :: fvl(:,:)
        real(dp)               , intent(inout)         :: pl(:)
        integer                , intent(in   )         :: probe_sign

        ! Local variables
        integer                        :: Nfe, Nsp, l, e, ie(3), i, j, k, q
        real(dp)                       :: delta, h, dpdeta_e, dpdeta_m, pe, ue(3), t(3), n(3)
        real(dp)                       :: duxsideta_e, d2uxsideta2_e, d3uxsideta3_e, visc
        real(dp)                       :: xs(Ndim,Ne), sk(Ne), ds(Ndim,Ne), phi(m,Ne)
        real(dp)                       :: Xe(3), Fp(Ndim+1)
        real(dp)                       :: Xe_s(3,-1:1,-1:1,-1:1), Pe_s(-1:1,-1:1,-1:1)
        real(dp)                       :: Ve_s(3,-1:1,-1:1,-1:1)
        real(dp)                       :: d1P(3), d2P(3,3), d1V(3,3), d2V(3,3,3)

        type(forcing_element), pointer :: forcElem

        ! probe distance
        delta = v%G%delta
        h = delta

        ! Number of forcing elements
        Nfe = size(obj%forcing_elements)

        ! Points array size
        Nsp = merge(27,9,Ndim==3)

        ! local viscosity (constant for now)
        visc = mu%f(p%G%lo(1), p%G%lo(2), p%G%lo(3))

        ! Set stresses to zero
        pl = 0.0_dp
        fvl = 0.0_dp

        ! Cycle over all the forcing elements of the solid
        do l = 1,Nfe

            ! Select the local forcing elemeng
            forcElem => obj%forcing_elements(l)
            n = 0.0_dp
            n(1:Ndim) = forcElem%n(1:Ndim)

            ! Generate the probes in the given normal direction
            Xe(1:Ndim) = forcElem%C%X + forcElem%n*h*probe_sign 

            ! Normal pressure gradient on the lagrangian marker, eq 19 of Wang et al JCP 2019
#ifdef MPI
            ! select the proper rank with probe inside
            if ( (ie(2) >= p%G%lo(2) .and. ie(2) <= p%G%hi(2)) .and. &
                 (ie(3) >= p%G%lo(3) .and. ie(3) <= p%G%hi(3))) then
#endif
                dpdeta_m = -rho%f(ie(1),ie(2),ie(3))*((forcElem%C%A(1) - g(1))*n(1)  + &
                                                      (forcElem%C%A(2) - g(2))*n(2))*probe_sign
#if DIM==3
                dpdeta_m = dpdeta_m - rho%f(ie(1),ie(2),ie(3))*((forcElem%C%A(3) - g(3))*n(3))*probe_sign
#endif
            else
                dpdeta_m = 0.0_dp
            endif

            ! build the points array for the probe
            k = 0
#if DIM==3
            do k = -1,1
#endif
                do j = -1,1
                    do i = -1,1
                        Xe_s(:,i,j,k) = Xe + [i*delta, j*delta, k*delta]
                    end do
                end do
#if DIM==3
            end do
#endif

            ! Now for each point inside the points array we need to evaluate pressure 
            ! and velocity
            Pe_s = 0.0_dp
            Ve_s = 0.0_dp
            k = 0
#if DIM==3
            do k = -1,1
#endif
                do j = -1,1
                    do i = -1,1

                        ! **** Pressure ****
                        ie = p%G%closest_grid_node(Xe_s(:,i,j,k), 0)
                        ! select the proper rank
                        if ( (ie(2) >= p%G%lo(2) .and. ie(2) <= p%G%hi(2)) .and. &
                             (ie(3) >= p%G%lo(3) .and. ie(3) <= p%G%hi(3))) then

                            call build_support_domain(p, ie, 0, xs, sk)                         

                            ! Evaluate shape function for the probe
                            ds = alpha_s*delta
                            phi = get_Phi(Xe_s(:,i,j,k), xs, ds, .false.)

                            ! Interpolate pressure
                            do q = 1, Ne
                                Pe_s(i,j,k) = Pe_s(i,j,k) + phi(1,q)*sk(q)
                            end do
                        endif

                        ! **** Vx ****
                        ie = p%G%closest_grid_node(Xe_s(:,i,j,k), 1)
                        ! select the proper rank
                        if ( (ie(2) >= p%G%lo(2) .and. ie(2) <= p%G%hi(2)) .and. &
                             (ie(3) >= p%G%lo(3) .and. ie(3) <= p%G%hi(3))) then

                            call build_support_domain(v%x, ie, 1, xs, sk)                         

                            ! Evaluate shape function for the probe
                            ds = alpha_s*delta
                            phi = get_Phi(Xe_s(:,i,j,k), xs, ds, .false.)

                            ! Interpolate velocity
                            do q = 1, Ne
                                Ve_s(1,i,j,k) = Ve_s(1,i,j,k) + phi(1,q)*sk(q)
                            end do
                        endif

                        ! **** Vy ****
                        ie = p%G%closest_grid_node(Xe_s(:,i,j,k), 2)
                        ! select the proper rank
                        if ( (ie(2) >= p%G%lo(2) .and. ie(2) <= p%G%hi(2)) .and. &
                             (ie(3) >= p%G%lo(3) .and. ie(3) <= p%G%hi(3))) then

                            call build_support_domain(v%y, ie, 2, xs, sk)

                            ! Evaluate shape function for the probe
                            ds = alpha_s*delta
                            phi = get_Phi(Xe_s(:,i,j,k), xs, ds, .false.)

                            ! Interpolate velocity
                            do q = 1, Ne
                                Ve_s(2,i,j,k) = Ve_s(2,i,j,k) + phi(1,q)*sk(q)
                            end do
                        endif
#if DIM==3

                        ! **** Vz ****
                        ie = p%G%closest_grid_node(Xe_s(:,i,j,k), 3)
                        ! select the proper rank
                        if ( (ie(2) >= p%G%lo(2) .and. ie(2) <= p%G%hi(2)) .and. &
                             (ie(3) >= p%G%lo(3) .and. ie(3) <= p%G%hi(3))) then

                            call build_support_domain(v%z, ie, 3, xs, sk)

                            ! Evaluate shape function for the probe
                            ds = alpha_s*delta
                            phi = get_Phi(Xe_s(:,i,j,k), xs, ds, .false.)

                            ! Interpolate velocity
                            do q = 1, Ne
                                Ve_s(3,i,j,k) = Ve_s(3,i,j,k) + phi(1,q)*sk(q)
                            end do
                        endif
#endif

                    end do
                end do
#if DIM==3
            end do
#endif

            call mpi_allreduce(mpi_in_place, dpdeta_m, 1, mpi_real8, mpi_sum, mpi_comm_world, ierror)
            call mpi_allreduce(mpi_in_place,  Pe_s, Ne, mpi_real8, mpi_sum, mpi_comm_world, ierror)
            call mpi_allreduce(mpi_in_place, Ve_s(1,:,:,:), Ne, mpi_real8, mpi_sum, mpi_comm_world, ierror)
            call mpi_allreduce(mpi_in_place, Ve_s(2,:,:,:), Ne, mpi_real8, mpi_sum, mpi_comm_world, ierror)
#if DIM==3
            call mpi_allreduce(mpi_in_place, Ve_s(3,:,:,:), Ne, mpi_real8, mpi_sum, mpi_comm_world, ierror)
#endif

            Pe = Pe_s(0,0,0)
            ue = Ve_s(:,0,0,0)

            ! Pressure derivatives in e
            d1P(1)   = 0.5_dp*(Pe_s( 1, 0, 0) - Pe_s(-1, 0, 0))/delta ! dpdx
            d1P(2)   = 0.5_dp*(Pe_s( 0, 1, 0) - Pe_s( 0,-1, 0))/delta ! dpdy
            d2P(1,1) = (Pe_s(1,0,0) + Pe_s(-1,0,0) - 2.0_dp*Pe_s(0,0,0))/delta**2 ! d2pdx2
            d2P(2,1) = 0.25_dp*(Pe_s(1,1,0) + Pe_s(-1,1,0) - Pe_s(1,-1,0) - Pe_s(-1,1,0))/delta**2 ! d2pdxdy
            d2P(1,2) = d2P(2,1)
            d2P(2,2) = (Pe_s(0,1,0) + Pe_s(0,-1,0) - 2.0_dp*Pe_s(0,0,0))/delta**2 !d2pdy2
#if DIM==3
            d1P(3)   = 0.5_dp*(Pe_s(0,0,1) - Pe_s(0,0,-1)) ! dpdz
            d2P(1,3) = 0.25_dp*(Pe_s(1,0,1) + Pe_s(-1,0,-1) - Pe_s(1,0,-1) - Pe_s(-1,0,1))/delta**2 ! d2pdxdz
            d2P(2,3) = 0.25_dp*(Pe_s(0,1,1) + Pe_s(0,-1,-1) - Pe_s(0,-1,1) - Pe_s(0,1,-1))/delta**2 ! d2pdydz
            d2P(3,1) = d2P(1,3)
            d2P(3,2) = d2P(2,3)
            d2P(3,3) = (Pe_s(0,0,1) + Pe_s(0,0,-1) - 2.0_dp*Pe_s(0,0,0))/delta**2 !d2pdz2
#endif

            ! Velocity derivatives in e
            do i = 1,Ndim
                d1V(i,1)   = 0.5_dp*(Ve_s(i,1,0,0) - Ve_s(i,-1,0,0))/delta ! duidx
                d1V(i,2)   = 0.5_dp*(Ve_s(i,0,1,0) - Ve_s(i,0,-1,0))/delta ! duidy
                d2V(i,1,1) = (Ve_s(i,1,0,0) + Ve_s(i,-1,0,0) - 2.0_dp*Ve_s(i,0,0,0))/delta**2 ! d2uidx2
                d2V(i,2,2) = (Ve_s(i,0,1,0) + Ve_s(i,0,-1,0) - 2.0_dp*Ve_s(i,0,0,0))/delta**2 ! d2uidy2
                d2V(i,1,2) = 0.25_dp*(Ve_s(i,1,1,0) + Ve_s(i,-1,-1,0) - Ve_s(i,-1,1,0) - Ve_s(i,1,-1,0))/delta**2 ! d2uidxdy
                d2V(i,2,1) = d2V(i,1,2)
#if DIM==3
                d1V(i,3)   = 0.5_dp*(Ve_s(i,0,0,1) - Ve_s(i,0,0,-1))/delta ! duidz
                d2V(i,3,3) = (Ve_s(i,0,0,1) + Ve_s(i,0,0,-1) - 2.0_dp*Ve_s(i,0,0,0))/delta**2 ! d2uidz2
                d2V(i,1,3) = 0.25_dp*(Ve_s(i,1,0,1) + Ve_s(i,-1,0,-1) - Ve_s(i,-1,0,1) - Ve_s(i,1,0,-1))/delta**2 ! d2uidxdz
                d2V(i,2,3) = 0.25_dp*(Ve_s(i,0,1,1) + Ve_s(i,0,-1,-1) - Ve_s(i,0,-1,1) - Ve_s(i,0,1,-1))/delta**2 ! d2uidydz
                d2V(i,3,1) = d2V(i,1,3)
                d2V(i,3,2) = d2V(i,2,3)
#endif
            end do

            ! Tangent vector
#if DIM==3
            t = ue - forcElem%C%V - dotProduct(ue - forcElem%C%V, forcElem%n)*n
            !t = ue - dotProduct(ue, forcElem%n)*n
            t = t / (distance(t, ZEROV) + 1.0e-12_dp)
#else
            t(1) = -n(2)
            t(2) =  n(1)
            t(3) = 0.0_dp
#endif

            ! Approximated derivatives:
            dpdeta_e      = 0.0_dp
            duxsideta_e   = 0.0_dp
            d2uxsideta2_e = 0.0_dp
            d3uxsideta3_e = 0.0_dp
            do i = 1,Ndim
                ! eq 27a
                dpdeta_e = dpdeta_e + d1P(i)*n(i)

                do j = 1,Ndim
                    ! eq 27b
                    duxsideta_e = duxsideta_e + d1V(i,j)*t(i)*n(j)
            
                    ! eq 27 d
                    d3uxsideta3_e = d3uxsideta3_e + d2P(i,j)*n(i)*n(j)/visc

                    do k = 1,Ndim
                        ! eq 27c
                        d2uxsideta2_e = d2uxsideta2_e + d2V(i,j,k)*t(i)*n(j)*n(k)
                    end do
                end do
            end do

            ! Surface normal stress, equation 20 of Wang et al JCP 2019
            pl(l) = pe - 0.5_dp*(dpdeta_m + dpdeta_e)*h

            ! Surface viscous stress, equation 26 of Wang et al JCP 2019 (with the proper sign 
            ! for the last term)
            fvl(:,l) = visc*(duxsideta_e - d2uxsideta2_e*h + 0.5_dp*d3uxsideta3_e*h**2)*t

        end do


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

end module
