module mls_mod

    ! This moudule contains all the procedures for the implementation of the
    ! Moving Least Square method.
    
    use precision_mod, only : dp
    use global_mod   , only : Ndim

    implicit none

    ! Number of points in the support domain
    integer, parameter :: Ne = merge(9, 27, Ndim == 2)

    ! Dimension of the basis vector for the interpolation
    integer, parameter :: m = merge(3, 4, Ndim == 2)
    
    ! MLS Shape function
    real(dp) :: phi(m,Ne)

    ! Procedure for weight computation. Available functions are:
    ! - Weight_W1: cubic spline (default)
    ! - Weight_W2: quartic spline
    ! - Weight_W3: exponential
    procedure(weight_function), pointer :: compute_weight => Weight_W1

    interface
        subroutine weight_function(dif, ds, w)
            use precision_mod, only : dp
            use global_mod   , only : Ndim
            import Ne, m
            real(dp), intent(in ) :: dif(Ndim, Ne), ds(Ndim, Ne)
            real(dp), intent(out) :: w(Ne,m)
        end subroutine weight_function
    end interface

contains

    !===============================================================================================
    function interpolate(f, xl, ie, ind) result(fl)

        ! This function interpolate the value of the scalar field f defined on 
        ! the eulerian grid in the point xl on the lagrangian grid.

        use mpi
        use scalar_mod
    
        ! In/Out variables
        integer     , intent(in) :: ie(3), ind
        real(dp)    , intent(in) :: xl(Ndim)
        type(scalar), intent(in) :: f
        
        ! Local variables
        integer  :: si, sj, q, ii, jj, kk, ierror
        real(dp) :: xe, ye, xs(Ndim,Ne), fl(m), ds(Ndim,Ne), fk(Ne), stagger(3)
#if DIM==3
        integer  :: sk
        real(dp) :: ze
#endif
        select case (ind)
        case(0)
             stagger = 0.0_dp
        case(1)
            stagger = [f%G%delta*0.5_dp, 0.0_dp, 0.0_dp]
        case(2)
            stagger = [0.0_dp, f%G%delta*0.5_dp, 0.0_dp]
        case(3)
            stagger = [0.0_dp, 0.0_dp, f%G%delta*0.5_dp]
        end select

        ! Select the closest Eulerian point
        xe = f%G%x(ie(1)) + stagger(1)
        ye = f%G%y(ie(2)) + stagger(2)
#if DIM==3
        ze = f%G%z(ie(3)) + stagger(3)
#endif
        ! Select the rank containing the Eulerian point
        rank_with_point: if ( (ie(2) >= f%G%lo(2) .and. ie(2) <= f%G%hi(2)) .and. &
                              (ie(3) >= f%G%lo(3) .and. ie(3) <= f%G%hi(3))) then 
            
            ! Build the array of positions and f values in the support domain
            q = 1
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
                        xs(1,q) = xe + si*f%G%delta
                        xs(2,q) = ye + sj*f%G%delta
                        fk(q) = f%f(ii,jj,kk) 
                        ds(1,q) = merge(f%G%delta, 2.0_dp*f%G%delta, si == 0)
                        ds(2,q) = merge(f%G%delta, 2.0_dp*f%G%delta, sj == 0)
#if DIM==3
                        xs(3,q) = ze + sk*f%G%delta
                        ds(3,q) = merge(f%G%delta, 2.0_dp*f%G%delta, sk == 0)
#endif
                        q = q + 1
                    end do
                end do
#if DIM==3
            end do
#endif      

            ! Size of the support domain
            ds = 1.60_dp*f%G%delta
            
            ! Compute shape function phi
            call compute_phi(xl, xs, ds)

            ! Compute interpolated value in xl
            fl = 0.0_dp
            do q = 1, Ne
                fl(1) = fl(1) + phi(1,q)*fk(q)
                fl(2) = fl(2) + phi(2,q)*fk(q)
                fl(3) = fl(3) + phi(3,q)*fk(q)
#if DIM==3
                fl(4) = fl(4) + phi(4,q)*fk(q)
#endif
            end do

        else

            ! Set to zero on other ranks
            fl = 0.0_dp

        end if rank_with_point

        ! Comunicate the interpolated value
        call mpi_allreduce(mpi_in_place, fl, m, mpi_real8, mpi_sum, &
                            mpi_comm_world, ierror)
        
    end function interpolate
    !===============================================================================================

    !===============================================================================================
    subroutine compute_phi(gpos, x, ds)

        ! Compute Moving-Least-Square shape function and its derivatives. 
        ! Routine adapted from: "An introduction to meshfree methods and their 
        ! programming" by Liu and Gu, chapter 3 appendix.    
        
        ! In/Out variables
        real(dp), intent(in) :: gpos(Ndim), x(Ndim,Ne), ds(Ndim,Ne)

        ! Local variables
        integer  :: i, j, k
        real(dp) :: gp(m,m), A(m,m,m), B(m,Ne,m), c(m), aa(m,m), gam(m,m), ep

        ! First compute the basis vector in the interpolation point gpos
        call Compute_Basis(gpos, gp)

        ! Then compute the matrices A and B
        call Compute_AB(gpos, x, ds, A, B)
        ep = 1.0e-20_dp

        !**** Compute gamma ****************************************************
        gam = 0.0_dp

        c = gp(1,:)
        do i = 1,m
            do j = 1,m
                aa(i,j) = A(i,j,1)
            end do
        end do
        
        c = solve_wbs(ge_wpp(aa,c))
        gam(:,1) = c
        !***********************************************************************

        !**** Compute dgamdx ***************************************************
        do i = 1,m
            c(i) = 0.0_dp
            do j = 1,m
                c(i) = c(i) + A(i,j,2)*gam(j,1)
            end do
        end do
        do k = 1,m
            c(k) = gp(2,k) - c(k)
        end do
        do i = 1,m
            do j = 1,m
                aa(i,j) = A(i,j,1)
            end do
        end do

        c = solve_wbs(ge_wpp(aa,c))
        gam(:,2) = c
        !***********************************************************************

        !**** Compute dgamdy ***************************************************
        do i = 1,m
            c(i) = 0.0_dp
            do j = 1,m
                c(i) = c(i) + A(i,j,3)*gam(j,1)
            end do
        end do
        do k = 1,m
            c(k) = gp(3,k) - c(k)
        end do
        do i = 1,m
            do j = 1,m
                aa(i,j) = A(i,j,1)
            end do
        end do
        c = solve_wbs(ge_wpp(aa,c))
        gam(:,3) = c
        !***********************************************************************

#if DIM==3
        !**** Compute dgammadz *************************************************
        do i = 1,m
            c(i) = 0.0_dp
            do j = 1,m
                c(i) = c(i) + A(i,j,4)*gam(j,1)
            end do
        end do
        do k = 1,m
            c(k) = gp(4,k) - c(k)
        end do
        do i = 1,m
            do j = 1,m
                aa(i,j) = A(i,j,1)
            end do
        end do

        c = solve_wbs(ge_wpp(aa,c))
        gam(:,4) = c
        !***********************************************************************
#endif

        !**** Compute Phi and its derivatives **********************************
        do i = 1,Ne
            do j = 1,m
                phi(j,i) = 0.0_dp
            enddo
            do j = 1,m
                phi(1,i) = phi(1,i) + gam(j,1)*b(j,i,1)
                phi(2,i) = phi(2,i) + gam(j,2)*b(j,i,1) + gam(j,1)*b(j,i,2)
                phi(3,i) = phi(3,i) + gam(j,3)*b(j,i,1) + gam(j,1)*b(j,i,3)
#if DIM==3
                phi(4,i) = phi(4,i) + gam(j,4)*b(j,i,1) + gam(j,1)*b(j,i,4)
#endif   
            end do
        end do
        
    end subroutine compute_phi
    !===============================================================================================

    !===============================================================================================
    subroutine Compute_Basis(gpos, gp)

        ! Evaluate the basis vector gp in the interpolation point gpos
        ! The base vector is gp = [1.0, x, y] (in 2D) and  
        ! gp = [1.0, x, y, z] (in 3D).
        ! Routine adapted from: "An introduction to meshfree methods and their 
        ! programming" by Liu and Gu, chapter 3 appendix.

        ! In/Out variables
        real(dp), intent(in   ) :: gpos(Ndim)
        real(dp), intent(inout) :: gp(m,m)

        gp = 0.0_dp
        gp(1,1) = 1.0_dp
        gp(1,2) = gpos(1)
        gp(1,3) = gpos(2)
        gp(2,2) = 1.0_dp
        gp(3,3) = 1.0_dp
#if DIM==3
        gp(1,4) = gpos(3)
        gp(4,4) = 1.0_dp
#endif
        
    end subroutine Compute_Basis
    !===============================================================================================

    !===============================================================================================
    subroutine Compute_AB(gpos, x, ds, A, B)

        ! Compute A matrix and B matrix and their derivatives. 
        ! Routine adapted from: "An introduction to meshfree methods and their 
        ! programming" by Liu and Gu, chapter 3 appendix.

        ! In/Out varialbes
        real(dp), intent(in ) :: gpos(Ndim), x(Ndim,Ne), ds(Ndim,Ne)
        real(dp), intent(out) :: A(m,m,m), B(m,Ne,m)

        ! Local variables
        integer :: i, j, k, ik, jk, kk
        real(dp) :: W(Ne,m), dif(Ndim,Ne), p(m,Ne), pp(m,m)

        do i = 1,Ne
            p(1,i) = 1.0_dp
            p(2,i) = x(1,i)
            p(3,i) = x(2,i)
            dif(1,i) = gpos(1) - x(1,i)
            dif(2,i) = gpos(2) - x(2,i)
#if DIM==3
            p(4,i) = x(3,i)
            dif(3,i) = gpos(3) - x(3,i)
#endif
        enddo

        ! Compute the weight function the support domain
        call compute_weight(dif, ds, W)

        ! ************* Compute B and its derivatives
        do i = 1,m
            do j = 1,Ne
                do k = 1,m
                    B(i,j,k) = p(i,j)*W(j,k)
                end do
            end do
        end do

        ! ************* Compute A and its derivatives
        A = 0.0_dp
        do i = 1,Ne
            do ik = 1,m
                do jk = 1,m
                    pp(ik,jk) = p(ik,i)*p(jk,i)
                end do
            end do
            do ik = 1,m
                do jk = 1,m
                    do kk = 1,m
                        A(ik,jk,kk) = A(ik,jk,kk) + W(i,kk)*pp(ik,jk)
                    end do
                end do
            end do
        end do

    end subroutine Compute_AB
    !===============================================================================================

    !===============================================================================================
    subroutine Weight_W1(dif, ds, w)

        ! Cubic spline weight function. 
        ! Routine adapted from "An introduction to meshfree methods and their 
        ! programming" by Liu and Gu, chapter 3 appendix.

        ! In/Out variables
        real(dp), intent(in ) :: dif(Ndim, Ne), ds(Ndim, Ne)
        real(dp), intent(out) :: w(Ne,m)

        ! Local variables
        integer  :: i
        real(dp) :: ep, difx, dify, drdx, drdy, rx, ry, wx, wy, dwxdx, dwydy
#if DIM==3
        real(dp) :: difz, drdz, rz, wz, dwzdz
#endif
        
        ep = 1.0e-20_dp

        do i = 1,Ne
            difx = dif(1,i)
            dify = dif(2,i)
            if (dabs(difx) .le. ep) then
                drdx = 0.0_dp
            else
                drdx = (difx/dabs(difx))/ds(1,i)
            end if
            if (dabs(dify).le.ep) then
                drdy = 0.0_dp
            else
                drdy = (dify/dabs(dify))/ds(2,i)
            end if
            
            rx = dabs(dif(1,i))/ds(1,i)
            ry = dabs(dif(2,i))/ds(2,i)

            if (rx > 0.5_dp) then
                wx = (4.0_dp/3.0_dp) - 4.0_dp*rx + 4.0_dp*rx*rx - &
                        (4.0_dp/3.0_dp)*rx*rx*rx
                dwxdx=(-4.0_dp + 8.0_dp*rx - 4.0_dp*rx*rx)*drdx
            else if (rx <= 0.5_dp) then
                wx = (2.0_dp/3.0_dp) - 4.0_dp*rx*rx + 4.0_dp*rx*rx*rx
                dwxdx = (-8.0_dp*rx + 12.0_dp*rx*rx)*drdx
            end if

            if (ry > 0.5_dp) then
                wy = (4.0_dp/3.0_dp) - 4.0_dp*ry + 4.0_dp*ry*ry - &
                        (4.0_dp/3.0_dp)*ry*ry*ry
                dwydy = (-4.0_dp + 8.0_dp*ry - 4.0_dp*ry*ry)*drdy
            else if (ry <= 0.5_dp) then
                wy = (2.0_dp/3.0_dp) - 4.0_dp*ry*ry + 4.0_dp*ry*ry*ry
                dwydy = (-8.0_dp*ry + 12.0_dp*ry*ry)*drdy
            end if

            w(i,1) = wx*wy
            w(i,2) = wy*dwxdx
            w(i,3) = wx*dwydy
#if DIM==3
            difz = dif(3,i)
            if (dabs(difz) .le. ep) then
                drdz = 0.0_dp
            else
                drdz = (difz/dabs(difz))/ds(3,i)
            end if
            rz = dabs(dif(3,i))/ds(3,i)
            if (rz > 0.5_dp) then
                wz = (4.0_dp/3.0_dp) - 4.0_dp*rz + 4.0_dp*rz*rz - &
                        (4.0_dp/3.0_dp)*rz*rz*rz
                dwzdz = (-4.0_dp + 8.0_dp*rz - 4.0_dp*rz*rz)*drdz
            else if (rz .le. 0.5_dp) then
                wz = (2.0_dp/3.0_dp) - 4.0_dp*rz*rz + 4.0_dp*rz*rz*rz
                dwzdz = (-8.0_dp*rz + 12.0_dp*rz*rz)*drdz
            end if
            w(i,1) = wx*wy*wz
            w(i,4) = wz*dwzdz
#endif
        end do
        
    end subroutine Weight_W1
    !===============================================================================================

    !===============================================================================================
    subroutine Weight_W2(dif, ds, w)

        ! Quartic spline weight function. 
        ! Routine adapted from "An introduction to meshfree methods and their 
        ! programming" by Liu and Gu, chapter 3.

        ! In/Out variables
        real(dp), intent(in ) :: dif(Ndim, Ne), ds(Ndim, Ne)
        real(dp), intent(out) :: w(Ne,m)

        ! Local variables
        integer  :: i
        real(dp) :: ep, difx, dify, drdx, drdy, rx, ry, wx, wy, dwxdx, dwydy
#if DIM==3
        real(dp) :: difz, drdz, rz, wz, dwzdz
#endif
        
        ep = 1.0e-20_dp

        do i = 1,Ne
            difx = dif(1,i)
            dify = dif(2,i)
            if (dabs(difx) .le. ep) then
                drdx = 0.0_dp
            else
                drdx = (difx/dabs(difx))/ds(1,i)
            end if
            if (dabs(dify).le.ep) then
                drdy = 0.0_dp
            else
                drdy = (dify/dabs(dify))/ds(2,i)
            end if
            
            rx = dabs(dif(1,i))/ds(1,i)
            ry = dabs(dif(2,i))/ds(2,i)

            if (rx <= 1.0_dp) then
                wx = 1.0_dp - 6.0_dp*rx**2 + 8.0_dp*rx**3 - 3.0_dp*rx**4
                dwxdx = (-12.0_dp*rx + 24.0_dp*rx**2 - 12.0_dp*rx**3)*drdx
            else
                wx = 0.0_dp
                dwxdx = 0.0_dp
            end if

            if (ry <= 1.0_dp) then
                wy = 1.0_dp - 6.0_dp*ry**2 + 8.0_dp*ry**3 - 3.0_dp*ry**4
                dwydy = (-12.0_dp*ry + 24.0_dp*ry**2 - 12.0_dp*ry**3)*drdy
            else
                wy = 0.0_dp
                dwydy = 0.0_dp
            end if

            w(i,1) = wx*wy
            w(i,2) = wy*dwxdx
            w(i,3) = wx*dwydy
#if DIM==3
            difz = dif(3,i)
            if (dabs(difz) .le. ep) then
                drdz = 0.0_dp
            else
                drdz = (difz/dabs(difz))/ds(3,i)
            end if
            rz = dabs(dif(3,i))/ds(3,i)
            if (rz <= 1.0_dp) then
                wz = 1.0_dp - 6.0_dp*rz**2 + 8.0_dp*rz**3 - 3.0_dp*rz**4
                dwzdz = (-12.0_dp*rz + 24.0_dp*rz**2 - 12.0_dp*rz**3)*drdz
            else
                wz = 0.0_dp
                dwzdz = 0.0_dp
            end if
            w(i,1) = wx*wy*wz
            w(i,4) = wz*dwzdz
#endif
        end do
        
    end subroutine Weight_W2
    !===============================================================================================

    !===============================================================================================
    subroutine Weight_W3(dif, ds, w)

        ! Exponential weight function. 
        ! Routine adapted from "An introduction to meshfree methods and their 
        ! programming" by Liu and Gu, chapter 3.

        ! In/Out variables
        real(dp), intent(in ) :: dif(Ndim, Ne), ds(Ndim, Ne)
        real(dp), intent(out) :: w(Ne,m)

        ! Local variables
        integer  :: i
        real(dp) :: ep, difx, dify, drdx, drdy, rx, ry, wx, wy, dwxdx, dwydy, alpha
#if DIM==3
        real(dp) :: difz, drdz, rz, wz, dwzdz
#endif
        
        ep = 1.0e-20_dp
        alpha = 0.3_dp

        do i = 1,Ne
            difx = dif(1,i)
            dify = dif(2,i)
            if (dabs(difx) .le. ep) then
                drdx = 0.0_dp
            else
                drdx = (difx/dabs(difx))/ds(1,i)
            end if
            if (dabs(dify).le.ep) then
                drdy = 0.0_dp
            else
                drdy = (dify/dabs(dify))/ds(2,i)
            end if

            rx = dabs(dif(1,i))/ds(1,i)
            ry = dabs(dif(2,i))/ds(2,i)

            if (rx <= 1.0_dp) then
                wx = exp(-rx/alpha)**2
                dwxdx = -2.0_dp*rx*wx*drdx/alpha**2
            else
                wx = 0.0_dp
                dwxdx = 0.0_dp
            end if

            if (ry <= 1.0_dp) then
                wy = exp(-ry/alpha)**2
                dwydy = -2.0_dp*ry*wy*drdy/alpha**2
            else
                wy = 0.0_dp
                dwydy = 0.0_dp
            end if

            w(i,1) = wx*wy
            w(i,2) = wy*dwxdx
            w(i,3) = wx*dwydy
#if DIM==3
            difz = dif(3,i)
            if (dabs(difz) .le. ep) then
                drdz = 0.0_dp
            else
                drdz = (difz/dabs(difz))/ds(3,i)
            end if
            rz = dabs(dif(3,i))/ds(3,i)
            if (rz <= 1.0_dp) then
                wz = exp(-rz/alpha)**2
                dwzdz = -2.0_dp*rz*wz*drdz/alpha**2
            else
                wz = 0.0_dp
                dwzdz = 0.0_dp
            end if
            w(i,1) = wx*wy*wz
            w(i,4) = wz*dwzdz
#endif
        end do
        
    end subroutine Weight_W3
    !===============================================================================================

    !===============================================================================================
    function solve_wbs(u) result(x) ! solve with backward substitution

        ! This routine has been taken from Rosetta Code @
        ! https://rosettacode.org/wiki/Gaussian_elimination#Fortran
        
        real(dp)              :: u(:,:)
        integer               :: i,n
        real(dp), allocatable :: x(:)

        n = size(u,1)
        allocate(x(n))
        forall (i=n:1:-1) x(i) = ( u(i,n+1) - sum(u(i,i+1:n)*x(i+1:n)) )/u(i,i)

    end function solve_wbs
    !===============================================================================================
    
    !===============================================================================================
    function ge_wpp(a,b) result(u) ! gaussian eliminate with partial pivoting

        ! This routine has been taken from Rosetta Code @
        ! https://rosettacode.org/wiki/Gaussian_elimination#Fortran
        
        real(dp)              :: a(:,:),b(:),upi
        integer               :: i,j,n,p
        real(dp), allocatable :: u(:,:)

        n = size(a,1)
        u = reshape( [a,b], [n,n+1] )

        do j = 1,n
            ! maxloc returns indices between (1,n-j+1)
            p = maxloc(abs(u(j:n,j)),1) + j-1 
            if (p /= j) u([p,j],j) = u([j,p],j)
            u(j+1:,j) = u(j+1:,j)/u(j,j)
            do i = j+1,n+1
                upi = u(p,i)
                if (p /= j) u([p,j],i) = u([j,p],i)
                u(j+1:n,i) = u(j+1:n,i) - upi*u(j+1:n,j)
            end do
        end do

    end function ge_wpp
    !===============================================================================================

end module
