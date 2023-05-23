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

contains

    !========================================================================================
    function interpolate(f, xl, ie, ind) result(fl)

        ! This function interpolate the value of the scalar field f in the point xl
        use mpi
        use scalar_mod
    
        ! In/Out variables
        integer     , intent(in) :: ie(3), ind
        real(dp)    , intent(in) :: xl(Ndim)
        type(scalar), intent(in) :: f
        
        ! Local variables
        integer  :: si, sj, q, ii, jj, ierror
        real(dp) :: xe, ye, xs(Ndim,Ne), fl(m), ds(Ndim,Ne), fk(Ne), stagger(3)
#if DIM==3
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
        rank_with_point: if (ie(2) >= f%G%lo(2) .and. ie(2) <= f%G%hi(2)) then 
            
            ! Build the array of positions and f values in the support domain
            q = 1
            fk = 0.0_dp
            do sj = -1,1
                jj = ie(2) + sj
                do si = -1,1
                    ii = ie(1) + si
                    xs(1,q) = xe + si*f%G%delta
                    xs(2,q) = ye + sj*f%G%delta
                    fk(q) = f%f(ii,jj,1)
                    ds(1,q) = merge(f%G%delta, 2.0_dp*f%G%delta, si == 0)
                    ds(2,q) = merge(f%G%delta, 2.0_dp*f%G%delta, sj == 0)
                    q = q + 1
                end do
            end do
        
            ! Size of the support domain
            ds = 1.60_dp*f%G%delta
            
            ! Compute shape function phi
            call MLS_ShapeFunc_2D(xl, xs, ds)
            
            ! Compute interpolated value in (xl, yl)
            fl = 0.0_dp
            do q = 1, Ne
                fl(1) = fl(1) + phi(1,q)*fk(q)
                fl(2) = fl(2) + phi(2,q)*fk(q)
                fl(3) = fl(3) + phi(3,q)*fk(q)
            end do

        else

            ! Set to zero on other ranks
            fl = 0.0_dp

        end if rank_with_point

        ! Comunicate the interpolated value
        call mpi_allreduce(mpi_in_place,fl,m,mpi_real8,mpi_sum,mpi_comm_world,ierror)
        
    end function interpolate
    !========================================================================================

    !========================================================================================
    subroutine MLS_ShapeFunc_2D(gpos, x, ds)

        ! Compute MLS shape functions and their derivatives. Routine adapted from:
        ! "An introduction to meshfree methods and their programming" by Liu and Gu
        ! chapter 3 appendix.    
        
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

        ! Compute gamma
        gam = 0.0_dp

        c = gp(1,:)
        do i = 1,m
            do j = 1,m
                aa(i,j) = A(i,j,1)
            end do
        end do
        
        c = solve_wbs(ge_wpp(aa,c))
        gam(:,1) = c

        ! Compute dgamdx
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

        ! ************* Compute dgamdy
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
        
        ! ************* Compute Phi and their derivatives
        do i = 1,Ne
            do j = 1,m
                phi(j,i) = 0.0_dp
            enddo
            do j = 1,m
                phi(1,i) = phi(1,i) + gam(j,1)*b(j,i,1)
                phi(2,i) = phi(2,i) + gam(j,2)*b(j,i,1) + gam(j,1)*b(j,i,2)
                phi(3,i) = phi(3,i) + gam(j,3)*b(j,i,1) + gam(j,1)*b(j,i,3)
            end do
        end do
        
    end subroutine MLS_ShapeFunc_2D
    !==============================================================================================

    !==============================================================================================
    subroutine Compute_Basis(gpos, gp)

        ! Evaluate the basis vector gp in the interpolation point gpos
        ! The base vector is gp = [1.0, x, y] (in 2D), gp = [1.0, x, y, z] (in 3D).
        ! Routine adapted from:
        ! "An introduction to meshfree methods and their programming" by Liu and Gu
        ! chapter 3 appendix.

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
    !==============================================================================================

    !==============================================================================================
    subroutine Compute_AB(gpos, x, ds, A, B)

        ! Compute A matrix and B matrix and their derivatives. Routine adapted from:
        ! "An introduction to meshfree methods and their programming" by Liu and Gu
        ! chapter 3 appendix.

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
            dif(3,i) = gpos(3) - x(3,i)
#endif
        enddo

        ! Compute the weight function the support domain
        call Weight_W1(dif, ds, W)

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
    !==============================================================================================

    !==============================================================================================
    subroutine Weight_W1(dif, ds, w)

        ! Cubic spline weight function. Routine adapted from
        ! "An introduction to meshfree methods and their programming" by Liu and Gu
        ! chapter 3 appendix.

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
                wx = (4.0_dp/3.0_dp) - 4.0_dp*rx + 4.0_dp*rx*rx - (4.0_dp/3.0_dp)*rx*rx*rx
                dwxdx=(-4.0_dp + 8.0_dp*rx - 4.0_dp*rx*rx)*drdx
            else if (rx <= 0.5_dp) then
                wx = (2.0_dp/3.0_dp) - 4.0_dp*rx*rx + 4.0_dp*rx*rx*rx
                dwxdx = (-8.0_dp*rx + 12.0_dp*rx*rx)*drdx
            end if

            if (ry > 0.5_dp) then
                wy = (4.0_dp/3.0_dp) - 4.0_dp*ry + 4.0_dp*ry*ry - (4.0_dp/3.0_dp)*ry*ry*ry
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
                wz = (4.0_dp/3.0_dp) - 4.0_dp*rz + 4.0_dp*rz*rz - (4.0_dp/3.0_dp)*rz*rz*rz
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
    !==============================================================================================
    
    !==============================================================================================
    function solve_wbs(u) result(x) ! solve with backward substitution

        ! This routine has been taken from Rosetta Code @
        ! https://rosettacode.org/wiki/Gaussian_elimination#Fortran
        
        real(dp)              :: u(:,:)
        integer               :: i,n
        real(dp), allocatable :: x(:)

        n = size(u,1)
        allocate(x(n))
        forall (i=n:1:-1) x(i) = ( u(i,n+1) - sum(u(i,i+1:n)*x(i+1:n)) ) / u(i,i)

    end function solve_wbs
    !==============================================================================================
    
    !==============================================================================================
    function ge_wpp(a,b) result(u) ! gaussian eliminate with partial pivoting

        ! This routine has been taken from Rosetta Code @
        ! https://rosettacode.org/wiki/Gaussian_elimination#Fortran
        
        real(dp)              :: a(:,:),b(:),upi
        integer               :: i,j,n,p
        real(dp), allocatable :: u(:,:)

        n = size(a,1)
        u = reshape( [a,b], [n,n+1] )

        do j = 1,n
            p = maxloc(abs(u(j:n,j)),1) + j-1 ! maxloc returns indices between (1,n-j+1)
            if (p /= j) u([p,j],j) = u([j,p],j)
            u(j+1:,j) = u(j+1:,j)/u(j,j)
            do i = j+1,n+1
                upi = u(p,i)
                if (p /= j) u([p,j],i) = u([j,p],i)
                u(j+1:n,i) = u(j+1:n,i) - upi*u(j+1:n,j)
            end do
        end do

    end function ge_wpp
    !=============================================================================================

end module
