module volume_of_fluid

   ! This module contains all the procedures for the solution of the MTHINC method by
   ! Ii et al JCP 2012.

   use precision , only : dp
   use class_Grid, only : base_grid
   use class_Scalar
   use class_Vector

   implicit none

   ! Declare all fields requested by the VoF solver:
   ! vof is the volume of fluid function
   ! h is the hyperbolic tangent
   ! d is the normalization parameter
   ! curv is the curvature field
   ! norm is the local norm to the interface
   ! l is the curvature tensor
   type(scalar) :: vof, h, d, curv
   type(vector) :: norm, l

   ! Sharpness parameter
   real(dp) :: beta = 1.0_dp

   ! Switch for the quadratic reconstruction
   logical :: quadratic = .true.

   ! Switch for split advection order
   logical :: x_first = .true.

   ! Two points Gaussian quadrature in the interval [0,1]
   real(dp), parameter :: rp = 0.5_dp*(1.0_dp + 1.0_dp/sqrt(3.0_dp))
   real(dp), parameter :: rm = 0.5_dp*(1.0_dp - 1.0_dp/sqrt(3.0_dp))

   ! Cut parameter for the reconstruction
   real(dp) :: cut = 1.0d-8

   ! Coefficients for the reconstruction
   real(dp) :: cx, cy, a10, a01, a20, a02

   interface
      function distance_from_interface(x,y) result(d)
         use precision, only : dp
         real(dp), intent(in) :: x, y
         real(dp) :: d
      end function distance_from_interface
   end interface
   procedure(distance_from_interface), pointer :: distance => Null()

contains

   !========================================================================================
   subroutine allocate_vof_fields

      ! Allocate memory requested for the volume of fluid solver.
      use io

      call vof%allocate(1)
      call h%allocate(1)
      call d%allocate(1)
      call curv%allocate(1)
      call norm%allocate(1)
      call l%allocate(1)

      ! Based on the boundary conditions on the physical domain select boundary conditions
      ! for all fields
      ! Left boundary
      if (base_grid%boundary_conditions(1)%s == 'Periodic') then
         vof%bc%type_l = 0
         h%bc%type_l = 0
         d%bc%type_l = 0
         curv%bc%type_l = 0
         norm%x%bc%type_l = 0
         norm%y%bc%type_l = 0
         l%x%bc%type_l = 0
         l%y%bc%type_l = 0
#if DIM==3
         norm%z%bc%type_l = 0
         l%z%bc%type_l = 0
#endif
      elseif (base_grid%boundary_conditions(1)%s == 'Wall') then
         vof%bc%type_l = 2
         h%bc%type_l = 2
         d%bc%type_l = 2
         curv%bc%type_l = 2
         norm%x%bc%type_l = 2
         norm%y%bc%type_l = 2
         l%x%bc%type_l = 2
         l%y%bc%type_l = 2
#if DIM==3
         norm%z%bc%type_l = 2
         l%z%bc%type_l = 2
#endif
      else
         call print_error_message('ERROR: wrong bc on left boundary')
      endif

      ! Right boundary
      if (base_grid%boundary_conditions(2)%s == 'Periodic') then
         vof%bc%type_r = 0
         h%bc%type_r = 0
         d%bc%type_r = 0
         curv%bc%type_r = 0
         norm%x%bc%type_r = 0
         norm%y%bc%type_r = 0
         l%x%bc%type_r = 0
         l%y%bc%type_r = 0
#if DIM==3
         norm%z%bc%type_r = 0
         l%z%bc%type_r = 0
#endif
      elseif (base_grid%boundary_conditions(2)%s == 'Wall') then
         vof%bc%type_r = 2
         h%bc%type_r = 2
         d%bc%type_r = 2
         curv%bc%type_r = 2
         norm%x%bc%type_r = 2
         norm%y%bc%type_r = 2
         l%x%bc%type_r = 2
         l%y%bc%type_r = 2
#if DIM==3
         norm%z%bc%type_r = 2
         l%z%bc%type_r = 2
#endif
      else
         call print_error_message('ERROR: wrong bc on right boundary')
      endif

      ! Bottom boundary
      if (base_grid%boundary_conditions(3)%s == 'Periodic') then
         vof%bc%type_b = 0
         h%bc%type_b = 0
         d%bc%type_b = 0
         curv%bc%type_b = 0
         norm%x%bc%type_b = 0
         norm%y%bc%type_b = 0
         l%x%bc%type_b = 0
         l%y%bc%type_b = 0
#if DIM==3
         norm%z%bc%type_b = 0
         l%z%bc%type_b = 0
#endif
      elseif (base_grid%boundary_conditions(3)%s == 'Wall') then
         vof%bc%type_b = 2
         h%bc%type_b = 2
         d%bc%type_b = 2
         curv%bc%type_b = 2
         norm%x%bc%type_b = 2
         norm%y%bc%type_b = 2
         l%x%bc%type_b = 2
         l%y%bc%type_b = 2
#if DIM==3
         norm%z%bc%type_b = 2
         l%z%bc%type_b = 2
#endif
      else
         call print_error_message('ERROR: wrong bc on bottom boundary')
      endif

      ! Top boundary
      if (base_grid%boundary_conditions(4)%s == 'Periodic') then
         vof%bc%type_t = 0
         h%bc%type_t = 0
         d%bc%type_t = 0
         curv%bc%type_t = 0
         norm%x%bc%type_t = 0
         norm%y%bc%type_t = 0
         l%x%bc%type_t = 0
         l%y%bc%type_t = 0
#if DIM==3
         norm%z%bc%type_t = 0
         l%z%bc%type_t = 0
#endif
      elseif (base_grid%boundary_conditions(4)%s == 'Wall') then
         vof%bc%type_t = 2
         h%bc%type_t = 2
         d%bc%type_t = 2
         curv%bc%type_t = 2
         norm%x%bc%type_t = 2
         norm%y%bc%type_t = 2
         l%x%bc%type_t = 2
         l%y%bc%type_t = 2
#if DIM==3
         norm%z%bc%type_t = 2
         l%z%bc%type_t = 2
#endif
      else
         call print_error_message('ERROR: wrong bc on top boundary')
      endif

   end subroutine allocate_vof_fields
   !========================================================================================

   !========================================================================================
   subroutine get_h_from_vof()

      ! This subrouinte compute the hyperbolic tangent h from the vof field

      implicit none

      ! Local variables
      integer  :: i, j, k
      real(dp) :: A, Bp, Bm, Q, aa, bb, cc

      ! In order to reconstruct the hyperbolic tangent h it is necessary to compute the
      ! local norm to the interface
      call compute_norm

      do k = base_grid%lo(3),base_grid%hi(3)
         do j = base_grid%lo(2),base_grid%hi(2)
            do i = base_grid%lo(1),base_grid%hi(1)

               ! Reconstruct only for values larger than the cut value
               if (vof%f(i,j,k) <= cut .or. vof%f(i,j,k) >= (1.0_dp - cut)) then

                  h%f(i,j,k) = vof%f(i,j,k)
                  d%f(i,j,k) = 0.0_dp

               else
                  ! Compute coefficients based on the maximum norm
                  if (abs(norm%x%f(i,j,k)) == max(abs(norm%x%f(i,j,k)),abs(norm%y%f(i,j,k)))) then
                     cx = 0.0_dp
                     cy = 1.0_dp
                  else
                     cx = 1.0_dp
                     cy = 0.0_dp
                  endif

                  ! Eq 12 of Ii et al JCP 2012
                  a10 = norm%x%f(i,j,k) - 0.5_dp*cx*l%x%f(i,j,k)
                  a01 = norm%y%f(i,j,k) - 0.5_dp*cy*l%y%f(i,j,k)
                  a20 = 0.5_dp*cx*l%x%f(i,j,k)
                  a02 = 0.5_dp*cy*l%y%f(i,j,k)

                  A = (1.0_dp - cx)*exp(2.0_dp*beta*a10) + &
                      (1.0_dp - cy)*exp(2.0_dp*beta*a01)

                  Bp = (1.0_dp - cx)*exp(2.0_dp*beta*P(0.0_dp,rp)) + &
                       (1.0_dp - cy)*exp(2.0_dp*beta*P(rp,0.0_dp))

                  Bm = (1.0_dp - cx)*exp(2.0_dp*beta*P(0.0_dp,rm)) + &
                       (1.0_dp - cy)*exp(2.0_dp*beta*P(rm,0.0_dp))

                  Q = (1.0_dp - cx)*exp(2.0_dp*beta*a10*(2.0_dp*vof%f(i,j,k) - 1.0_dp)) + &
                      (1.0_dp - cy)*exp(2.0_dp*beta*a01*(2.0_dp*vof%f(i,j,k) - 1.0_dp))

                  aa = A*Bm*Bp*(A - Q)
                  bb = A*(Bp + Bm)*(1.0_dp - Q)
                  cc = 1.0_dp - A*Q

                  ! Compute the normalization parameter
                  d%f(i,j,k) = log(solve_quadratic(aa,bb,cc))/(2.0_dp*beta)

                  ! Eq 7 of Ii et al JCP 2012
                  h%f(i,j,k) = 0.5_dp*(1.0_dp + dtanh(beta*(P(0.5_dp,0.5_dp) + d%f(i,j,k))))

               endif

            end do
         end do
      end do

      ! Update boundary conditions for h and d
      call h%apply_bc()
      call d%apply_bc()

   end subroutine get_h_from_vof
   !========================================================================================

   !========================================================================================
   subroutine compute_norm()

      ! Compute the normal vector to the interface using the Young method

      use constants, only : small
      use class_Scalar
      
      implicit none

      ! Local variables
      integer  :: i, j, ip, jp, im, jm, k, c
      real(dp) :: mx(4), mxc, my(4), myc, normx(4), normy(4), delta, idelta, idelta2
      
      delta = base_grid%delta
      idelta = 1.0_dp/delta
      idelta2 = 1.0_dp/delta**2
      
      do k = base_grid%lo(3),base_grid%hi(3)
         do j = base_grid%lo(2),base_grid%hi(2)
            jp = j + 1
            jm = j - 1
            do i = base_grid%lo(1),base_grid%hi(1)
               ip = i + 1
               im = i - 1

               ! X component
               ! i-1/2,j-1/2
               mx(1) = 0.5_dp*(vof%f(i,jm,k) + vof%f(i,j,k) - vof%f(im,jm,k) - vof%f(im,j,k))*idelta

               ! i-1/2,j+1/2
               mx(2) = 0.5_dp*(vof%f(i,j,k) + vof%f(i,jp,k) - vof%f(im,j,k) - vof%f(im,jp,k))*idelta

               ! i+1/2,j+1/2
               mx(3) = 0.5_dp*(vof%f(ip,j,k) + vof%f(ip,jp,k) - vof%f(i,j,k) - vof%f(i,jp,k))*idelta

               ! i+1/2,j-1/2
               mx(4) = 0.5_dp*(vof%f(ip,jm,k) + vof%f(ip,j,k) - vof%f(i,jm,k) - vof%f(i,j,k))*idelta

               ! Cell center value
               mxc = 0.25_dp*(mx(1) + mx(2) + mx(3) + mx(4))

               ! Y component
               ! i-1/2,j-1/2
               my(1) = 0.5_dp*(vof%f(im,j,k) + vof%f(i,j,k) - vof%f(im,jm,k) - vof%f(i,jm,k))*idelta

               ! i-1/2,j+1/2
               my(2) = 0.5_dp*(vof%f(im,jp,k) + vof%f(i,jp,k) - vof%f(im,j,k) - vof%f(i,j,k))*idelta

               ! i+1/2,j+1/2
               my(3) = 0.5_dp*(vof%f(i,jp,k) + vof%f(ip,jp,k) - vof%f(i,j,k) - vof%f(ip,j,k))*idelta

               ! i+1/2,j-1/2
               my(4) = 0.5_dp*(vof%f(i,j,k) + vof%f(ip,j,k) - vof%f(i,jm,k) - vof%f(ip,jm,k))*idelta

               ! Cell center value
               myc = 0.25_dp*(my(1) + my(2) + my(3) + my(4))

               ! Norm vector components
               do c = 1,4
                  normx(c) = mx(c)/sqrt(mx(c)**2 + my(c)**2 + small)
                  normy(c) = my(c)/sqrt(mx(c)**2 + my(c)**2 + small)
               end do
               norm%x%f(i,j,k) = mxc/sqrt(mxc**2 + myc**2 + small)
               norm%y%f(i,j,k) = myc/sqrt(mxc**2 + myc**2 + small)

               ! Curvature components
               if (quadratic) then
                  l%x%f(i,j,k) = 0.5_dp*delta*(normx(4) + normx(3) - normx(2) - normx(1))
                  l%y%f(i,j,k) = 0.5_dp*delta*(normy(2) + normy(3) - normy(1) - normy(4))
                  curv%f(i,j,k) = -(l%x%f(i,j,k) + l%y%f(i,j,k))*idelta2
               else
                  l%x%f(i,j,k) = 0.0_dp
                  l%y%f(i,j,k) = 0.0_dp
                  curv%f(i,j,k) = -(l%x%f(i,j,k) + l%y%f(i,j,k))*idelta2
               endif
               
            end do
         end do
      end do

      ! Update halo and bc
      call curv%apply_bc()
      call norm%apply_bc()
      call l%apply_bc()
      
   end subroutine compute_norm
   !========================================================================================

   !========================================================================================
   real(dp) pure function P(x, y)

      ! Quadratic surface function.

      implicit none

      real(dp), intent(in) :: x, y

      ! Eq. 9 of Ii et al JCP 2012
      P = cx*a20*x**2 + cy*a02*y**2 + a10*x + a01*y

   end function P
   !========================================================================================

   !========================================================================================
   pure function solve_quadratic(a, b, c) result(x)

      ! Find root of quadratic equation

      implicit none

      real(dp), intent(in) :: a, b, c
      real(dp) :: x1, x2, x

      x1 = (-b + dsqrt(b**2 - 4.0_dp*a*c))/(2.0_dp*a)
      x2 = (-b - dsqrt(b**2 - 4.0_dp*a*c))/(2.0_dp*a)

      x = max(x1,x2)

   end function solve_quadratic
   !========================================================================================

   !========================================================================================
   subroutine advect_vof(v, dt)

      ! This subroutine solve the advection equation for the volume of fluid function
      ! with a given velocity fields v

      implicit none

      ! In/Out variables
      real(dp)    , intent(in) :: dt
      type(vector), intent(in) :: v

      ! Local variables
      integer  :: i, j, k, lo(3), hi(3)
      real(dp) :: fp, fm, delta
      type(scalar) :: vof1, vof2

      call vof1%allocate(1)
      call vof2%allocate(1)

      lo = base_grid%lo
      hi = base_grid%hi
      delta = base_grid%delta

      ! First reconstruct the hyperbolic tangent function
      call get_h_from_vof

      if (x_first) then

         ! X direction
         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)
                  fp = compute_flux(1, i  , j, k, v%x%f(i  ,j,k), dt, delta)
                  fm = compute_flux(1, i-1, j, k, v%x%f(i-1,j,k), dt, delta)
                  vof1%f(i,j,k) = (vof%f(i,j,k) - (fp - fm)/delta)/ &
                                  (1.0_dp - dt*(v%x%f(i,j,k) - v%x%f(i-1,j,k))/delta)
               end do
            end do
         end do

         ! Get h from vof1
         vof = vof1
         call vof%apply_bc()
         call get_h_from_vof

         ! y direction
         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)
                  fp = compute_flux(2, i, j  , k, v%y%f(i,j  ,k), dt, delta)
                  fm = compute_flux(2, i, j-1, k, v%y%f(i,j-1,k), dt, delta)
                  vof2%f(i,j,k) = (vof1%f(i,j,k) - (fp - fm)/delta)/ &
                                   (1.0_dp - dt*(v%y%f(i,j,k) - v%y%f(i,j-1,k))/delta)
               end do
            end do
         end do

         ! Compute vof at time n+1
         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)
                  vof%f(i,j,k) = vof2%f(i,j,k) - dt*(                        &
                     vof1%f(i,j,k)*(v%x%f(i,j,k) - v%x%f(i-1,j,k))/delta + &
                     vof2%f(i,j,k)*(v%y%f(i,j,k) - v%y%f(i,j-1,k))/delta)
               end do
            end do
         end do

         x_first = .false.

      else

         ! y direction
         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)
                  fp = compute_flux(2, i, j  , k, v%y%f(i,j  ,k), dt, delta)
                  fm = compute_flux(2, i, j-1, k, v%y%f(i,j-1,k), dt, delta)
                  vof1%f(i,j,k) = (vof%f(i,j,k) - (fp - fm)/delta)/ &
                     (1.0_dp - dt*(v%y%f(i,j,k) - v%y%f(i,j-1,k))/delta)
               end do
            end do
         end do

         ! Get h from vof1
         vof = vof1
         call vof%apply_bc()
         call get_h_from_vof

         ! X direction
         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)
                  fp = compute_flux(1, i  , j, k, v%x%f(i  ,j,k), dt, delta)
                  fm = compute_flux(1, i-1, j, k, v%x%f(i-1,j,k), dt, delta)
                  vof2%f(i,j,k) = (vof1%f(i,j,k) - (fp - fm)/delta)/ &
                     (1.0_dp - dt*(v%x%f(i,j,k) - v%x%f(i-1,j,k))/delta)
               end do
            end do
         end do

         ! Compute vof at time n+1
         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)
                  vof%f(i,j,k) = vof2%f(i,j,k) - dt*(                        &
                     vof2%f(i,j,k)*(v%x%f(i,j,k) - v%x%f(i-1,j,k))/delta + &
                     vof1%f(i,j,k)*(v%y%f(i,j,k) - v%y%f(i,j-1,k))/delta)
               end do
            end do
         end do
      
         x_first = .true.

      endif


      ! Update halo and bc
      call vof%apply_bc()

      ! Free the memory for the temporary vof fields
      call vof1%destroy()
      call vof2%destroy()

   end subroutine advect_vof
   !========================================================================================

   !========================================================================================
   function compute_flux(dir, i, j, k, u, dt, delta) result(f)

      implicit none

      ! In/Out variables
      integer , intent(in) :: dir, i, j, k
      real(dp), intent(in) :: u, dt, delta

      ! Local variable
      integer  :: sgn, ii, jj, kk
      real(dp) :: f, xa, xb, ya, yb

      ! Select integration intervals
      if (dir == 1) then
         if (u >= 0.0_dp)then
            xa = 1.0_dp - dt*u/delta
            xb = 1.0_dp
            sgn = 1
            ii = i
            jj = j
            kk = k
         else
            xa = 0.0_dp
            xb = -dt*u/delta
            sgn = -1
            ii = i + 1
            jj = j
            kk = k
         endif
         ya = 0.0_dp
         yb = 1.0_dp
      elseif (dir == 2) then
         xa = 0.0_dp
         xb = 1.0_dp
         if (u >= 0.0_dp) then
            ya = 1.0_dp - dt*u/delta
            yb = 1.0_dp
            sgn = 1
            ii = i
            jj = j
            kk = k
         else
            ya = 0.0_dp
            yb = -dt*u/delta
            sgn = -1
            ii = i
            jj = j + 1
            kk = k
         endif
      endif

      if (vof%f(ii,jj,kk) <= cut .or. vof%f(ii,jj,kk) >= (1.0_dp - cut)) then
         f = sgn*delta*vof%f(ii,jj,kk)*(xb - xa)*(yb - ya)
      else
         ! Select type of integration based on the local normal vector
         if (abs(norm%x%f(ii,jj,kk)) == max(abs(norm%x%f(ii,jj,kk)),dabs(norm%y%f(ii,jj,kk)))) then
            cx = 0.0_dp
            cy = 1.0_dp
            a10 = norm%x%f(ii,jj,kk) - 0.5_dp*cx*l%x%f(ii,jj,kk)
            a01 = norm%y%f(ii,jj,kk) - 0.5_dp*cy*l%y%f(ii,jj,kk)
            a20 = 0.5_dp*cx*l%x%f(ii,jj,kk)
            a02 = 0.5_dp*cy*l%y%f(ii,jj,kk)

            ! Analytical integration in x, numerical integration in y
            f = sgn*delta*Num_Int(ya, yb, An_Int(xa, xb, xa, xb, rm*(ya + yb), rm*(yb + ya), &
               a10, d%f(ii,jj,kk)), &
               An_Int(xa, xb, xa, xb, rp*(ya + yb), rp*(yb + ya), &
               a10, d%f(ii,jj,kk)))
         else
            cx = 1.0_dp
            cy = 0.0_dp
            a10 = norm%x%f(ii,jj,kk) - 0.5_dp*cx*l%x%f(ii,jj,kk)
            a01 = norm%y%f(ii,jj,kk) - 0.5_dp*cy*l%y%f(ii,jj,kk)
            a20 = 0.5_dp*cx*l%x%f(ii,jj,kk)
            a02 = 0.5_dp*cy*l%y%f(ii,jj,kk)

            ! Numerical integration in x, analytical integration in y
            f = sgn*delta*Num_Int(xa, xb, An_Int(ya, yb, rm*(xa + xb), rm*(xa + xb), ya, yb, &
               a01, d%f(ii,jj,kk)), &
               An_Int(ya, yb, rp*(xa + xb), rp*(xa + xb), ya, yb, &
               a01, d%f(ii,jj,kk)))
         endif
      endif

   end function compute_flux
   !========================================================================================

   !========================================================================================
   real(dp) pure function Num_Int(a, b, qrm, qrp)

      ! Numerical integration of the hyperbolic tangent h in the interval [a, b]
      ! equation 14 of Ii et al JCP 2012.

      ! In/Out variables
      real(dp), intent(in) :: a, b, qrm, qrp

      Num_Int = 0.5_dp*(qrm + qrp)*(b - a)

   end function Num_Int
   !========================================================================================

   !========================================================================================
   real(dp) pure function An_Int(a, b, xa, xb, ya, yb, c, d)

      ! Analytical integration of the hyperbolic tangent h in the interval [a, b] equation 14
      ! of Ii et al JCP 2012, evaluating the surface P in (xb,yb) and (xa,ya)

      ! In/Out variables
      real(dp), intent(in) :: a, b, xa, xb, ya, yb, c, d

      An_Int = 0.5_dp*(b - a + 1.0_dp/(c*beta)*log(cosh(beta*(P(xb,yb) + d))/&
         cosh(beta*(P(xa,ya) + d))))

   end function An_Int
   !========================================================================================

   !========================================================================================
   subroutine get_vof_from_distance()

      ! This subroutine compute the vof function from a distance function distance

      use io    , only : print_error_message

      ! Local variables
      integer  :: i, j, k
      real(dp) :: delta, y, yp, ym, x, xp, xm

      if (associated(distance) .eqv. .false.) then
         call print_error_message('ERROR: distance function not defined.')
      endif

      delta = base_grid%delta

      do k = base_grid%lo(3),base_grid%hi(3)
         do j = base_grid%lo(2),base_grid%hi(2)
            y = base_grid%y(j)
            yp = y + delta*(rp - 0.5_dp)
            ym = y + delta*(rm - 0.5_dp)
            do i = base_grid%lo(1),base_grid%hi(1)
               x = base_grid%x(i)
               xp = x + delta*(rp - 0.5_dp)
               xm = x + delta*(rm - 0.5_dp)

               vof%f(i,j,k) = 0.5_dp*(0.5_dp*(0.5_dp*(1.0_dp + tanh(beta*distance(xm,ym)/delta))  + &
                                              0.5_dp*(1.0_dp + tanh(beta*distance(xp,ym)/delta))) + &
                                      0.5_dp*(0.5_dp*(1.0_dp + tanh(beta*distance(xm,yp)/delta))  + &
                                              0.5_dp*(1.0_dp + tanh(beta*distance(xp,yp)/delta))))

               h%f(i,j,k) = 0.5_dp*(1.0_dp + tanh(beta*distance(x,y)/delta))
               
            end do
         end do
      end do

      ! Update boundary conditions for the vof
      call vof%apply_bc()
      call h%apply_bc()

   end subroutine get_vof_from_distance
   !========================================================================================

   !========================================================================================
   subroutine check_vof_integral(int_phase_1, int_phase_2)

      use mpi

      implicit none

      ! In/Out variables
      real(dp), intent(out) :: int_phase_1, int_phase_2

      ! Local variables
      integer :: i, j, k, ierror

      int_phase_1 = 0.0_dp
      int_phase_2 = 0.0_dp

      do k = base_grid%lo(3),base_grid%hi(3)
         do j = base_grid%lo(2),base_grid%hi(2)
            do i = base_grid%lo(1),base_grid%hi(1)
               int_phase_1 = int_phase_1 + vof%f(i,j,k)
               int_phase_2 = int_phase_2 + (1.0_dp - vof%f(i,j,k))
            end do
         end do
      end do

      call mpi_allreduce(mpi_in_place,int_phase_1,1,mpi_real8,mpi_sum,mpi_comm_world,ierror)
      call mpi_allreduce(mpi_in_place,int_phase_2,1,mpi_real8,mpi_sum,mpi_comm_world,ierror)

      int_phase_1 = int_phase_1*base_grid%delta*base_grid%delta*base_grid%delta
      int_phase_2 = int_phase_2*base_grid%delta*base_grid%delta*base_grid%delta

   end subroutine check_vof_integral
   !========================================================================================

   !========================================================================================
   subroutine destroy_vof

      call vof%destroy()
      call h%destroy()
      call d%destroy()
      call curv%destroy()
      call norm%destroy()
      call l%destroy()

   end subroutine destroy_vof
   !========================================================================================

end module volume_of_fluid
