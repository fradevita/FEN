module navier_stokes

   ! This module contains all procedures for the solution of the Navier-Stokes equation.
   ! Variables are: pressure(p), velocity (v), density (rho), viscosity (mu)
   ! rhs of the equation (dv), old rhs (dv_o) for AB2, grad_p to store the pressure 
   ! gradient, S for volumetric source terms. rhof is an additional field to 
   ! interpolate density on the cell face.

   use precision   , only : dp
   use class_Scalar
   use class_Vector
   use class_Tensor

   implicit none

   ! Single fase material properties
   real(dp) :: density = 1.0_dp, viscosity = 1.0_dp

   ! Gravity vector
   real(dp) :: g(3) = [0.0_dp, 0.0_dp, 0.0_dp]

   ! Maximum divergence maximum CFL, and CFL limit
   real(dp) :: maxdiv, maxCFL, CFL = 1.0_dp

   ! Timestep limits
   real(dp) :: dt_visc, dt_conv
#ifdef MF
   real(dp) :: dt_surf, CFL_rho
#endif
   
   ! Old timestep for constant CFL simulations
   real(dp) :: dt_o

   ! Navier-Stokes fields
   type(scalar) :: p, phi, rho, mu
   type(vector) :: v, dv, dv_o, grad_p, S, rhof
   type(tensor) :: D

   ! Flag to select between constant CFL or constant dt simulation
   ! by default uses constant dt
   logical :: constant_CFL = .false.
   
   ! Flag for constant or variable viscosity, ture by default
   logical :: constant_viscosity = .true.

contains

   !========================================================================================
   subroutine navier_stokes_solver(step, dt)

      ! The time marching scheme is a classical fractional step. The steps are three:
      ! 1): compute the predicted velocity field
      ! 2): solve Poisson equation for the pressure increment
      ! 3): update velocity and pressure.

      use Poisson        , only : solve_Poisson
      use fields         , only : divergence
      use class_Grid     , only : base_grid
#ifdef MF
      use volume_of_fluid, only : vof
      use multiphase     , only : advect_interface, update_material_properties
      use multiphase     , only : p_hat, p_o, rhomin
#endif
#ifdef IBM
      use ibm
      use eulerian_ibm
#endif

      ! In/Out variables
      integer , intent(in   ) :: step
      real(dp), intent(inout) :: dt

      ! Local variables
      integer :: i, j, k

      ! Check if the simulation is at constant CFL or constant timestep and
      ! in case update dt
      if (constant_CFL) call update_timestep(dt)

#ifdef MF
      ! Advect the fluid-fluid interface with the velocity field v over a timestep dt
      call advect_interface(v, dt)

      ! Update material properties based on the new location of the interface
      call update_material_properties(rho, mu, vof)

      ! Compute approximated pressure p_hat following Dodd & Ferrante JCP 2014
      if (constant_CFL) then
         p_hat%f = p_o%f + (dt + dt_o)*(p%f - p_o%f)/dt_o
      else
         p_hat%f = 2.0_dp*p%f - p_o%f
      endif
      
      ! Update boundary condtiions for the approximate pressure
      call p_hat%apply_bc()
#endif

! #ifdef FSI
!       ! Pressure extrapolation inside the solid phase is necessary only for FSI problems
!       ! when using the eulerian ibm.
!       call pressure_extrapolation(p, rho, g, eulerian_solid_list)
! #endif

      ! Compute the predicted velocity field
      call predicted_velocity_field(dt)

      ! Solve Poisson equation
      call divergence(v, phi)
#ifdef MF
      phi%f = phi%f*rhomin/dt
#else
      do k = base_grid%lo(3),base_grid%hi(3)
         do j = base_grid%lo(2),base_grid%hi(2)
            do i = base_grid%lo(1),base_grid%hi(1)
               phi%f(i,j,k) = phi%f(i,j,k)*rho%f(i,j,k)/dt
            end do
         end do
      end do
#endif
      call solve_Poisson(phi)
      call phi%apply_bc

      ! Correct the velocity field
      call correct_velocity_field(dt)

      ! Update the pressure
      call update_pressure

      ! Performs two checks: compute the maximum divergence of the new velocity field
      ! and the current maximum CFL
      call checks(dt, step)

   end subroutine navier_stokes_solver
   !========================================================================================

   !========================================================================================
   subroutine predicted_velocity_field(dt)

      use class_Grid, only : base_grid
      use fields    , only : gradient, center_to_face
#ifdef MF
      use multiphase, only : p_hat, grad_p_hat, irhomin
      use class_Grid, only : base_grid
#endif
#ifdef IBM
      use ibm       , only : compute_ibm_forcing, Fe
#endif 

      ! In/Out variables
      real(dp), intent(in) :: dt

      ! Local variables
      integer      :: i, j, k
      real(dp)     :: A, B
      type(vector) :: RHS

      ! Coefficient for variable timestep Adams-Bashfort 2nd Order time marching scheme
      A = 1.0_dp + 0.5_dp*dt/dt_o
      B = -0.5_dp*dt/dt_o

      ! Evaluate density on the cell face
      call center_to_face(rho, rhof)

      ! Compute RHS of momentum equations
      call compute_explicit_terms(dv)
      call gradient(p, grad_p)

      ! Assemble RHS
      call RHS%allocate()
      RHS%x%f = -grad_p%x%f/rhof%x%f + A*dv%x%f + B*dv_o%x%f + g(1)
      RHS%y%f = -grad_p%y%f/rhof%y%f + A*dv%y%f + B*dv_o%y%f + g(2)
#if DIM==3
      RHS%z%f = -grad_p%z%f/rhof%z%f + A*dv%z%f + B*dv_o%z%f + g(3)
#endif
#ifdef MF
      call gradient(p_hat, grad_p_hat)
      RHS%x%f = RHS%x%f + grad_p%x%f/rhof%x%f -irhomin*grad_p%x%f - &
               (1.0_dp/rhof%x%f - irhomin)*grad_p_hat%x%f
      RHS%y%f = RHS%y%f + grad_p%y%f/rhof%y%f -irhomin*grad_p%y%f - &
               (1.0_dp/rhof%y%f - irhomin)*grad_p_hat%y%f
#if DIM==3
      RHS%z%f = RHS%z%f + grad_p%z%f/rhof%z%f -irhomin*grad_p%z%f - &
               (1.0_dp/rhof%z%f - irhomin)*grad_p_hat%z%f
#endif
#endif

#ifdef IBM
      ! If using the Immersed Boundary Method evaluate the forcing vector field
      call compute_ibm_forcing(v, RHS, dt, Fe)

      ! and add it to the RHS of momentum
      RHS%x%f = RHS%x%f + Fe%x%f
      RHS%y%f = RHS%y%f + Fe%y%f
#if DIM==3
      RHS%z%f = RHS%z%f + Fe%z%f
#endif
#endif

      ! Compute intermediate velocity field
      do k = base_grid%lo(3),base_grid%hi(3)
         do j = base_grid%lo(2),base_grid%hi(2)
            do i = base_grid%lo(1),base_grid%hi(1)

               v%x%f(i,j,k) = v%x%f(i,j,k) + dt*RHS%x%f(i,j,k)
               v%y%f(i,j,k) = v%y%f(i,j,k) + dt*RHS%y%f(i,j,k)
#if DIM==3
               v%z%f(i,j,k) = v%z%f(i,j,k) + dt*RHS%z%f(i,j,k)
#endif
            enddo
         end do
      end do

      ! Save old RHS
      dv_o%x%f = dv%x%f
      dv_o%y%f = dv%y%f
#if DIM==3
      dv_o%z%f = dv%z%f
#endif

      ! Apply phisical BC on the predicted velocity field
      call v%apply_bc()

      ! Free RHS memory
      call RHS%destroy()

   end subroutine predicted_velocity_field
   !========================================================================================

   !========================================================================================
   subroutine compute_explicit_terms(RHS)

      use class_Grid, only : base_grid

      ! In/Out variables
      type(vector), intent(inout) :: RHS

      ! Local variables
      integer :: i, j, k

      ! Set RHS to zero
      RHS%x%f = 0.0_dp
      RHS%y%f = 0.0_dp
#if DIM==3
      RHS%z%f = 0.0_dp
#endif

      ! Add convective terms
      call add_advection(RHS)

      ! Add viscous terms
      call add_diffusion(RHS)

      ! For multiphase flow add the surface tension term
#ifdef MF
      call add_surface_tension(RHS)
#endif

      ! Add body forces
      do k = base_grid%lo(3),base_grid%hi(3)
         do j = base_grid%lo(2),base_grid%hi(2)
            do i = base_grid%lo(1),base_grid%hi(1)
               RHS%x%f(i,j,k) = RHS%x%f(i,j,k) + S%x%f(i,j,k)/rhof%x%f(i,j,k)
               RHS%y%f(i,j,k) = RHS%y%f(i,j,k) + S%y%f(i,j,k)/rhof%y%f(i,j,k)
#if DIM==3
               RHS%z%f(i,j,k) = RHS%z%f(i,j,k) + S%z%f(i,j,k)/rhof%z%f(i,j,k)
#endif
            end do
         end do
      end do

   end subroutine compute_explicit_terms
   !========================================================================================

   !========================================================================================
   subroutine add_advection(RHS)

      ! Convective terms computation.
      use class_Grid, only : base_grid

      ! In/Out variables
      type(vector), intent(inout) :: RHS

      ! Local variables
      integer  :: i, j, k, ip, im, jp, jm, lo(3), hi(3)
      real(dp) :: idelta
      real(dp) :: uuip, uuim, uvjp, uvjm
      real(dp) :: vvjp, vvjm, vuip, vuim
#if DIM==3
      integer  :: kp, km
      real(dp) :: uwkp, uwkm, vwkp, vwkm, wuip, wuim, wvjp, wvjm, wwkp, wwkm
#endif

      lo = base_grid%lo
      hi = base_grid%hi
      idelta  = 1.0_dp/base_grid%delta

      do k = lo(3),hi(3)
#if DIM==3
         kp = k + 1
         km = k - 1
#endif
         do j = lo(2),hi(2)
            jp = j + 1
            jm = j - 1
            do i = lo(1),hi(1)
               ip = i + 1
               im = i - 1
               !==== u terms ================================================================
               ! Advection: -d(uu)/dx - d(uv)/dy - d(uw)/dz
               ! d(uu)/dx = (uuip - uuim)/dx
               uuip = 0.25_dp*(v%x%f(ip,j,k) + v%x%f(i,j,k))**2
               uuim = 0.25_dp*(v%x%f(im,j,k) + v%x%f(i,j,k))**2

               ! d(uv)/dy = (uvjp - uvjm)/dy
               uvjp = (v%x%f(i,jp,k) + v%x%f(i,j ,k))*(v%y%f(ip,j ,k) + v%y%f(i,j ,k))*0.25_dp
               uvjm = (v%x%f(i,j ,k) + v%x%f(i,jm,k))*(v%y%f(ip,jm,k) + v%y%f(i,jm,k))*0.25_dp

               ! Add advection terms in x to dv
               RHS%x%f(i,j,k) = RHS%x%f(i,j,k) - (uuip - uuim)*idelta - (uvjp - uvjm)*idelta
#if DIM==3
               ! d(uw)/dz = (uwkp - uwkm)/dz
               uwkp = (v%x%f(i,j,kp) + v%x%f(i,j ,k))*(v%z%f(ip,j ,k) + v%z%f(i,j ,k))*0.25_dp
               uwkm = (v%x%f(i,j ,k) + v%x%f(i,j,km))*(v%z%f(ip,j,km) + v%z%f(i,j,km))*0.25_dp
               RHS%x%f(i,j,k) = RHS%x%f(i,j,k) - (uwkp - uwkm)*idelta
#endif
               !==== v terms ===================================================================
               ! Advection: -d(vu)/dx - d(vv)/dy - d(vw)/dz
               ! d(vu)/dx = (vuip - vuim)/delta
               vuip = (v%y%f(ip,j,k) + v%y%f(i ,j,k))*(v%x%f(i ,jp,k) + v%x%f(i ,j,k))*0.25_dp
               vuim = (v%y%f(i ,j,k) + v%y%f(im,j,k))*(v%x%f(im,jp,k) + v%x%f(im,j,k))*0.25_dp

               ! d(vv)/dy = (vvjp - vvjm)/delta
               vvjp = 0.25_dp*(v%y%f(i,jp,k) + v%y%f(i,j,k))**2
               vvjm = 0.25_dp*(v%y%f(i,jm,k) + v%y%f(i,j,k))**2

               ! Add advection terms in y to dv
               RHS%y%f(i,j,k) = RHS%y%f(i,j,k) - (vuip - vuim)*idelta - (vvjp - vvjm)*idelta
#if DIM==3
               vwkp = (v%y%f(i,j,kp) + v%y%f(i,j,k ))*(v%z%f(i,jp,k ) + v%z%f(i,j,k ))*0.25_dp
               vwkm = (v%y%f(i,j,k ) + v%y%f(i,j,km))*(v%z%f(i,jp,km) + v%z%f(i,j,km))*0.25_dp
               RHS%y%f(i,j,k) = RHS%y%f(i,j,k) - (vwkp - vwkm)*idelta
#endif

#if DIM==3
               !==== w terms ===================================================================
               ! Advection: -d(wu)/dx - d(wv)/dy - d(ww)/dz
               ! d(wu)/dx = (wuip - wuim)/delta
               wuip = (v%z%f(i,j,k) + v%z%f(ip,j,k))*(v%x%f(i ,j,k) + v%x%f(i ,j,kp))*0.25_dp
               wuim = (v%z%f(i,j,k) + v%z%f(im,j,k))*(v%x%f(im,j,k) + v%x%f(im,j,kp))*0.25_dp

               ! d(wv)/dy = (wvjp - wvjm)/delta
               wvjp = (v%z%f(i,j,k) + v%z%f(i,jp,k))*(v%y%f(i ,j,k) + v%y%f(i,j ,kp))*0.25_dp
               wvjm = (v%z%f(i,j,k) + v%z%f(i,jm,k))*(v%y%f(i,jm,k) + v%y%f(i,jm,kp))*0.25_dp

               ! d(ww)/dz = (wwkp - wwkm)/delta
               wwkp = (v%z%f(i,j,k) + v%z%f(i,j,kp))*(v%z%f(i ,j,k) + v%z%f(i,j,kp ))*0.25_dp
               wwkm = (v%z%f(i,j,k) + v%z%f(i,j,km))*(v%z%f(i ,j,k) + v%z%f(i,j,km ))*0.25_dp

               ! Add advection terms in z to dv
               RHS%z%f(i,j,k) = RHS%z%f(i,j,k) - (wuip - wuim)*idelta - (wvjp - wvjm)*idelta - &
                  (wwkp - wwkm)*idelta
#endif
            end do
         end do
      end do

   end subroutine add_advection
   !========================================================================================

   !========================================================================================
   subroutine add_diffusion(RHS)

      ! Add diffusion term to the RHS of Navier-Stokes equations. Diffusion term is
      ! div(2 mu D) = 2 D grad(mu) + 2 mu div(D), D being the simmetric part of the velocity
      ! gradient. If viscosity is constant it is equal to = mu laplacian(v).
      ! We write it as div(2 mu D) = 2 D grad(mu) + mu laplacian(v)

      use class_Grid , only : base_grid
      use fields, only : laplacian
#ifdef MF
      use fields, only : gradient
#endif
#ifdef NN
      use fields       , only : gradient
      use non_newtonian, only : compute_viscosity
#endif

      ! Local variables
      type(vector), intent(inout) :: RHS

      ! Local variables
      integer  :: i, j, k, ip, im, jp, jm
      real(dp) :: idelta
      real(dp) :: tauxxip, tauxxim, dtauxxdx, tauxyjp, tauxyjm, dtauxydy
      real(dp) :: tauyxip, tauyxim, tauyyjp, tauyyjm, dtauyxdx, dtauyydy
      type(vector) :: lap_v

      if (constant_viscosity) then

         ! First compute the laplacian of v
         call lap_v%allocate(1)
         call laplacian(v, lap_v)

         ! Diffusion term is only mu*lap_v
         do k = base_grid%lo(3),base_grid%hi(3)
            do j = base_grid%lo(2),base_grid%hi(2)
               do i = base_grid%lo(1),base_grid%hi(1)
                  RHS%x%f(i,j,k) = RHS%x%f(i,j,k) + mu%f(i,j,k)*lap_v%x%f(i,j,k)/rhof%x%f(i,j,k)
                  RHS%y%f(i,j,k) = RHS%y%f(i,j,k) + mu%f(i,j,k)*lap_v%y%f(i,j,k)/rhof%y%f(i,j,k)
#if DIM==3
                  RHS%z%f(i,j,k) = RHS%z%f(i,j,k) + mu%f(i,j,k)*lap_v%z%f(i,j,k)/rhof%z%f(i,j,k)
#endif
               end do
            end do
         end do

         ! Free memory associated to laplacian of v
         call lap_v%destroy
      else

         ! For variable viscosity simulation viscosity can change due to multiphase or non-newtonain
         ! simulations.

#ifdef NN
         call gradient(v, D)
         call compute_viscosity(mu, D)
#endif

         idelta = 1.0_dp/base_grid%delta

         ! Add diffusion to RHS
         do k = base_grid%lo(3),base_grid%hi(3)
            do j = base_grid%lo(2),base_grid%hi(2)
               jm = j - 1
               jp = j + 1
               do i = base_grid%lo(1),base_grid%hi(1)
                  im = i - 1
                  ip = i + 1

                  tauxxip = 2.0_dp*mu%f(ip,j,k)*(v%x%f(ip,j,k) - v%x%f(i,j,k))*idelta
                  tauxxim = 2.0_dp*mu%f(i,j,k)*(v%x%f(i,j,k) - v%x%f(im,j,k))*idelta
                  dtauxxdx = (tauxxip - tauxxim)*idelta

                  tauxyjp = 0.25_dp*(mu%f(i,j,k) + mu%f(ip,j,k) + mu%f(i,jp,k) + mu%f(ip,jp,k))* &
                     ( (v%x%f(i,jp,k) - v%x%f(i,j,k))*idelta + (v%y%f(ip,j,k) - v%y%f(i,j,k))*idelta )
                  tauxyjm = 0.25_dp*(mu%f(i,jm,k) + mu%f(ip,jm,k) + mu%f(i,j,k) + mu%f(ip,j,k))* &
                     ( (v%x%f(i,j,k) - v%x%f(i,jm,k))*idelta + (v%y%f(ip,jm,k) - v%y%f(i,jm,k))*idelta )
                  dtauxydy = (tauxyjp - tauxyjm)*idelta

                  RHS%x%f(i,j,k) = RHS%x%f(i,j,k) + (dtauxxdx + dtauxydy)/rhof%x%f(i,j,k)

                  tauyxip = tauxyjp
                  tauyxim = 0.25_dp*(mu%f(im,j,k) + mu%f(i,j,k) + mu%f(im,jp,k) + mu%f(i,jp,k))* &
                     ( (v%x%f(im,jp,k) - v%x%f(im,j,k))*idelta + (v%y%f(i,j,k) - v%y%f(im,j,k))*idelta )
                  dtauyxdx = (tauyxip - tauyxim)*idelta

                  tauyyjp = 2.0_dp*mu%f(i,jp,k)*(v%y%f(i,jp,k) - v%y%f(i,j,k))*idelta
                  tauyyjm = 2.0_dp*mu%f(i,j,k)*(v%y%f(i,j,k) - v%y%f(i,jm,k))*idelta
                  dtauyydy = (tauyyjp - tauyyjm)*idelta

                  RHS%y%f(i,j,k) = RHS%y%f(i,j,k) + (dtauyxdx + dtauyydy)/rhof%y%f(i,j,k)
#if DIM==3
#endif
               end do
            end do
         end do
      end if

   end subroutine add_diffusion
   !========================================================================================
#ifdef MF
   !========================================================================================
   subroutine add_surface_tension(RHS)

      ! Add the surface tension force to the RHS of momentum using the CSF method.

      use class_Grid           , only : base_grid
      use volume_of_fluid , only : vof, curv
      use multiphase      , only : sigma

      ! In/Out variables
      type(vector), intent(inout) :: RHS

      ! Local variables
      integer  :: i, ip, j, jp, k, lo(3), hi(3)
#if DIM==3
      integer :: kp
#endif
      real(dp) :: idelta

      lo = base_grid%lo
      hi = base_grid%hi
      idelta = 1.0_dp/base_grid%delta

      do k = lo(3),hi(3)
#if DIM==3
         kp = k + 1
#endif
         do j = lo(2),hi(2)
            jp = j + 1
            do i = lo(1),hi(1)
               ip = i + 1

               RHS%x%f(i,j,k) = RHS%x%f(i,j,k) + sigma*0.5_dp*(curv%f(ip,j,k) + curv%f(i,j,k))*&
                  (vof%f(ip,j,k) - vof%f(i,j,k))*idelta/rhof%x%f(i,j,k)
               RHS%y%f(i,j,k) = RHS%y%f(i,j,k) + sigma*0.5_dp*(curv%f(i,jp,k) + curv%f(i,j,k))*&
                  (vof%f(i,jp,k) - vof%f(i,j,k))*idelta/rhof%y%f(i,j,k)
#if DIM==3
               RHS%z%f(i,j,k) = RHS%z%f(i,j,k) + sigma*0.5_dp*(curv%f(i,j,kp) + curv%f(i,j,k))*&
                  (vof%f(i,j,kp) - vof%f(i,j,k))*idelta/rhof%z%f(i,j,k)
#endif
            end do
         end do
      end do

   end subroutine add_surface_tension
   !========================================================================================
#endif
   !========================================================================================
   subroutine correct_velocity_field(dt)

      use class_Grid, only : base_grid
      use fields    , only : gradient
#ifdef MF
      use multiphase, only : irhomin
#endif

      ! In/Out variables
      real(dp), intent(in) :: dt

      ! Local variables
      integer :: i, j ,k

      ! Compute the gradient of the projector operator phi, using the fields dpdx dpdy to save
      ! memory
      call gradient(phi, grad_p)

      do k = base_grid%lo(3),base_grid%hi(3)
         do j = base_grid%lo(2),base_grid%hi(2)
            do i = base_grid%lo(1),base_grid%hi(1)
#ifdef MF
               v%x%f(i,j,k) = v%x%f(i,j,k) - grad_p%x%f(i,j,k)*dt*irhomin
               v%y%f(i,j,k) = v%y%f(i,j,k) - grad_p%y%f(i,j,k)*dt*irhomin
#if DIM == 3
               v%z%f(i,j,k) = v%z%f(i,j,k) - grad_p%z%f(i,j,k)*dt*irhomin
#endif
#else
               v%x%f(i,j,k) = v%x%f(i,j,k) - grad_p%x%f(i,j,k)*dt/rhof%x%f(i,j,k)
               v%y%f(i,j,k) = v%y%f(i,j,k) - grad_p%y%f(i,j,k)*dt/rhof%y%f(i,j,k)
#if DIM == 3
               v%z%f(i,j,k) = v%z%f(i,j,k) - grad_p%z%f(i,j,k)*dt/rhof%z%f(i,j,k)
#endif
#endif
            end do
         end do
      end do

      ! Update physical bc
      call v%apply_bc()

   end subroutine correct_velocity_field
   !========================================================================================

   !========================================================================================
   subroutine update_pressure()

     use mpi
     use class_Grid
#ifdef MF
      use multiphase, only : p_o

      ! If using the Dodd & Ferrante approximation for pressure, update pressure at 
      ! old timestep.
      p_o%f = p%f
#endif

      ! Update pressure
      p%f = p%f + phi%f
      
      ! Update boundary conditions
      call p%apply_bc()

   end subroutine update_pressure
   !========================================================================================

   !========================================================================================
   subroutine checks(dt, step)

      use mpi
      use class_Grid , only : base_grid
      use fields, only : divergence

      ! In/Out variables
      integer , intent(in) :: step
      real(dp), intent(in) :: dt

      ! Local variables
      integer      :: ierror, i, j, k
      real(dp)     :: vel, max_vel
      type(scalar) :: div

      call div%allocate()

      ! Compute the divergence of the velocity field
      call divergence(v, div)

      ! Compute the maximum value on each rank
      maxdiv = div%max_value()
      
      ! Destroy the divergence field
      call div%destroy()

      ! Compute maximum CFL
      max_vel = 0.0_dp
      do k = base_grid%lo(3),base_grid%hi(3)
         do j = base_grid%lo(2),base_grid%hi(2)
            do i = base_grid%lo(1),base_grid%hi(1)
#if DIM==3
               vel = abs(v%x%f(i,j,k)) + abs(v%y%f(i,j,k)) + abs(v%z%f(i,j,k))
#else
               vel = abs(v%x%f(i,j,k)) + abs(v%y%f(i,j,k))
#endif
               if (vel > max_vel) max_vel = vel
            end do
         end do
      end do

      ! Find maximum between all ranks
      call mpi_allreduce(mpi_in_place,max_vel,1,mpi_real8,mpi_max,mpi_comm_world,ierror)

      ! Actual CFL
      maxCFL = dt*max_vel/base_grid%delta

   end subroutine checks
   !========================================================================================

   !========================================================================================
   subroutine set_timestep(dt, U)
     
      ! Compute the timestep

      use class_Grid     , only : base_grid
#ifdef MF
      use constants , only : pi
      use multiphase, only : rho_0, rho_1, mu_0, mu_1, sigma
#endif
#ifdef NN
      use non_newtonian, only : mu_max
#endif

      ! In/Out variables
      real(dp), intent(in   ) :: U
      real(dp), intent(inout) :: dt

      ! Convective timestep
      dt_conv = CFL*base_grid%delta/U

      ! Viscous timestep
      dt_visc = 0.125_dp*base_grid%delta*base_grid%delta*density/viscosity
#ifdef NN
      dt_visc = 0.125_dp*base_grid%delta*base_grid%delta*density/mu_max
#endif

#if DIM==3
      dt_visc = 0.01_dp*base_grid%delta*base_grid%delta*base_grid%delta*density/viscosity
#endif

#ifdef MF
      dt_visc = 0.125_dp*base_grid%delta*base_grid%delta*min(rho_0/mu_0,rho_1/mu_1)
      dt = min(dt_conv, dt_visc)
      if (sigma > 0.0_dp) then
         dt_surf = sqrt(0.5_dp*(rho_0 + rho_1)*base_grid%delta**3/(pi*sigma + 1.0d-16))
         dt = min(dt, dt_surf)
      endif
#else
      ! Set the timestep as the minimum between the viscous
      dt = min(dt_conv, dt_visc)
#endif
      dt_o = dt

   end subroutine set_timestep
   !========================================================================================

   !========================================================================================
   subroutine update_timestep(dt)

      ! Compute the timestep

      use mpi
      use class_Grid     , only : base_grid
#ifdef MF
      use constants , only : pi
      use multiphase, only : rho_0, rho_1, mu_0, mu_1, sigma
#endif
#ifdef NN
      use non_newtonian, only : mu_max
#endif

      ! In/Out variables
      real(dp), intent(inout) :: dt

      ! Local variables
      integer :: i, j, k, ierror
      real(dp) :: vel, max_vel

      ! Save old timestep
      dt_o = dt

      ! Compute convective timestep
      max_vel = 0.0_dp
      do k = base_grid%lo(3),base_grid%hi(3)
         do j = base_grid%lo(2),base_grid%hi(2)
            do i = base_grid%lo(1),base_grid%hi(1)
#if DIM==3
               vel = abs(v%x%f(i,j,k)) + abs(v%y%f(i,j,k)) + abs(v%z%f(i,j,k))
#else
               vel = abs(v%x%f(i,j,k)) + abs(v%y%f(i,j,k))
#endif
               if (vel > max_vel) max_vel = vel
            end do
         end do
      end do

      ! Find maximum between all ranks
      call mpi_allreduce(mpi_in_place,max_vel,1,mpi_real8,mpi_max,mpi_comm_world,ierror)

      ! Compute convective timetep
      dt_conv = CFL*base_grid%delta/max_vel

      ! Set actual timestep
#ifdef MF
      dt = min(dt_conv,dt_visc,dt_surf)
#else
      dt = min(dt_conv,dt_visc)
#endif
      if (dt > 1.1_dp*dt_o) dt = 1.1_dp*dt_o

   end subroutine update_timestep
   !========================================================================================
   
   !========================================================================================
   subroutine print_navier_stokes_solver_status(log_id, step, time, dt)

      use class_Grid, only : base_grid

      ! In/Out variables
      integer , intent(in) :: log_id, step
      real(dp), intent(in) :: time, dt

      ! Print solver informations
      if (base_grid%rank == 0) write(log_id,10) 'step: ', step, 'time: ', time, 'dt: ', dt, &
         'maxdiv: ', maxdiv, 'maxCFL: ', maxCFL

10    format(A6,I7,1x,A6,E13.6,1x,A4,E13.6,1x,A8,E13.6,1x,A9,1x,E13.6)

   end subroutine print_navier_stokes_solver_status
   !========================================================================================

   !========================================================================================
   subroutine allocate_navier_stokes_fields()

      use class_Grid, only : base_grid
      use io

      ! Allocate all fields for solving Navier-Stokes equations.
      call p%allocate(1)
      call phi%allocate(1)
      call rho%allocate(1)
      call rhof%allocate(0)
      call mu%allocate(1)
      call v%allocate(1)
      call dv%allocate()
      call dv_o%allocate()
      call grad_p%allocate()
      call S%allocate()
#ifdef NN
      call D%allocate()
#endif

      ! Set material properties
      rho%f = density
      mu%f  = viscosity

      ! Based on the boundary conditions on the physical domain select boundary conditions
      ! for all fields
      ! Left boundary
      if (base_grid%boundary_conditions(1)%s == 'Periodic') then
         p%bc%type_l = 0
         phi%bc%type_l = 0
         rho%bc%type_l = 0
         mu%bc%type_l = 0
         v%x%bc%type_l = 0
         v%y%bc%type_l = 0
#if DIM==3
         v%z%bc%type_l = 0
#endif
      elseif (base_grid%boundary_conditions(1)%s == 'Wall') then
         p%bc%type_l = 2
         phi%bc%type_l = 2
         rho%bc%type_l = 2
         mu%bc%type_l = 2
         v%x%bc%type_l = 1
         v%y%bc%type_l = 1
#if DIM==3
         v%z%bc%type_l = 1
#endif
      elseif (base_grid%boundary_conditions(1)%s == 'Inflow') then
         p%bc%type_l = 2
         phi%bc%type_l = 2
         rho%bc%type_l = 2
         mu%bc%type_l = 2
         v%x%bc%type_l = 1
         v%y%bc%type_l = 2
#if DIM==3
         v%z%bc%type_l = 2
#endif
      elseif (base_grid%boundary_conditions(1)%s == 'Outflow') then
         p%bc%type_l = 1
         phi%bc%type_l = 1
         rho%bc%type_l = 2
         mu%bc%type_l = 2
         v%x%bc%type_l = 2
         v%y%bc%type_l = 2
#if DIM==3
         v%z%bc%type_l = 2
#endif
      else
         call print_error_message('ERROR: wrong bc on left boundary')
      endif

      ! Right boundary
      if (base_grid%boundary_conditions(2)%s == 'Periodic') then
         p%bc%type_r = 0
         phi%bc%type_r = 0
         rho%bc%type_r = 0
         mu%bc%type_r = 0
         v%x%bc%type_r = 0
         v%y%bc%type_r = 0
#if DIM==3
         v%z%bc%type_r = 0
#endif
      elseif (base_grid%boundary_conditions(2)%s == 'Wall') then
         p%bc%type_r = 2
         phi%bc%type_r = 2
         rho%bc%type_r = 2
         mu%bc%type_r = 2
         v%x%bc%type_r = 1
         v%y%bc%type_r = 1
#if DIM==3
         v%z%bc%type_r = 1
#endif
      elseif (base_grid%boundary_conditions(2)%s == 'Inflow') then
         p%bc%type_r = 2
         phi%bc%type_r = 2
         rho%bc%type_r = 2
         mu%bc%type_r = 2
         v%x%bc%type_r = 1
         v%y%bc%type_r = 2
#if DIM==3
         v%z%bc%type_r = 2
#endif
      elseif (base_grid%boundary_conditions(2)%s == 'Outflow') then
         p%bc%type_r = 1
         phi%bc%type_r = 1
         rho%bc%type_r = 2
         mu%bc%type_r = 2
         v%x%bc%type_r = 2
         v%y%bc%type_r = 2
#if DIM==3
         v%z%bc%type_r = 2
#endif
      else
         call print_error_message('ERROR: wrong bc on right boundary')
      endif

      ! Bottom boundary
      if (base_grid%boundary_conditions(3)%s == 'Periodic') then
         p%bc%type_b = 0
         phi%bc%type_b = 0
         rho%bc%type_b = 0
         mu%bc%type_b = 0
         v%x%bc%type_b = 0
         v%y%bc%type_b = 0
#if DIM==3
         v%z%bc%type_b = 0
#endif
      elseif (base_grid%boundary_conditions(3)%s == 'Wall') then
         p%bc%type_b = 2
         phi%bc%type_b = 2
         rho%bc%type_b = 2
         mu%bc%type_b = 2
         v%x%bc%type_b = 1
         v%y%bc%type_b = 1
#if DIM==3
         v%z%bc%type_b = 1
#endif
      elseif (base_grid%boundary_conditions(3)%s == 'Inflow') then
         p%bc%type_b = 2
         phi%bc%type_b = 2
         rho%bc%type_b = 2
         mu%bc%type_b = 2
         v%x%bc%type_b = 1
         v%y%bc%type_b = 1
#if DIM==3
         v%z%bc%type_b = 1
#endif
      elseif (base_grid%boundary_conditions(3)%s == 'Outflow') then
         p%bc%type_b = 1
         phi%bc%type_b = 1
         rho%bc%type_b = 2
         mu%bc%type_b = 2
         v%x%bc%type_b = 2
         v%y%bc%type_b = 2
#if DIM==3
         v%z%bc%type_b = 2
#endif
      else
         call print_error_message('ERROR: wrong bc on bottom boundary')
      endif

      ! top boundary
      if (base_grid%boundary_conditions(4)%s == 'Periodic') then
         p%bc%type_t = 0
         phi%bc%type_t = 0
         rho%bc%type_t = 0
         mu%bc%type_t = 0
         v%x%bc%type_t = 0
         v%y%bc%type_t = 0
#if DIM==3
         v%z%bc%type_t = 0
#endif
      elseif (base_grid%boundary_conditions(4)%s == 'Wall') then
         p%bc%type_t = 2
         phi%bc%type_t = 2
         rho%bc%type_t = 2
         mu%bc%type_t = 2
         v%x%bc%type_t = 1
         v%y%bc%type_t = 1
#if DIM==3
         v%z%bc%type_t = 1
#endif
      elseif (base_grid%boundary_conditions(4)%s == 'Inflow') then
         p%bc%type_t = 2
         phi%bc%type_t = 2
         rho%bc%type_t = 2
         mu%bc%type_t = 2
         v%x%bc%type_t = 1
         v%y%bc%type_t = 1
#if DIM==3
         v%z%bc%type_t = 1
#endif
      elseif (base_grid%boundary_conditions(4)%s == 'Outflow') then
         p%bc%type_t = 1
         phi%bc%type_t = 1
         rho%bc%type_t = 2
         mu%bc%type_t = 2
         v%x%bc%type_t = 2
         v%y%bc%type_t = 2
#if DIM==3
         v%z%bc%type_t = 2
#endif
      else
         call print_error_message('ERROR: wrong bc on top boundary')
      endif

#if DIM==3
      ! FOR NOW 3D Simultaions supports only periodic BC in z
      ! Front boundary
      if (base_grid%boundary_conditions(5)%s == 'Periodic') then
         p%bc%type_f = 0
         phi%bc%type_f = 0
         rho%bc%type_f = 0
         mu%bc%type_f = 0
         v%x%bc%type_f = 0
         v%y%bc%type_f = 0
         v%z%bc%type_f = 0
      else
         call print_error_message('ERROR: wrong bc on top boundary')
      endif

      ! End boundary
      if (base_grid%boundary_conditions(6)%s == 'Periodic') then
         p%bc%type_e = 0
         phi%bc%type_e = 0
         rho%bc%type_e = 0
         mu%bc%type_e = 0
         v%x%bc%type_e = 0
         v%y%bc%type_e = 0
         v%z%bc%type_e = 0
      else
         call print_error_message('ERROR: wrong bc on top boundary')
      endif
#endif

   end subroutine allocate_navier_stokes_fields
   !========================================================================================

   !========================================================================================
   subroutine destroy_navier_stokes_solver()

      ! Free the memory allocated by the Navier-Stokes solver.
      call p%destroy()
      call phi%destroy()
      call rho%destroy()
      call mu%destroy()
      call v%destroy()
      call dv%destroy()
      call dv_o%destroy()
      call grad_p%destroy()
      call S%destroy()
      call rhof%destroy

   end subroutine destroy_navier_stokes_solver
   !========================================================================================

end module navier_stokes
