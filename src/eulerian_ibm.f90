!> This module contains all procedures for the Eulerian Immersed Boundary Method.
module eulerian_ibm

    use precision, only : dp

    implicit none

    ! The Eulerian fields are:
    ! closest: an integer field which gives the closest solid body
    ! ibm_index: an integer field which gives the tagging of the point.
    ! Tag legend
    !  0 = solid point
    !  1 = interface interpolation point
    !  2 = fluid point
    ! Indexes are:
    ! 1 = i
    ! 2 = j
    ! 3 = k
    ! 4 = location (1 = cell center, 2 = x face, 3 = y face, 4 = z face)
    integer, dimension(:,:,:,:), allocatable :: closest
    integer, dimension(:,:,:,:), allocatable :: ibm_index

    ! Staggering of the variables
    real(dp), dimension(3,4), parameter :: stagger = &
        reshape([0.5_dp, 0.5_dp, 0.5_dp, &
                0.0_dp, 0.5_dp, 0.5_dp, &
                0.5_dp, 0.0_dp, 0.5_dp, &
                0.5_dp, 0.5_dp, 0.0_dp], shape(stagger))

    ! Define the interpolation procedures
    abstract interface
        function velocity_interpolation(i, j, k, f, solid, dir) result(vs)
            use precision           , only : dp
            use class_Grid          , only : bg => base_grid
            use class_eulerian_solid, only : eulerian_solid
            integer              , intent(in) :: i, j, k, dir
            real(dp)             , intent(in) :: f(bg%lo(1)-1:bg%hi(1)+1, &
                                                bg%lo(2)-1:bg%hi(2)+1, &
                                                bg%lo(3)-1:bg%hi(3)+1)
            class(Eulerian_Solid), intent(in) :: solid
            real(dp)                               :: vs
        end function velocity_interpolation
    end interface

    ! Define procedure pointer
    procedure(velocity_interpolation), pointer :: interpolate_velocity => Null()

contains

    !========================================================================================
    subroutine init_eulerian_ibm(solid_list)

        ! This subroutine initialize all eulerian fields.

        use class_Grid          , only : bg => base_grid
        use class_Eulerian_Solid, only : Eulerian_Solid_pointer

        ! In/Out variables
        type(Eulerian_Solid_pointer), intent(in) :: solid_list(:)

        ! Allocate memoroy for the fields
#if DIM==3
        allocate(closest(bg%lo(1):bg%hi(1),bg%lo(2):bg%hi(2),bg%lo(3):bg%hi(3),4))
        allocate(ibm_index(bg%lo(1)-1:bg%hi(1)+1,bg%lo(2)-1:bg%hi(2)+1,bg%lo(3)-1:bg%hi(3)+1,4))
#else
        allocate(closest(bg%lo(1):bg%hi(1),bg%lo(2):bg%hi(2),bg%lo(3):bg%hi(3),3))
        allocate(ibm_index(bg%lo(1)-1:bg%hi(1)+1,bg%lo(2)-1:bg%hi(2)+1,bg%lo(3)-1:bg%hi(3)+1,3))
#endif

        ! Initialize the eulerian fields
        call tag_cells(solid_list)

#if DIM==3
        interpolate_velocity => velocity_interpolation_3D
#else
        interpolate_velocity => velocity_interpolation_2D
#endif

    end subroutine init_eulerian_ibm
    !========================================================================================

    !========================================================================================
    subroutine tag_cells(solid_list)

        ! This subroutine compute the eulerian fields phi, ibm_index, norm and cbl
        ! for the identification of solid bodies inside the Eulerian grid

        use class_Grid          , only : bg => base_grid
        use class_Eulerian_Solid, only : Eulerian_Solid_pointer

        ! In/Out variables
        type(Eulerian_Solid_pointer), intent(in) :: solid_list(:)

        ! Local variables
        integer  :: nb, i, j, k, dir, im, ip, jm, jp, b, clb_i, ndir
        real(dp) :: x, y, z, delta
        real(dp), dimension(:), allocatable :: distance
        real(dp), dimension(:,:,:,:), allocatable :: phi
#if DIM==3
        integer  :: kp, km
#endif

        delta = bg%delta

        nb = size(solid_list)
        allocate(distance(nb))

        ! Set the number of directions
#if DIM==3
        ndir = 4
#else
        ndir = 3
#endif

        allocate(phi(bg%lo(1)-1:bg%hi(1)+1,bg%lo(2)-1:bg%hi(2)+1,bg%lo(3)-1:bg%hi(3)+1,ndir))
        dir_cycle: do dir = 1,ndir

            do k = bg%lo(3),bg%hi(3)
                do j = bg%lo(2),bg%hi(2)
                    do i = bg%lo(1),bg%hi(1)

                        ! Local coordinates
                        x = (i - stagger(1,dir))*delta
                        y = (j - stagger(2,dir))*delta
                        z = (k - stagger(3,dir))*delta

                        ! Compute distance from every solid body
                        do b = 1,nb
                            distance(b) = Solid_list(b)%pS%distance([x,y,z])
                        end do
                        phi(i,j,k,dir) = minval(distance)

                        ! Select the closest solid body clb to the local point
                        clb_i = minloc(distance,1)
                        closest(i,j,k,dir) = clb_i
                    end do
                end do
            end do

            ! Update boundary conditions and halo
            call update_halo_bc_solid(phi(:,:,:,dir))

            do k = bg%lo(3),bg%hi(3)
                do j = bg%lo(2),bg%hi(2)
                    jp = j + 1
                    jm = j - 1
                    do i = bg%lo(1),bg%hi(1)
                        ip = i + 1
                        im = i - 1

                        ! Tag cell
                        if (phi(i,j,k,dir) <= 0.0_dp) then
                            ! Solid point
                            ibm_index(i,j,k,dir) = 0
                        else
#if DIM==3
                        if (phi(ip,j,k,dir) < 0.0_dp .or. phi(im,j,k,dir) < 0.0_dp .or. &
                            phi(i,jp,k,dir) < 0.0_dp .or. phi(i,jm,k,dir) < 0.0_dp .or. &
                            phi(i,j,kp,dir) < 0.0_dp .or. phi(i,j,km,dir) < 0.0_dp) then
                            ! interface point
                            ibm_index(i,j,k,dir) = 1
                        else
                            ! fluid point
                            ibm_index(i,j,k,dir) = 2
                        endif
#else
                        if (phi(ip,j,k,dir) < 0.0_dp .or. phi(im,j,k,dir) < 0.0_dp .or. &
                            phi(i,jp,k,dir) < 0.0_dp .or. phi(i,jm,k,dir) < 0.0_dp) then
                            ! interface point
                            ibm_index(i,j,k,dir) = 1
                        else
                            ! fluid point
                            ibm_index(i,j,k,dir) = 2
                        endif
#endif
                    endif
                  end do
                end do
            end do

            ! BC and ghost
            !call update_halo_bc_ibm_index(ibm_index(:,:,:,dir))

        end do dir_cycle

   end subroutine tag_cells
   !========================================================================================

    !========================================================================================
    subroutine compute_ibm_forcing(v, RHS, solid_list, dt, F)

        !> Compute the force field due to all eulerian solid.

        use mpi
        use class_Grid          , only : base_grid
        use class_Vector        , only : vector
        use class_Eulerian_Solid, only : Eulerian_Solid_pointer

        ! In/Out variables
        type(vector)                , intent(in   ) :: v
        type(vector)                , intent(in   ) :: RHS
        type(Eulerian_Solid_pointer), intent(in   ) :: solid_list(:)
        real(dp)                    , intent(in   ) :: dt
        type(vector)                , intent(inout) :: F

        ! Local variables
        integer  :: i, j, k, b
        real(dp) :: delta, x, y, z, vs

        delta = base_grid%delta

        do k = base_grid%lo(3),base_grid%hi(3)
            do j = base_grid%lo(2),base_grid%hi(2)
                do i = base_grid%lo(1),base_grid%hi(1)

                ! X component of the forcing
                if (ibm_index(i,j,k,2) == 2) then
                    ! Fluid point do nothing
                    F%x%f(i,j,k) = 0.0_dp
                elseif (ibm_index(i,j,k,2) == 1) then
                    ! Interface point
                    b = closest(i,j,k,2)
                    vs = interpolate_velocity(i, j, k, v%x%f, solid_list(b)%pS, 2)
                    F%x%f(i,j,k) = (vs - v%x%f(i,j,k))/dt - RHS%x%f(i,j,k)
                else
                    ! Solid point
                    x = (i - 0.0_dp)*delta
                    y = (j - 0.5_dp)*delta
                    z = (k - 0.5_dp)*delta
                    b = closest(i,j,k,2)
                    vs = solid_list(b)%pS%velocity([x, y, z], 1)
                    F%x%f(i,j,k) = (vs - v%x%f(i,j,k))/dt - RHS%x%f(i,j,k)
                endif

                ! Y component of the forcing
                if (ibm_index(i,j,k,3) == 2) then
                    ! Fluid point, do nothing
                    F%y%f(i,j,k) = 0.0_dp
                elseif (ibm_index(i,j,k,3) == 1) then
                    ! Interface point
                    b = closest(i,j,k,3)
                    vs = interpolate_velocity(i, j, k, v%y%f, solid_list(b)%pS, 3)
                    F%y%f(i,j,k) = (vs - v%y%f(i,j,k))/dt - RHS%y%f(i,j,k)
                else
                    ! Solid point
                    x = (i - 0.5_dp)*delta
                    y = (j - 0.0_dp)*delta
                    z = (k - 0.5_dp)*delta
                    b = closest(i,j,k,3)
                    vs = solid_list(b)%pS%velocity([x, y, z], 2)
                    F%y%f(i,j,k) = (vs - v%y%f(i,j,k))/dt - RHS%y%f(i,j,k)
                endif

!#if DIM==3
!                ! Force z component of velocity
!                if (ibm_index(i,j,k,4) == 2) then
!                   ! Fluid point, do nothing
!                elseif (ibm_index(i,j,k,4) == 1) then
!                   ! Interface point
!                   v%z%f(i,j,k) = interpolate_velocity(i, j, k, v%z%f, 4)
!                else
!                   ! Solid point
!                   x = (i - 0.5_dp)*delta
!                   y = (j - 0.5_dp)*delta
!                   z = (k - 0.0_dp)*delta
!                   v%z%f(i,j,k) = rigid_body_velocity(x, y, z, 4, closest(i,j,k,4))
!                endif
!#endif
                end do
            end do
        end do

    end subroutine compute_ibm_forcing
    !======================================================================================

   !======================================================================================
   function velocity_interpolation_2D(i, j, k, f, solid, dir) result(fl)

      ! This function compute the interpolated velocity component (f) on the forcing point
      ! with coordinates (i,j,k)

      use class_Grid          , only : base_grid
      use class_Eulerian_Solid, only : eulerian_solid
      use utils               , only : bilinear

      ! In/Out variables
      integer , intent(in) :: i, j, k, dir
      real(dp), intent(in) :: f(base_grid%lo(1)-1:base_grid%hi(1)+1, &
                                base_grid%lo(2)-1:base_grid%hi(2)+1, &
                                base_grid%lo(3)-1:base_grid%hi(3)+1)
      class(eulerian_solid), intent(in) :: solid

      ! Local variables
      integer  :: i2, j2
      real(dp) :: delta, x, y, s, a, b, q, fl, xl, yl, xx, yy, velb, xb, yb, nx, ny, x2, y2, dist, nn(2)

      delta = base_grid%delta

      ! Physical coordinates
      x = (i - stagger(1,dir))*delta
      y = (j - stagger(2,dir))*delta

      ! Local normal vector
      nn = solid%norm([x, y, 0.0_dp])
      nx = nn(1)
      ny = nn(2)

      ! Local distance
      s = solid%distance([x, y, 0.0_dp])

      ! Physical coordiantes on solid boundary
      xb = x - nx*s
      yb = y - ny*s

      ! Velocity on solid boundary
      velb = solid%velocity([xb, yb, 0.0_dp], dir)

      ! Select nodes for interpolation based on local norm
      i2 = i + int(sign(1.0_dp,nx))
      j2 = j + int(sign(1.0_dp,ny))
      x2 = (i2 - stagger(1,dir))*delta
      y2 = (j2 - stagger(2,dir))*delta

      ! If one norm component is zero interpolate along cartesian directions
      if (nx == 0.0_dp) then
         ! Interpolate in y
         q = f(i,j2,k)
         fl = velb + (q - velb)*s/(s + delta)
      elseif (ny == 0) then
         ! Interpolate in x
         q = f(i2,j,k)
         fl = velb + (q - velb)*s/(s + delta)
      else

         ! Norm line equation
         a = ny/nx
         b = y - a*x

         if (ibm_index(i2,j,k,dir) < 2) then
            ! Virtual point with interpolation in y

            ! Intersezione tra retta normale e retta a z = z2 ha coordinate yy e z2
            xx = (y2 - b)/a

            ! Valore della velocità nel punto virtuale
            q = f(min(i,i2),j2,k) + (f(max(i,i2),j2,k) - f(min(i,i2),j2,k))*(xx - min(x,x2))/delta

            ! Interpolazione nel punto di forzaggio
            dist = sqrt((xx - xb)**2 + (y2 - yb)**2)
            fl = velb + (q - velb)*s/dist

         elseif (ibm_index(i,j2,k,dir) < 2) then
            ! Punto virutale con interpolazione in y

            ! Intersezione tra retta normale e retta a y = y2 ha coordinate y2 e zz
            yy = a*x2 + b

            ! Valore della velocità nel punto virtuale
            q = f(i2,min(j,j2),k) + (f(i2,max(j,j2),k) - f(i2,min(j,j2),k))*(yy - min(y,y2))/delta

            ! Interpolazione nel punto di forzaggio
            dist = sqrt( (x2 - xb)**2 + (yy - yb)**2 )
            fl = velb + (q - velb)*s/dist

         else
            ! Bilinear interpolation

            ! Virtual point in
            xl = x + nx*s
            yl = y + ny*s

            ! Normalization
            xl = (xl - min(x,x2))/delta
            yl = (yl - min(y,y2))/delta

            q = bilinear(xl, yl, f(min(i,i2),min(j,j2),k), f(max(i,i2),min(j,j2),k), &
               f(min(i,i2),max(j,j2),k), f(max(i,i2),max(j,j2),k))

            fl = velb + (q - velb)*0.5_dp
         end if

      endif

   end function velocity_interpolation_2D
   !========================================================================================

#if DIM==3
   !========================================================================================
   function velocity_interpolation_3D(i,j,k,f,dir) result(fl)

      ! This function compute the interpolated velocity component (f) on the forcing point
      ! with coordinates (i,j,k)

      use class_Grid, only : bg => base_grid
      use utils     , only : bilinear, trilinear

      ! In/Out variables
      integer , intent(in) :: i, j, k, dir
      real(dp), intent(in) :: f(bg%lo(1)-1:bg%hi(1)+1,bg%lo(2)-1:bg%hi(2)+1,bg%lo(3)-1:bg%hi(3)+1)

      ! Varibili locali
      integer :: i2, j2, k2
      real(dp) :: x, y, z, nx, ny, nz, s, xb, yb, zb, velb, x2, y2, z2, q, fl
      real(dP) :: xv, yv, zv, xl, yl, zl, d, delta

      delta = bg%delta

      ! Local coordinates
      x = (i - stagger(1,dir))*delta
      y = (j - stagger(2,dir))*delta
      z = (k - stagger(3,dir))*delta

      ! Local normal
      nx = norm(i,j,k,dir,1)
      ny = norm(i,j,k,dir,2)
      nz = norm(i,j,k,dir,3)

      ! Local distance
      s = phi(i,j,k,dir)

      ! Coordinates on the solid surface
      xb = x - nx*s
      yb = y - ny*s
      zb = z - nz*s

      ! Solid body velocity on the surface
      velb = rigid_body_velocity(xb, yb, zb, dir, closest(i,j,k,dir))

      ! Select neighbours based on local norm
      i2 = i + int(sign(1.0_dp,nx))
      j2 = j + int(sign(1.0_dp,ny))
      k2 = k + int(sign(1.0_dp,nz))
      x2 = (i2 - stagger(1,dir))*delta
      y2 = (j2 - stagger(2,dir))*delta
      z2 = (k2 - stagger(3,dir))*delta

      ! Check if some of the neighbours is a forcing point
      if (ibm_index(i2,j,k,dir) < 2 ) then
         if (abs(ny) >= abs(nz)) then
            ! Virtual point in the xz plane at y = y2
            xv = x + nx*delta/abs(ny)
            yv = y + ny*delta/abs(ny)
            zv = z + nz*delta/abs(ny)
            xl = (xv - min(x,x2))/delta
            zl = (zv - min(z,z2))/delta

            ! Bilinear interpolation in the virtual point
            q = bilinear(xl, zl, f(min(i,i2),j2,min(k,k2)), f(max(i,i2),j2,min(k,k2)), &
               f(min(i,i2),j2,max(k,k2)), f(max(i,i2),j2,max(k,k2)))

            ! Distance between solid surface and virtual point
            d = sqrt((xv - xb)**2 + (yv - yb)**2 + (zv - zb)**2)

            ! Interpolated velocity in the forcing point
            fl = velb + (q - velb)*s/d
         else
            ! Virtual point in the xy plane at z = z2
            xv = x + nx*delta/abs(nz)
            yv = y + ny*delta/abs(nz)
            zv = z + nz*delta/abs(nz)
            xl = (xv - min(x,x2))/delta
            yl = (yv - min(y,y2))/delta

            ! Bilinear interpolation in the virtual point
            q = bilinear(xl, yl, f(min(i,i2),min(j,j2),k2), f(max(i,i2),min(j,j2),k2), &
               f(min(i,i2),max(j,j2),k2), f(max(i,i2),max(j,j2),k2))

            ! Distance between solid surface and virtual point
            d = sqrt((xv - xb)**2 + (yv - yb)**2 + (zv - zb)**2)

            ! Interpolated velocity in the forcing point
            fl = velb + (q - velb)*s/d
         endif
      elseif (ibm_index(i,j2,k,dir) < 2) then
         if (abs(nx) >= abs(nz)) then
            ! Virtual point in the yz plane at x = x2
            xv = x + nx*delta/abs(nx)
            yv = y + ny*delta/abs(nx)
            zv = z + nz*delta/abs(nx)
            yl = (yv - min(y,y2))/delta
            zl = (zv - min(z,z2))/delta

            ! Bilinear interpolation in the virtual point
            q = bilinear(yl, zl, f(i2,min(j,j2),min(k,k2)), f(i2,max(j,j2),min(k,k2)), &
               f(i2,min(j,j2),max(k,k2)), f(i2,max(j,j2),max(k,k2)))

            ! Distance between solid surface and virtual point
            d = sqrt((xv - xb)**2 + (yv - yb)**2 + (zv - zb)**2)

            ! Interpolated velocity in the forcing point
            fl = velb + (q - velb)*s/d
         else
            ! Virtual point in the xy plane at z = z2
            xv = x + nx*delta/abs(nz)
            yv = y + ny*delta/abs(nz)
            zv = z + nz*delta/abs(nz)
            xl = (xv - min(x,x2))/delta
            yl = (yv - min(y,y2))/delta

            ! Bilinear interpolation in the virtual point
            q = bilinear(xl, yl, f(min(i,i2),min(j,j2),k2), f(max(i,i2),min(j,j2),k2), &
               f(min(i,i2),max(j,j2),k2), f(max(i,i2),max(j,j2),k2))

            ! Distance between solid surface and virtual point
            d = sqrt((xv - xb)**2 + (yv - yb)**2 + (zv - zb)**2)

            ! Interpolated velocity in the forcing point
            fl = velb + (q - velb)*s/d
         endif
      elseif (ibm_index(i,j,k2,dir) < 2) then
         if (abs(nx) >= abs(ny)) then
            ! Virtual point in the yz plane at x = x2
            xv = x + nx*delta/abs(nx)
            yv = y + ny*delta/abs(nx)
            zv = z + nz*delta/abs(nx)
            yl = (yv - min(y,y2))/delta
            zl = (zv - min(z,z2))/delta
            ! Bilinear interpolation in the virtual point
            q = bilinear(yl, zl, f(i2,min(j,j2),min(k,k2)), f(i2,max(j,j2),min(k,k2)), &
               f(i2,min(j,j2),max(k,k2)), f(i2,max(j,j2),max(k,k2)))

            ! Distance between solid surface and virtual point
            d = sqrt((xv - xb)**2 + (yv - yb)**2 + (zv - zb)**2)

            ! Interpolated velocity in the forcing point
            fl = velb + (q - velb)*s/d
         else
            ! Virtual point in the xz at y = y2
            xv = x + nx*delta/abs(ny)
            yv = y + ny*delta/abs(ny)
            zv = z + nz*delta/abs(ny)
            xl = (xv - min(x,x2))/delta
            zl = (zv - min(z,z2))/delta

            ! Bilinear interpolation in the virtual point
            q = bilinear(xl, zl, f(min(i,i2),j2,min(k,k2)), f(max(i,i2),j2,min(k,k2)), &
               f(min(i,i2),j2,max(k,k2)), f(max(i,i2),j2,max(k,k2)))

            ! Distance between solid surface and virtual point
            d = sqrt((xv - xb)**2 + (yv - yb)**2 + (zv - zb)**2)

            ! Interpolated velocity in the forcing point
            fl = velb + (q - velb)*s/d
         endif
      else
         ! If all neighbours are fluid points perform trilinear interpolation
         ! in the virual point:
         xv = x + nx*s
         yv = y + ny*s
         zv = z + nz*s

         ! Normalization
         xl = (xv - min(x,x2))/delta
         yl = (yv - min(y,y2))/delta
         zl = (zv - min(z,z2))/delta

         ! Trilinear interpolation
         q = trilinear(xl,yl,zl,f(min(i,i2),min(j,j2),min(k,k2)), &
            f(max(i,i2),min(j,j2),min(k,k2)), &
            f(min(i,i2),max(j,j2),min(k,k2)), &
            f(min(i,i2),min(j,j2),max(k,k2)), &
            f(max(i,i2),max(j,j2),min(k,k2)), &
            f(max(i,i2),min(j,j2),max(k,k2)), &
            f(min(i,i2),max(j,j2),max(k,k2)), &
            f(max(i,i2),max(j,j2),max(k,k2)))

         ! Interpolated velocity in the forcing point
         fl = velb + (q - velb)*0.5_dp
      end if

   end function velocity_interpolation_3D
   !========================================================================================
#endif

   !========================================================================================
   subroutine update_halo_bc_solid(f)

      use decomp_2d
      use class_Grid, only : base_grid

      ! In/Out variables
      real(mytype), intent(inout) :: &
         f(base_grid%lo(1)-1:base_grid%hi(1)+1,base_grid%lo(2)-1:base_grid%hi(2)+1,base_grid%lo(3)-1:base_grid%hi(3)+1)

      ! Local variables
      real(mytype), dimension(:,:,:), allocatable :: fh

      ! Call decomp_2d function to update halos
      call update_halo(f(base_grid%lo(1):base_grid%hi(1),base_grid%lo(2):base_grid%hi(2),base_grid%lo(3):base_grid%hi(3)), &
         fh, level = 1, opt_global = .true.)

      ! Copy into f
      f(base_grid%lo(1):base_grid%hi(1),base_grid%lo(2)-1:base_grid%hi(2)+1,base_grid%lo(3)-1:base_grid%hi(3)+1) = &
         fh(base_grid%lo(1):base_grid%hi(1),base_grid%lo(2)-1:base_grid%hi(2)+1,base_grid%lo(3)-1:base_grid%hi(3)+1)

      ! Free memroy
      deallocate(fh)

      ! X direction
      if (base_grid%periodic_bc(1)) then
         f(base_grid%lo(1)-1,:,:) = f(base_grid%hi(1),:,:)
         f(base_grid%hi(1)+1,:,:) = f(base_grid%lo(1),:,:)
      else
         f(base_grid%lo(1)-1,:,:) = f(base_grid%lo(1),:,:)
         f(base_grid%hi(1)+1,:,:) = f(base_grid%hi(1),:,:)
      endif

      ! If using periodic bc in y and 1 proc need to overwrite the physical bc
      if (base_grid%nranks == 1 .and. base_grid%periodic_bc(2) .eqv. .true.) then
         f(:,base_grid%lo(2)-1,:) = f(:,base_grid%hi(2),:)
         f(:,base_grid%hi(2)+1,:) = f(:,base_grid%lo(2),:)
      endif

      ! If non periodic in y select physical bc
      if (base_grid%periodic_bc(2) .eqv. .false.) then
         if (base_grid%rank == 0) then
            f(:,base_grid%lo(2)-1,:) = f(:,base_grid%lo(2),:)
         endif
         if (base_grid%rank == base_grid%nranks-1) then
            f(:,base_grid%hi(2)+1,:) = f(:,base_grid%hi(2),:)
         endif
      endif

   end subroutine update_halo_bc_solid
   !========================================================================================

   !========================================================================================
   subroutine update_halo_bc_ibm_index(ff)

      use decomp_2d
      use class_Grid, only : base_grid

      ! In/Out variables
      integer, intent(inout) :: &
         ff(base_grid%lo(1)-2:base_grid%hi(1)+2,base_grid%lo(2)-2:base_grid%hi(2)+2,base_grid%lo(3)-2:base_grid%hi(3)+2)

      real(dp) :: &
         f(base_grid%lo(1)-2:base_grid%hi(1)+2,base_grid%lo(2)-2:base_grid%hi(2)+2,base_grid%lo(3)-2:base_grid%hi(3)+2)

      ! Local variables
      real(dp), dimension(:,:,:), allocatable :: fh

      f = ff*1.0_dp

      ! Call decomp_2d function to update halos
      call update_halo(f(base_grid%lo(1):base_grid%hi(1),base_grid%lo(2):base_grid%hi(2),base_grid%lo(3):base_grid%hi(3)), &
         fh, level = 1, opt_global = .true.)

      ! Copy into f
      f(base_grid%lo(1):base_grid%hi(1),base_grid%lo(2)-1:base_grid%hi(2)+1,base_grid%lo(3)-1:base_grid%hi(3)+1) = &
         fh(base_grid%lo(1):base_grid%hi(1),base_grid%lo(2)-1:base_grid%hi(2)+1,base_grid%lo(3)-1:base_grid%hi(3)+1)

      ! Free memroy
      deallocate(fh)

      ! X direction
      if (base_grid%periodic_bc(1)) then
         f(base_grid%lo(1)-1,:,:) = f(base_grid%hi(1),:,:)
         !f(base_grid%lo(1)-2,:,:) = f(base_grid%hi(1)-1,:,:)
         f(base_grid%hi(1)+1,:,:) = f(base_grid%lo(1),:,:)
         !f(base_grid%hi(1)+2,:,:) = f(base_grid%lo(1)+1,:,:)
      else
         f(base_grid%lo(1)-1,:,:) = f(base_grid%lo(1),:,:)
         !f(base_grid%lo(1)-1,:,:) = f(base_grid%lo(1)+1,:,:)
         f(base_grid%hi(1)+1,:,:) = f(base_grid%hi(1),:,:)
         !f(base_grid%hi(1)+1,:,:) = f(base_grid%hi(1)-1,:,:)
      endif

      ! If using periodic bc in y and 1 proc need to overwrite the physical bc
      if (base_grid%nranks == 1 .and. base_grid%periodic_bc(2) .eqv. .true.) then
         f(:,base_grid%lo(2)-1,:) = f(:,base_grid%hi(2),:)
         !f(:,base_grid%lo(2)-2,:) = f(:,base_grid%hi(2)-1,:)
         f(:,base_grid%hi(2)+1,:) = f(:,base_grid%lo(2),:)
         !f(:,base_grid%hi(2)+2,:) = f(:,base_grid%lo(2)+1,:)
      endif

      ! If non periodic in y select physical bc
      if (base_grid%periodic_bc(2) .eqv. .false.) then
         if (nrank == 0) then
            f(:,base_grid%lo(2)-1,:) = f(:,base_grid%lo(2),:)
            !f(:,base_grid%lo(2)-2,:) = f(:,base_grid%lo(2)+1,:)
         endif
         if (base_grid%rank == base_grid%nranks - 1) then
            f(:,base_grid%hi(2)+1,:) = f(:,base_grid%hi(2),:)
            !f(:,base_grid%hi(2)+2,:) = f(:,base_grid%hi(2)-1,:)
         endif
      endif

      ff = int(f)

   end subroutine update_halo_bc_ibm_index
   !========================================================================================

   !========================================================================================
   subroutine destroy_ibm

        ! Free the allocated memory
        deallocate(closest, ibm_index)

   end subroutine destroy_ibm
   !========================================================================================

end module eulerian_ibm
