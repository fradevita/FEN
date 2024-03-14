!> This module contains all procedures for the Eulerian Immersed Boundary Method.
module eulerian_ibm_mod

    use precision_mod     , only : dp
    use grid_mod          , only : grid
    use eulerian_solid_mod, only : eulerian_solid_pointer, eulerian_solid

    implicit none

    ! The Eulerian fields are:
    ! closest: an integer field which gives the closest solid body
    ! ibm_index: an integer field which gives the tagging of the point.
    ! Tag legend
    ! -1 = extrapolation point
    !  0 = solid point
    !  1 = interface interpolation point
    !  2 = fluid point
    ! Indexes are:
    ! 1 = i
    ! 2 = j
    ! 3 = k
    ! 4 = location (0 = cell center, 1 = x face, 2 = y face, 3 = z face)
    integer, dimension(:,:,:,:), allocatable :: closest
    integer, dimension(:,:,:,:), allocatable :: ibm_index

    ! Define the interpolation procedures
    abstract interface
        function velocity_interpolation(i, j, k, v, solid, dir) result(vs)
            use precision_mod     , only : dp
            use scalar_mod          
            use eulerian_solid_mod, only : eulerian_solid
            integer              , intent(in) :: i, j, k, dir
            type(scalar)         , intent(in) :: v
            class(eulerian_solid), intent(in) :: solid
            real(dp)                          :: vs
        end function velocity_interpolation
    end interface

    ! Define procedure pointer
    procedure(velocity_interpolation), pointer :: interpolate_velocity => Null()

contains

    !===============================================================================================
    subroutine init_eulerian_ibm(solid_list, G)

        ! This subroutine initialize all eulerian fields.

        ! In/Out variables
        type(Eulerian_Solid_pointer), intent(in) :: solid_list(:)
        type(grid)                  , intent(in) :: G

        ! Allocate memoroy for the fields
#if DIM==3
        allocate(  closest(G%lo(1)-1:G%hi(1)+1,G%lo(2)-1:G%hi(2)+1,G%lo(3)-1:G%hi(3)+1,0:3))
        allocate(ibm_index(G%lo(1)-1:G%hi(1)+1,G%lo(2)-1:G%hi(2)+1,G%lo(3)-1:G%hi(3)+1,0:3))
#else
        allocate(  closest(G%lo(1)-1:G%hi(1)+1,G%lo(2)-1:G%hi(2)+1,G%lo(3)-1:G%hi(3)+1,0:2))
        allocate(ibm_index(G%lo(1)-1:G%hi(1)+1,G%lo(2)-1:G%hi(2)+1,G%lo(3)-1:G%hi(3)+1,0:2))
#endif

        ! Initialize the eulerian fields
        call tag_cells(solid_list)

#if DIM==3
        interpolate_velocity => velocity_interpolation_3D
#else
        interpolate_velocity => velocity_interpolation_2D
#endif

    end subroutine init_eulerian_ibm
    !===============================================================================================

    !===============================================================================================
    subroutine tag_cells(solid_list)

        ! This subroutine compute the eulerian fields phi, ibm_index, norm and cbl
        ! for the identification of solid bodies inside the Eulerian grid

        use global_mod, only : istagger, stagger

        ! In/Out variables
        type(eulerian_solid_pointer), intent(in) :: solid_list(:)

        ! Local variables
        integer :: nb, i, j, k, dir, im, ip, jm, jp, b, clb_i, ndir
        real(dp) :: x, y, z, delta
        real(dp), dimension(:), allocatable :: distance
        real(dp), dimension(:,:,:,:), allocatable :: phi
        type(grid) :: bg
#if DIM==3
        integer  :: kp, km
#endif

        bg = solid_list(1)%pS%G
        delta = bg%delta

        nb = size(solid_list)
        allocate(distance(nb))

        ! Set the number of directions
#if DIM==3
        ndir = 3
#else
        ndir = 2
#endif

        allocate(phi(bg%lo(1)-1:bg%hi(1)+1,bg%lo(2)-1:bg%hi(2)+1,bg%lo(3)-1:bg%hi(3)+1,0:ndir))
        dir_cycle: do dir = 0,ndir

            do k = bg%lo(3)-1,bg%hi(3)+1
                do j = bg%lo(2)-1,bg%hi(2)+1
                    do i = bg%lo(1)-1,bg%hi(1)+1

                        ! Local coordinates
                        ! x = (i - istagger(1,dir))*delta
                        ! y = (j - istagger(2,dir))*delta
                        ! z = (k - istagger(3,dir))*delta

                        x = bg%x(i) + stagger(1, dir)*delta
                        y = bg%y(j) + stagger(2, dir)*delta
                        z = bg%z(k) + stagger(3, dir)*delta

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

            do k = bg%lo(3),bg%hi(3)
#if DIM==3
                kp = k + 1
                km = k - 1
#endif
                do j = bg%lo(2),bg%hi(2)
                    jp = j + 1
                    jm = j - 1
                    do i = bg%lo(1),bg%hi(1)
                        ip = i + 1
                        im = i - 1

                        ! Tag cell
                        if (phi(i,j,k,dir) <= 0.0_dp) then

                            !if (phi(ip,j,k,dir) > 0.0_dp .or. phi(im,j,k,dir) > 0.0_dp .or. &
                            !   phi(i,jp,k,dir) > 0.0_dp .or. phi(i,jm,k,dir) > 0.0_dp) then
                                ! extrapolation point
                            !    ibm_index(i,j,k,dir) = -1
                            !else
                                ! Solid point
                                ibm_index(i,j,k,dir) = 0
                            !endif
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
            call update_halo_bc_ibm_index(ibm_index(:,:,:,dir), bg)

        end do dir_cycle

    end subroutine tag_cells
    !===============================================================================================

    !===============================================================================================
    subroutine forcing_velocity(v, solid_list, Fe, dt)

        ! Force the velocity field v based on the location of the solid body.
        ! Use this function instead of compute_ibm_forcing
    
        use vector_mod
        use scalar_mod
        use global_mod, only : Ndim, istagger
    
        ! In/Out variables
        real(dp)                    , intent(in   ) :: dt
        type(vector), target        , intent(inout) :: v
        type(vector), target        , intent(inout) :: Fe
        type(Eulerian_Solid_pointer), intent(in   ) :: solid_list(:)

        ! Local variables
        integer                :: i, j, k, b, dir
        real(dp)               :: delta, x, y, z, vs
        type(scalar) , pointer :: vi, Fi
        
        delta = v%G%delta
        call Fe%setToValue(0.0_dp)

        do k = v%G%lo(3),v%G%hi(3)
            do j = v%G%lo(2),v%G%hi(2)
                do i = v%G%lo(1),v%G%hi(1)

                    velocity_component_cycle: do dir = 1,Ndim

                        ! Select scalar velocity component vi
                        if (dir == 1) then
                            vi => v%x
                            Fi => Fe%x
                        elseif (dir == 2) then
                            vi => v%y
                            Fi => Fe%y
                        elseif (dir ==3) then
                            vi => v%z
                            Fi => Fe%z
                        endif

                        if (ibm_index(i,j,k,dir) == 2) then
                            ! Fluid point, do nothing
                            Fi%f(i,j,k) = 0.0_dp
                        elseif (ibm_index(i,j,k,dir) == 1) then
                            ! Interface point
                            b = closest(i,j,k,dir)
                            vs = interpolate_velocity(i, j, k, vi, solid_list(b)%pS, dir)
                            ! Evaluate the integral of the forcing 
                            Fi%f(i,j,k) = (vs - vi%f(i,j,k))/dt
                            ! Force i-th velocity component
                            vi%f(i,j,k) = vs
                        else 
                            ! Solid point
                            x = (i - istagger(1,dir))*delta + v%G%origin(1)
                            y = (j - istagger(2,dir))*delta + v%G%origin(2)
                            z = (k - istagger(3,dir))*delta + v%G%origin(3)
                            b = closest(i,j,k,dir)
                            vs = solid_list(b)%pS%velocity([x, y, z], dir)
                            ! Evaluate the integral of the forcing 
                            Fi%f(i,j,k) = (vs - vi%f(i,j,k))/dt
                            ! Force i-th velocity component 
                            vi%f(i,j,k) = vs
                        endif

                    end do velocity_component_cycle
    
              end do
           end do
        end do
    
        ! Apply boundary conditions after forcing
        call v%update_ghost_nodes
        
    end subroutine forcing_velocity
    !===============================================================================================
    
    !===============================================================================================
    function velocity_interpolation_2D(i, j, k, v, solid, dir) result(fl)

        ! This function compute the interpolated velocity component (f) on the forcing point
        ! with coordinates (i,j,k)

        use global_mod, only : istagger
        use scalar_mod
        use utils_mod , only : linear_interpolation, bilinear_interpolation

        ! In/Out variables
        integer              , intent(in) :: i, j, k, dir
        type(scalar)         , intent(in) :: v
        class(eulerian_solid), intent(in) :: solid

        ! Local variables
        integer  :: i2, j2
        real(dp) :: delta, x, y, s, a, b, q, fl, xl, yl, xx, yy, velb, xb, yb, nx, ny, x2, y2, nn(3)

        delta = v%G%delta

        ! Physical coordinates
        x = (i - istagger(1,dir))*delta + v%G%origin(1)
        y = (j - istagger(2,dir))*delta + v%G%origin(2)

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
        x2 = (i2 - istagger(1,dir))*delta + v%G%origin(1)
        y2 = (j2 - istagger(2,dir))*delta + v%G%origin(2)

        ! If one norm component is zero interpolate along cartesian directions
        if (abs(nx) <= 1.0e-12_dp) then
            ! Auxiliary point in (x,y2)
            q = v%f(i,j2,k)
            ! Inerpolate in the forcing point
            fl = linear_interpolation(s, 0.0_dp, delta, velb, q)
        elseif (abs(ny) <= 1.0e-12_dp) then
            ! Auxiliary point in (x2,y)
            q = v%f(i2,j,k)
            ! Interpolate in the forcing point
            fl = linear_interpolation(s, 0.0_dp, delta, velb, q)
        else

            ! Norm line equation
            a = ny/nx
            b = y - a*x

            if (ibm_index(i2,j,k,dir) < 2) then
                ! Intersection between normal line and line at y = y2: (xx, y2)
                xx = (y2 - b)/a

                ! Valocity in the auxiliary point
                q = linear_interpolation(xx, min(x,x2), max(x,x2), v%f(min(i,i2),j2,k), v%f(max(i,i2),j2,k))

                ! Interpolation in forcing point
                fl = linear_interpolation(s, 0.0_dp, solid%distance([xx, y2, 0.0_dp]), velb, q)

            elseif (ibm_index(i,j2,k,dir) < 2) then
                ! Intersection between normal line and line at x = x2: (x2, yy)
                yy = a*x2 + b

                ! Velocity in the auxiliary point
                q = linear_interpolation(yy, min(y,y2), max(y,y2), v%f(i2,min(j,j2),k), v%f(i2,max(j,j2),k))

                ! Interpolation in the forcing point
                fl = linear_interpolation(s, 0.0_dp, solid%distance([x2, yy, 0.0_dp]), velb, q)

            else
                ! Bilinear interpolation

                ! Auxiliary point at distance s in normal direction
                xl = x + nx*s
                yl = y + ny*s

                ! Normalization
                xl = (xl - min(x,x2))/delta
                yl = (yl - min(y,y2))/delta

                q = bilinear_interpolation(xl, yl, v%f(min(i,i2),min(j,j2),k), &
                                                   v%f(max(i,i2),min(j,j2),k), &
                                                   v%f(min(i,i2),max(j,j2),k), &
                                                   v%f(max(i,i2),max(j,j2),k))

                ! Interpolation in the forcing point
                fl = linear_interpolation(s, 0.0_dp, 2.0_dp*s, velb, q)
            end if

        endif

    end function velocity_interpolation_2D
    !===============================================================================================

#if DIM==3
    !===============================================================================================
    function velocity_interpolation_3D(i, j, k, v, solid, dir) result(fl)

        ! This function compute the interpolated velocity component (f) on the forcing point
        ! with coordinates (i,j,k)

        use global_mod, only : istagger
        use utils_mod , only : bilinear_interpolation, trilinear_interpolation
        use scalar_mod 

        ! In/Out variables
        integer              , intent(in) :: i, j, k !< Cell index
        type(scalar)         , intent(in) :: v
        class(eulerian_solid), intent(in) :: solid   !< Solid object 
        integer              , intent(in) :: dir     !< velocity direction

        ! Varibili locali
        integer :: i2, j2, k2
        real(dp) :: x, y, z, nn(3), nx, ny, nz, s, xb, yb, zb, velb, x2, y2, z2, q, fl
        real(dP) :: xv, yv, zv, xl, yl, zl, d, delta
        type(grid) :: bg

        bg = v%G
        delta = bg%delta

        ! Local coordinates
        x = (i - istagger(1,dir))*delta + bg%origin(1)
        y = (j - istagger(2,dir))*delta + bg%origin(2)
        z = (k - istagger(3,dir))*delta + bg%origin(3)

        ! Local normal
        nn = solid%norm([x, y, z])
        nx = nn(1)
        ny = nn(2)
        nz = nn(3)

        ! Local distance
        s = solid%distance([x, y, z])

        ! Coordinates on the solid surface
        xb = x - nx*s
        yb = y - ny*s
        zb = z - nz*s

        ! Solid body velocity on the surface
        velb = solid%velocity([xb, yb, zb], dir)

        ! Select neighbours based on local norm
        i2 = i + int(sign(1.0_dp,nx))
        j2 = j + int(sign(1.0_dp,ny))
        k2 = k + int(sign(1.0_dp,nz))
        x2 = (i2 - istagger(1,dir))*delta + bg%origin(1)
        y2 = (j2 - istagger(2,dir))*delta + bg%origin(2)
        z2 = (k2 - istagger(3,dir))*delta + bg%origin(3)

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
                q = bilinear_interpolation(xl, zl, v%f(min(i,i2),j2,min(k,k2)), &
                                                   v%f(max(i,i2),j2,min(k,k2)), &
                                                   v%f(min(i,i2),j2,max(k,k2)), &
                                                   v%f(max(i,i2),j2,max(k,k2)))

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
                q = bilinear_interpolation(xl, yl, v%f(min(i,i2),min(j,j2),k2), &
                                                   v%f(max(i,i2),min(j,j2),k2), &
                                                   v%f(min(i,i2),max(j,j2),k2), &
                                                   v%f(max(i,i2),max(j,j2),k2))

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
                q = bilinear_interpolation(yl, zl, v%f(i2,min(j,j2),min(k,k2)), &
                                                   v%f(i2,max(j,j2),min(k,k2)), &
                                                   v%f(i2,min(j,j2),max(k,k2)), &
                                                   v%f(i2,max(j,j2),max(k,k2)))

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
                q = bilinear_interpolation(xl, yl, v%f(min(i,i2),min(j,j2),k2), &
                                                   v%f(max(i,i2),min(j,j2),k2), &
                                                   v%f(min(i,i2),max(j,j2),k2), &
                                                   v%f(max(i,i2),max(j,j2),k2))

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
                q = bilinear_interpolation(yl, zl, v%f(i2,min(j,j2),min(k,k2)), &
                                                   v%f(i2,max(j,j2),min(k,k2)), &
                                                   v%f(i2,min(j,j2),max(k,k2)), &
                                                   v%f(i2,max(j,j2),max(k,k2)))

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
                q = bilinear_interpolation(xl, zl, v%f(min(i,i2),j2,min(k,k2)), &
                                                   v%f(max(i,i2),j2,min(k,k2)), &
                                                   v%f(min(i,i2),j2,max(k,k2)), &
                                                   v%f(max(i,i2),j2,max(k,k2)))

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
            q = trilinear_interpolation(xl,yl,zl,v%f(min(i,i2),min(j,j2),min(k,k2)), &
                                                 v%f(max(i,i2),min(j,j2),min(k,k2)), &
                                                 v%f(min(i,i2),max(j,j2),min(k,k2)), &
                                                 v%f(min(i,i2),min(j,j2),max(k,k2)), &
                                                 v%f(max(i,i2),max(j,j2),min(k,k2)), &
                                                 v%f(max(i,i2),min(j,j2),max(k,k2)), &
                                                 v%f(min(i,i2),max(j,j2),max(k,k2)), &
                                                 v%f(max(i,i2),max(j,j2),max(k,k2)))

            ! Interpolated velocity in the forcing point
            fl = velb + (q - velb)*0.5_dp
        end if

    end function velocity_interpolation_3D
    !===============================================================================================
#endif

    !===============================================================================================
    subroutine compute_hydrodynamic_loads(solid, v, p, mu, rho, g, Fe, density)

        use scalar_mod
        use vector_mod
        use euclidean_mod
        use fields_mod

        ! In/Out variables
        real(dp)             , intent(in   ) :: g(:), density
        class(eulerian_solid), intent(inout) :: solid
        type(vector)         , intent(in   ) :: v, Fe
        type(scalar)         , intent(in   ) :: p, mu, rho

        ! Local variables
        integer      :: i, j, k
        real(dp)     :: T(3), r(3)
        type(scalar) :: Mx, My, Mz 

        if (solid%use_probes) then
            call loads_from_probes(solid, v, p, mu, rho, g)
        else
            ! TODO: for now single phase only, use density
            solid%center_of_mass%Fh(1) = -Fe%x%integral()*density + &
                                            solid%volume()*solid%center_of_mass%A(1)*density
            solid%center_of_mass%Fh(2) = -Fe%y%integral()*density + &
                                            solid%volume()*solid%center_of_mass%A(2)*density
#if DIM==3
            solid%center_of_mass%Fh(3) = -Fe%z%integral()*density + &
                                            solid%volume()*solid%center_of_mass%A(3)*density
#else
            solid%center_of_mass%Fh(3) = 0.0_dp
#endif
            ! For the torque must evaluate integral of local IBM force torque
            call Mx%allocate(solid%G)
            call My%allocate(solid%G)
            call Mz%allocate(solid%G)
            do k = solid%G%lo(3),solid%G%hi(3)
                do j = solid%G%lo(2),solid%G%hi(2)
                    do i = solid%G%lo(1),solid%G%hi(1)
#if DIM==3
                        r = [ solid%G%x(i) - solid%center_of_mass%X(1), &
                              solid%G%y(j) - solid%center_of_mass%X(2), &
                              solid%G%z(k) - solid%center_of_mass%X(3)]
                        T = crossProduct(r, [Fe%x%f(i,j,k), Fe%y%f(i,j,k), Fe%z%f(i,j,k)])
#else
                        r = [ solid%G%x(i) - solid%center_of_mass%X(1), &
                              solid%G%y(j) - solid%center_of_mass%X(2), &
                              0.0_dp]
                        T = crossProduct(r, [Fe%x%f(i,j,k), Fe%y%f(i,j,k), 0.0_dp])
#endif
                        Mx%f(i,j,k) = T(1)
                        My%f(i,j,k) = T(2)
                        Mz%f(i,j,k) = T(3)
                    end do
                end do
            end do
            solid%center_of_mass%Fh(4) = density*(-Mx%integral() + &
                                            solid%IM(4)*solid%center_of_mass%A(4)/solid%density)
            solid%center_of_mass%Fh(5) = density*(-My%integral() + &
                                            solid%IM(5)*solid%center_of_mass%A(5)/solid%density)
            solid%center_of_mass%Fh(6) = density*(-Mz%integral() + &
                                            solid%IM(6)*solid%center_of_mass%A(6)/solid%density)
            call Mx%destroy()
            call My%destroy()
            call Mz%destroy()
        endif

    end subroutine
    !===============================================================================================

    !===============================================================================================
    subroutine loads_from_probes(solid, v, p, mu, rho, g)

        ! This subroutine compute hydrodynamic forces for each solid body using the probes 
        ! method.

        use mpi
        use utils_mod           , only : bilinear_interpolation, integrate_1D
        use global_mod          , only : pi, Ndim
        use scalar_mod
        use vector_mod
        use tensor_mod

        ! In/Out variables
        real(dp)             , intent(in   ) :: g(2)
        class(eulerian_solid), intent(inout) :: solid
        type(vector)         , intent(in   ) :: v
        type(scalar)         , intent(in   ) :: p, mu, rho

        ! Local variables
        integer :: i, j, k, ip, im, jp, jm, n, error, nr, build_matrix
        integer , dimension(3) :: Ic, I1, I2
        real(dp) :: delta, idelta, ds, pl, Dxxl, Dxyl, Dyyl
        real(dp) :: Fvx, Fvy, Fpx, Fpy, accl_tot(2), rr(3)
        real(dp), dimension(:), allocatable :: Fvxl, Fvyl, Fpxl, Fpyl, Fxl, Fyl, s, rx, ry
        real(dp), dimension(3) :: loc_norm, X_probe, X_center, X, Xl
        type(scalar) :: pB
        type(tensor) :: D

        delta = v%G%delta
        idelta = 1.0_dp/delta

        ! Compute shear stress and pressure
        call pB%allocate(v%G, 1)
        call D%allocate(v%G, 1)

        do k = v%G%lo(3),v%G%hi(3)
            do j = v%G%lo(2),v%G%hi(2)
                jp = j + 1
                jm = j - 1
                do i = v%G%lo(1),v%G%hi(1)
                    ip = i + 1
                    im = i - 1
                    pB%f(i,j,k) = p%f(i,j,k)
                    D%x%x%f(i,j,k) = (v%x%f(i,j,k) - v%x%f(im,j,k))*idelta
                    D%x%y%f(i,j,k) = 0.5_dp*(0.25_dp*(v%x%f(i,jp,k) + v%x%f(im,jp,k) -     &
                        v%x%f(i,jm,k) - v%x%f(im,jm,k))*idelta + 0.25_dp*(v%y%f(ip,j,k) + &
                        v%y%f(ip,jm,k) - v%y%f(im,j,k) - v%y%f(im,jm,k))*idelta)
                    D%y%y%f(i,j,k) = (v%y%f(i,j,k) - v%y%f(i,jm,k))*idelta
                end do
            end do
        end do

        ! Update BC and ghost nodes
        call pB%update_ghost_nodes()
        call D%update_ghost_nodes()

        ! Set the number of auxiliary points
        nr = size(solid%surface_points(1,:))

        ! Allocate array for local stress
        allocate(Fvxl(nr))
        allocate(Fvyl(nr))
        allocate(Fpxl(nr))
        allocate(Fpyl(nr))
        allocate(Fxl(nr))
        allocate(Fyl(nr))

        ! Perform a cycle over the surface points nr
        cycle_nr: do n = 1,nr

            ! Local norm
            loc_norm = solid%norm(solid%surface_points(:,n))

            ! Initial position of the probe: 1.25 times grid spacing from the solid surface
            ds = 1.25_dp*delta
            X_probe = solid%surface_points(:,n) + ds*loc_norm

            ! Check if the probe is out of the domain and in case traslate it.
            X_probe = traslate(X_probe, v%G)

            ! Search neighbours for the interplation
            build_matrix = 1
            buid_matrix: do while (build_matrix == 1)

                ! Closest cell center to the probe
                Ic = v%G%closest_grid_node(X_probe, 0)

                ! Phyiscal coordinates in this cell
                X_center = (Ic - 0.5_dp)*delta

                ! Select neighbours based on the local norm
                if (X_probe(1) > X_center(1) .and. X_probe(2) > X_center(2)) then
                    I1 = Ic
                    I2 = Ic + [1, 1, 0]
                    X = X_center
                elseif (X_probe(1) > X_center(1) .and. X_probe(2) < X_center(2)) then
                    I1 = Ic + [0, -1, 0]
                    I2 = Ic + [1, 0, 0]
                    X = X_center + [0.0_dp, -delta, 0.0_dp]
                elseif (X_probe(1) < X_center(1) .and. X_probe(2) < X_center(2)) then
                    I1 = Ic + [-1, -1, 0]
                    I2 = Ic
                    X = X_center + [-delta, -delta, 0.0_dp]
                elseif (X_probe(1) < X_center(1) .and. X_probe(2) > X_center(2)) then
                    I1 = Ic + [-1, 0, 0]
                    I2 = Ic + [0, 1, 0]
                    X = X_center + [-delta, 0.0_dp, 0.0_dp]
                else
                    I1 = [-10, -10, 0]
                    I2 = [-10, -10, 0]
                    X = [-10.0_dp, -10.0_dp, 0.0_dp]
                    print *, 'Error while searching neighbours'
                endif

                ! Check that all neighbours are fluid points
                ! The check must be done on the proper rank
                build_matrix = 0
                if (Ic(2) >= v%G%lo(2) .and. Ic(2) <= v%G%hi(2)) then
                    if (ibm_index(I2(1),I1(2),1,0) < 2 .or. &
                        ibm_index(I1(1),I2(2),1,0) < 2 .or. ibm_index(I2(1),I2(2),1,0) < 2) then
                        ! Move the probe
                        ds = ds*1.25_dp
                        X_probe = solid%surface_points(:,n) + ds*loc_norm
                        ! Check if the probes is outside of the domain
                        X_probe = traslate(X_probe, v%G)
                        build_matrix = 1
                    else

                    endif
                else
                    ds = 0.0_dp
                    X_probe = 0.0_dp
                endif
                call mpi_allreduce(mpi_in_place,X_probe,3,mpi_real8,mpi_sum,mpi_comm_world,error)
                call mpi_allreduce(mpi_in_place,build_matrix,1,mpi_int,mpi_sum,mpi_comm_world,error)
                call mpi_allreduce(mpi_in_place,ds,1,mpi_real8,mpi_sum,mpi_comm_world,error)

            end do buid_matrix

            ! Interpolate shear tensor and pressure on the probe.
            ! Select the rank containing the probe.
            if (Ic(2) >= v%G%lo(2) .and. Ic(2) <= v%G%hi(2)) then
                Xl   = (X_probe - X)*idelta

                pl   = bilinear_interpolation(Xl(1), Xl(2),    pB%f(I1(1),I1(2),1),    pB%f(I2(1),I1(2),1), &
                                                               pB%f(I1(1),I2(2),1),    pB%f(I2(1),I2(2),1))
                Dxxl = bilinear_interpolation(Xl(1), Xl(2), D%x%x%f(I1(1),I1(2),1), D%x%x%f(I2(1),I1(2),1), &
                                                            D%x%x%f(I1(1),I2(2),1), D%x%x%f(I2(1),I2(2),1))
                Dxyl = bilinear_interpolation(Xl(1), Xl(2), D%x%y%f(I1(1),I1(2),1), D%x%y%f(I2(1),I1(2),1), &
                                                            D%x%y%f(I1(1),I2(2),1), D%x%y%f(I2(1),I2(2),1))
                Dyyl = bilinear_interpolation(Xl(1), Xl(2), D%y%y%f(I1(1),I1(2),1), D%y%y%f(I2(1),I1(2),1), &
                                                            D%y%y%f(I1(1),I2(2),1), D%y%y%f(I2(1),I2(2),1))

                ! The pressure on the solid body boundary depends on the acceleartion
                accl_tot(1) = solid%acceleration(solid%surface_points(:,n), 1)
                accl_tot(2) = solid%acceleration(solid%surface_points(:,n), 2)
                pl = pl + rho%f(Ic(1),Ic(2),1)*ds*((accl_tot(1) - g(1))*loc_norm(1) + &
                                                   (accl_tot(2) - g(2))*loc_norm(2))

                ! Compute forces components: Fv is the viscous force, Fp the pressure force
                Fvxl(n) = 2.0_dp*mu%f(Ic(1),Ic(2),1)*(Dxxl*loc_norm(1) + Dxyl*loc_norm(2))
                Fvyl(n) = 2.0_dp*mu%f(Ic(1),Ic(2),1)*(Dxyl*loc_norm(1) + Dyyl*loc_norm(2))
                Fpxl(n) = -pl*loc_norm(1)
                Fpyl(n) = -pl*loc_norm(2)

            else

                ! Set forces to zero on other ranks
                Fvxl(n) = 0.0_dp
                Fvyl(n) = 0.0_dp
                Fpxl(n) = 0.0_dp
                Fpyl(n) = 0.0_dp

            endif

        end do cycle_nr

        ! Broadcast forces
        call mpi_allreduce(mpi_in_place,Fvxl,nr,mpi_real8,mpi_sum,mpi_comm_world,error)
        call mpi_allreduce(mpi_in_place,Fvyl,nr,mpi_real8,mpi_sum,mpi_comm_world,error)
        call mpi_allreduce(mpi_in_place,Fpxl,nr,mpi_real8,mpi_sum,mpi_comm_world,error)
        call mpi_allreduce(mpi_in_place,Fpyl,nr,mpi_real8,mpi_sum,mpi_comm_world,error)

        ! Add viscous and pressure force
        Fxl = Fvxl + Fpxl
        Fyl = Fvyl + Fpyl

        ! *** Forces integral ***
        ! Init variables to zero
        solid%center_of_mass%Fh = 0.0_dp
        Fvx = 0.0_dp
        Fvy = 0.0_dp
        Fpx = 0.0_dp
        Fpy = 0.0_dp

        ! s is a curvilinear coordinate for integration
        allocate(s(nr))
        s = 0.0_dp
        do n = 2,nr
            s(n) = s(n-1) + sqrt((solid%surface_points(1,n) - solid%surface_points(1,n-1))**2 + &
                                 (solid%surface_points(2,n) - solid%surface_points(2,n-1))**2)
        end do

        ! [rx,ry] is the position vector of the force with respect to the rotation center
        ! used to compute the torque.
        allocate(rx(nr))
        allocate(ry(nr))
        rx = 0.0_dp
        ry = 0.0_dp
        do n = 1,nr
            rr = [solid%surface_points(1,n) - solid%rotation_center(1),           &
                  solid%surface_points(1,n) - solid%rotation_center(1) + v%G%Lx,  &
                  solid%surface_points(1,n) - solid%rotation_center(1) - v%G%Lx]
            rx(n) = rr(minloc(abs(rr),1))
            ry(n) = solid%surface_points(2,n) - solid%rotation_center(2)
        end do

        ! Integrate forces
        solid%center_of_mass%Fh(1) = integrate_1D(s, Fxl, .true.)
        solid%center_of_mass%Fh(2) = integrate_1D(s, Fyl, .true.)
        solid%center_of_mass%Fh(3:5) = 0.0_dp ! TODO: For now only 2D, extend to 3D
        solid%center_of_mass%Fh(6) = integrate_1D(s, rx*Fyl - ry*Fxl, .true.)

        ! Evalute also froces contribution (viscous and pressure)
        solid%center_of_mass%Fv(1) = integrate_1D(s, Fvxl, .true.)
        solid%center_of_mass%Fp(1) = integrate_1D(s, Fpxl, .true.)
        solid%center_of_mass%Fv(2) = integrate_1D(s, Fvyl, .true.)
        solid%center_of_mass%Fp(2) = integrate_1D(s, Fpyl, .true.)

        ! Free the memory
        deallocate(Fvxl, Fvyl, Fpxl, Fpyl, Fxl, Fyl, s, rx, ry)
        call pB%destroy()
        call D%destroy()

    end subroutine 
    !===============================================================================================

    !===============================================================================================
    function traslate(X, G) result(X1)
   
        ! In/out variables
        real(dp)  , intent(in) :: X(3)
        type(grid), intent(in) :: G
        real(dp)               :: X1(3)
    
        X1 = X
    
        if (G%periodic_bc(1) .eqv. .true.) then
           if (X(1) > G%Lx) then
              X1(1) = X(1) - G%Lx
           elseif (X(1) < 0.0_dp) then
              X1(1) = X(1) + G%Lx
           endif
        endif
    
    end function traslate
    !===============================================================================================

    !===============================================================================================
    subroutine update_halo_bc_ibm_index(ff, base_grid)

        use decomp_2d

        ! In/Out variables
        type(grid), intent(in) :: base_grid
        integer, intent(inout) :: &
            ff(base_grid%lo(1)-1:base_grid%hi(1)+1,base_grid%lo(2)-1:base_grid%hi(2)+1,base_grid%lo(3)-1:base_grid%hi(3)+1)

        real(dp) :: &
            f(base_grid%lo(1)-1:base_grid%hi(1)+1,base_grid%lo(2)-1:base_grid%hi(2)+1,base_grid%lo(3)-1:base_grid%hi(3)+1)

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
            if (nrank == 0) then
                f(:,base_grid%lo(2)-1,:) = f(:,base_grid%lo(2),:)
            endif
            if (base_grid%rank == base_grid%nranks - 1) then
                f(:,base_grid%hi(2)+1,:) = f(:,base_grid%hi(2),:)
            endif
        endif

        ff = int(f)

    end subroutine update_halo_bc_ibm_index
    !===============================================================================================

    !===============================================================================================
    subroutine destroy_ibm

            ! Free the allocated memory
            deallocate(closest, ibm_index)

    end subroutine destroy_ibm
    !===============================================================================================

end module
