program main

    ! Program to test the velocity interpolation routine of the Eulerian Immersed Boundary method.

    use mpi
    use constants            , only : pi
    use precision            , only : dp
    use class_Grid           , only : base_grid, bc_type
    use class_Vector
    use class_eulerian_solid
    use class_test_solid
    use eulerian_ibm

    implicit none

    ! Variables
    integer       :: ierror, out_id, res, Nx, Ny, Nz, l, i, j, k
    real(dp)      :: Lx, Ly, Lz, origin(3), e_u, e_v, e_u_max(100), e_v_max(100), vs
#if DIM==3
    real(dp)      :: e_w, e_w_max(100)
#endif
    type(vector)  :: v
#if DIM==3
    type(bc_type) :: bc(6)
#else
    type(bc_type) :: bc(4)
#endif

    ! Create one solid of type circle, it must be target
    type(test_solid), target :: C
    type(eulerian_solid_pointer) :: solid_list(1)

    ! Initialize MPI
    call mpi_init(ierror)

    ! The domain is a unit squared box
    Lx = 1.0_dp
    Ly = 1.0_dp
    origin = [0.0_dp, 0.0_dp, 0.0_dp]
    bc(1)%s = 'Periodic'
    bc(2)%s = 'Periodic'
    bc(3)%s = 'Periodic'
    bc(4)%s = 'Periodic'
#if DIM==3
    bc(5)%s = 'Periodic'
    bc(6)%s = 'Periodic'
#endif

    C = test_solid(R = 0.25_dp)
    solid_list(1)%pS => C

    ! Perform computation for increasing resolution
    open(newunit = out_id, file = 'out.txt')

    resolution_loop: do res = 3,7

        ! Set the resolution
        Nx = 2**(res)
        Ny = 2**(res)
#if DIM==3
        Nz = Nx
        Lx = 1.0_dp
#else
        Nz = 1
        Lz = Lx*float(Nz)/float(Nx)
#endif
        ! Create the grid
        call base_grid%setup(Nx, Ny, Nz, Lx, Ly, Lz, origin, 1, 1, bc)

        ! Set the velocity field
        call v%allocate(1)
        call init_velocity(v)

        ! Perform interpolation for several random initial position of the solid body
        e_u_max = 0.0_dp
        e_v_max = 0.0_dp
#if DIM==3
        e_w_max = 0.0_dp
#endif

        call init_eulerian_ibm(solid_list)

        circle_position: do l = 1,100

            ! Set the center of the circle
            C%X(1) = 0.5_dp + rand()*0.1_dp
            C%X(2) = 0.5_dp + rand()*0.1_dp
#if DIM==3
            C%X(3) = 0.5_dp + rand()*0.1_dp
#endif

            ! Tag eulerian cells
            call tag_cells(solid_list)
            
            ! Compute maximum error
            do k = base_grid%lo(3),base_grid%hi(3)
                do j = base_grid%lo(2),base_grid%hi(2)
                    do i = base_grid%lo(1),base_grid%hi(1)

                        if (ibm_index(i,j,k,1) == 1) then
                            vs = interpolate_velocity(i, j, k, v%x%f, C, 1)
                            e_u = abs(v%x%f(i,j,k) - vs)
                            if (e_u > e_u_max(l)) e_u_max(l) = e_u
                        endif

                        if (ibm_index(i,j,k,2) == 1) then
                            vs = interpolate_velocity(i, j, k, v%y%f, C, 2)
                            e_v = abs(v%y%f(i,j,k) - vs)
                            if (e_v > e_v_max(l)) e_v_max(l) = e_v    
                        endif                       
                    
#if DIM==3
                        if (ibm_index(i,j,k,3) == 1) then
                            vs = interpolate_velocity(i, j, k, v%z%f, C, 3)
                            e_w = abs(v%z%f(i,j,k) - vs)
                            if (e_w > e_w_max(l)) e_w_max(l) = e_w    
                        endif    
#endif

                    end do
                end do
            end do

        end do circle_position

        ! Print maximum error for the current resolution
#if DIM==3
        write(out_id, *) Nx, maxval(e_u_max), maxval(e_v_max), maxval(e_w_max)
#else
        write(out_id, *) Nx, maxval(e_u_max), maxval(e_v_max)
#endif

        ! free memory   
        call v%destroy()
        call base_grid%destroy()

        call destroy_ibm

    end do resolution_loop

   ! Finalize the simulation
   call MPI_FINALIZE(ierror)

contains

    !=======================================================================================
    subroutine init_velocity(vi)

        use constants, only : pi

        type(vector), intent(inout) :: vi

        integer :: i, j, k
        real(dp) :: x, y
#if DIM==3
        real(dp) :: z
#endif

        do k = base_grid%lo(3),base_grid%hi(3)
            do j = base_grid%lo(2),base_grid%hi(2)
                do i = base_grid%lo(1),base_grid%hi(1)
                    x = float(i)*base_grid%delta
                    y = (float(j) - 0.5_dp)*base_grid%delta
                    vi%x%f(i,j,k) = -cos(2*pi*x)*sin(2*pi*y)
                    x = (float(i) - 0.5_dp)*base_grid%delta
                    y = float(j)*base_grid%delta
                    vi%y%f(i,j,k) = sin(2*pi*x)*cos(2*pi*y)
#if DIM==3
                    x = (float(i) - 0.5_dp)*base_grid%delta
                    z = float(k)*base_grid%delta
                    vi%z%f(i,j,k) = cos(2*pi*x)*sin(2*pi*z)
#endif
                end do
            end do
        end do

        ! Update halos and bc
        call v%apply_bc()

    end subroutine init_velocity
    !=======================================================================================

end program main
