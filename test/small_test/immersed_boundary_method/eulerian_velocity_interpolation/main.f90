program main

    ! Program to test the velocity interpolation routine of the Eulerian Immersed Boundary method.

    use mpi
    use global_mod         , only : ierror, myrank, pi
    use precision_mod      , only : dp
    use grid_mod           
    use vector_mod
    use eulerian_solid_mod
    use test_solid_mod
    use eulerian_ibm_mod

    implicit none

    ! Variables
    integer       :: out_id, res, Nx, Ny, Nz, l, i, j, k
    real(dp)      :: Lx, Ly, Lz, origin(3), e_u, e_v, e_u_max(100), e_v_max(100), vs
#if DIM==3
    real(dp)      :: e_w, e_w_max(100)
#endif
    type(vector)       :: v
    type(grid), target :: comp_grid


    ! Create one solid of type circle, it must be target
    type(test_solid), target :: C
    type(eulerian_solid_pointer) :: solid_list(1)

    ! Initialize MPI
    call mpi_init(ierror)
    call mpi_comm_rank(mpi_comm_world, myrank, ierror)

    ! The domain is a unit squared box
    Lx = 1.0_dp
    Ly = 1.0_dp
    origin = [0.0_dp, 0.0_dp, 0.0_dp]

    C = test_solid(R = 0.25_dp)
    C%G => comp_grid
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
        call comp_grid%setup(Nx, Ny, Nz, Lx, Ly, Lz, origin, 1, 1)

        ! Set the velocity field
        call v%allocate(comp_grid, 1)
        call init_velocity(v)

        ! Perform interpolation for several random initial position of the solid body
        e_u_max = 0.0_dp
        e_v_max = 0.0_dp
#if DIM==3
        e_w_max = 0.0_dp
#endif

        call init_eulerian_ibm(solid_list, comp_grid)

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
            do k = comp_grid%lo(3),comp_grid%hi(3)
                do j = comp_grid%lo(2),comp_grid%hi(2)
                    do i = comp_grid%lo(1),comp_grid%hi(1)

                        if (ibm_index(i,j,k,1) == 1) then
                            vs = interpolate_velocity(i, j, k, v%x, C, 1)
                            e_u = abs(v%x%f(i,j,k) - vs)
                            if (e_u > e_u_max(l)) e_u_max(l) = e_u
                        endif

                        if (ibm_index(i,j,k,2) == 1) then
                            vs = interpolate_velocity(i, j, k, v%y, C, 2)
                            e_v = abs(v%y%f(i,j,k) - vs)
                            if (e_v > e_v_max(l)) e_v_max(l) = e_v    
                        endif                       
                    
#if DIM==3
                        if (ibm_index(i,j,k,3) == 1) then
                            vs = interpolate_velocity(i, j, k, v%z, C, 3)
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
        call comp_grid%destroy()

        call destroy_ibm

    end do resolution_loop

   ! Finalize the simulation
   call MPI_FINALIZE(ierror)

contains

    !==============================================================================================
    subroutine init_velocity(vi)

        type(vector), intent(inout) :: vi

        integer :: i, j, k
        real(dp) :: x, y
#if DIM==3
        real(dp) :: z
#endif

        do k = comp_grid%lo(3),comp_grid%hi(3)
            do j = comp_grid%lo(2),comp_grid%hi(2)
                do i = comp_grid%lo(1),comp_grid%hi(1)
                    x = float(i)*comp_grid%delta
                    y = (float(j) - 0.5_dp)*comp_grid%delta
                    vi%x%f(i,j,k) = -cos(2*pi*x)*sin(2*pi*y)
                    x = (float(i) - 0.5_dp)*comp_grid%delta
                    y = float(j)*comp_grid%delta
                    vi%y%f(i,j,k) = sin(2*pi*x)*cos(2*pi*y)
#if DIM==3
                    x = (float(i) - 0.5_dp)*comp_grid%delta
                    z = float(k)*comp_grid%delta
                    vi%z%f(i,j,k) = cos(2*pi*x)*sin(2*pi*z)
#endif
                end do
            end do
        end do

        ! Update halos and bc
        call v%update_ghost_nodes()

    end subroutine init_velocity
    !==============================================================================================

end program main
