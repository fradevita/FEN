program forces

    ! Starting from a given pressure and velocity fields, evaluate forces
    ! at given Lagrangian locations.

    use mpi
    use precision_mod
    use global_mod         , only : ierror, stagger, pi, myrank
    use grid_mod
    use ibm_mod
    use lagrangian_ibm_mod , only : compute_hydrodynamic_loads
    use navier_stokes_mod  , only : p, v, rho, mu, g 
    use navier_stokes_mod  , only : allocate_navier_stokes_fields
    use utils_mod

    implicit none

    real(dp), parameter :: L = 2.0_dp*pi

    integer, parameter  :: ns = 100
    integer             :: r, N, i, fid1, fid2, s
    real(dp)            :: Fv(2), Fp(2), theta, local_error, shift
    real(dp)            :: error(ns, 4)
    character(len=3)    :: args
    type(grid)          :: comp_grid
    type(lagrangian_solid), target :: C, C0
    
    ! Initialize MPI 
    call mpi_init(ierror)
    call mpi_comm_rank(mpi_comm_world, myrank, ierror)

    ! The test accept as input argument the resolution
    call get_command_argument(1, args)
    read(args,'(I03)') r
    
    !**** Setup the computational grid *****************************************
    N = 2**r
    call comp_grid%setup(N, N, 1, L, L, L/N, [0.0_dp, 0.0_dp, 0.0_dp], 1, 1)    
    !***************************************************************************

    !**** Init fields **********************************************************
    call allocate_navier_stokes_fields(comp_grid)
    call init_fields
    rho%f = 1.0_dp
    mu%f  = 1.0_dp
    g = [0.0_dp, 0.0_dp, 0.0_dp]
    !***************************************************************************

    !**** Create the lagrangian mesh *******************************************
#if DIM==3
    call C%create('mesh_3D.txt', 'C')
#else
    call C%create('mesh.txt', 'C')
#endif

    ! Make a copy of C
    C0 = C
    !***************************************************************************

    ! Initialise random number generator.
    call random_seed(size=r)               ! Get size of seed array.
    call random_seed(put=urandom_seed(r))  ! Put seed array into PRNG.

    ! Initialize max error to zero
    error = 0.0_dp

    testing_loop: do s = 1,ns

        ! Generate and output values.
        call random_number(shift)
        call mpi_allreduce(mpi_in_place, shift, 1, mpi_real8, mpi_sum, &
                                mpi_comm_world,ierror)
        shift = shift/4

        ! Apply the shift
        do i = 1,C%number_of_mass_points
            C%mass_points(i)%X(1:2) = C0%mass_points(i)%X(1:2) + shift
        end do
        call C%update_lagrangian_markers()
        call C%update_center_of_mass()
 
        ! Set forces on Lagrangian nodes equal to zero
        do i = 1,C%number_of_edges
            C%edges(i)%C%Fv = 0.0_dp
            C%edges(i)%C%Fp = 0.0_dp
        end do

        ! Evaluate hydrodinamic loads on C
        call compute_hydrodynamic_loads(C, v, p, mu, rho, g)

        error(s,:) = 0.0_dp
        do i = 1,C%number_of_edges

            ! Evaluate reference forces on edge i
            call reference_forces(C%edges(i), Fv, Fp)
            
            ! Evaluate errors
            local_error = abs(C%edges(i)%C%Fv(1) - Fv(1))
            if (local_error > error(s, 1)) error(s, 1) = local_error

            local_error = abs(C%edges(i)%C%Fv(2) - Fv(2))
            if (local_error > error(s, 2)) error(s, 2) = local_error

            local_error = abs(C%edges(i)%C%Fp(1) - Fp(1))
            if (local_error > error(s, 3)) error(s, 3) = local_error

            local_error = abs(C%edges(i)%C%Fp(2) - Fp(2))
            if (local_error > error(s, 4)) error(s, 4) = local_error
        end do
    
        ! call mpi_allreduce(mpi_in_place, e_Fvx_max, 1, mpi_real8, mpi_max, &
        !             mpi_comm_world,ierror)
        ! call mpi_allreduce(mpi_in_place, e_Fvy_max, 1, mpi_real8, mpi_max, &
        !             mpi_comm_world,ierror)
        ! call mpi_allreduce(mpi_in_place, e_Fpx_max, 1, mpi_real8, mpi_max, &
        !             mpi_comm_world,ierror)
        ! call mpi_allreduce(mpi_in_place, e_Fpy_max, 1, mpi_real8, mpi_max, &
        !             mpi_comm_world,ierror)
   
    end do testing_loop

    if (myrank == 0) then   
        print *, N, maxval(error(:,1)), maxval(error(:,2)), &
                    maxval(error(:,3)), maxval(error(:,4)) 
    endif

    call mpi_finalize(ierror)

contains

    !===========================================================================
    subroutine init_fields()

        integer  :: i, j, k
        real(dp) :: x, y

        do k = comp_grid%lo(3),comp_grid%hi(3)
            do j = comp_grid%lo(2),comp_grid%hi(2)
                do i = comp_grid%lo(1),comp_grid%hi(1)
                    !**** P ****************************************************
                    x = comp_grid%x(i) + stagger(1,0)*comp_grid%delta
                    y = comp_grid%y(j) + stagger(2,0)*comp_grid%delta
                    p%f(i,j,k) = test_field(x, y)

                    !**** U ****************************************************
                    x = comp_grid%x(i) + stagger(1,1)*comp_grid%delta
                    y = comp_grid%y(j) + stagger(2,1)*comp_grid%delta
                    v%x%f(i,j,k) = test_field(x, y)
                    
                    !**** V ****************************************************
                    x = comp_grid%x(i) + stagger(1,2)*comp_grid%delta
                    y = comp_grid%y(j) + stagger(2,2)*comp_grid%delta
                    v%y%f(i,j,k) = test_field(x, y)
                end do
            end do
        end do

        call p%update_ghost_nodes()
        call v%update_ghost_nodes()

    end subroutine
    !===========================================================================

    !===========================================================================
    real(dp)  function test_field(x, y)

        real(dp), intent(in) :: x, y

        test_field = sin(x) + cos(y)

    end function
    !===========================================================================

    !===========================================================================
    real(dp)  function test_field_dx(x, y)

        real(dp), intent(in) :: x, y

        test_field_dx = cos(x)

    end function
    !===========================================================================

    !===========================================================================
    real(dp)  function test_field_dy(x, y)

        real(dp), intent(in) :: x, y

        test_field_dy = -sin(y)

    end function
    !===========================================================================
        
    !===========================================================================
    subroutine reference_forces(l_edge, Fvl, Fpl)

        use lagrangian_solid_mod , only : edge

        ! In/Out variables
        type(edge), intent(in   ) :: l_edge
        real(dp)  , intent(inout) :: Fvl(2), Fpl(2)

        ! Local variables
        real(dp) :: tau11, tau12, tau22, pl

        tau11 = 2.0_dp*test_field_dx(l_edge%C%X(1), l_edge%C%X(2))
        tau22 = 2.0_dp*test_field_dy(l_edge%C%X(1), l_edge%C%X(2))
        tau12 = test_field_dx(l_edge%C%X(1), l_edge%C%X(2)) + &
                    test_field_dy(l_edge%C%X(1), l_edge%C%X(2))
        Fvl(1) = (tau11*l_edge%n(1) + tau12*l_edge%n(2))*l_edge%l
        Fvl(2) = (tau22*l_edge%n(2) + tau12*l_edge%n(1))*l_edge%l

        pl = test_field(l_edge%C%X(1), l_edge%C%X(2))
        Fpl(1) = -pl*l_edge%n(1)*l_edge%l
        Fpl(2) = -pl*l_edge%n(2)*l_edge%l

    end subroutine
    !===========================================================================
    
end program
