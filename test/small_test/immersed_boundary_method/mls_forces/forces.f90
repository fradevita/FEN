program forces

    use mpi
    use precision_mod
    use global_mod           , only : ierror, stagger, pi
    use grid_mod
    use lagrangian_solid_mod , only : solid
    use lagrangian_ibm_mod   , only : compute_hydrodynamic_loads
    use scalar_mod
    use vector_mod
    use utils_mod

    implicit none

    integer, parameter  :: ns = 100
    integer             :: r, N, l, fid1, fid2, s
    real(dp)            :: Fv(2), Fp(2), theta, e_Fvx, e_Fvy, e_Fpx, e_Fpy, shift
    real(dp)            :: e_Fvx_max(ns), e_Fvy_max(ns), e_Fpx_max(ns), e_Fpy_max(ns)
    character(len=3)    :: args
    type(grid)          :: comp_grid
    type(solid), target :: C, C1
    type(scalar)        :: p, rho, mu
    type(vector)        :: v
    
    call mpi_init(ierror)

    call get_command_argument(1, args)
    read(args,'(I03)') r
    N = 2**r
    call comp_grid%setup(N, N, 1, 2.0_dp*pi, 2.0_dp*pi, 2.0*pi/real(N, dp), [0.0_dp, 0.0_dp, 0.0_dp], 4, 1)    
    
    call p%allocate(comp_grid, 1)
    call mu%allocate(comp_grid, 1)
    call rho%allocate(comp_grid, 1)
    call v%allocate(comp_grid, 1)

    ! init fields
    call init_fields
    rho%f = 1.0_dp
    mu%f  = 1.0_dp
    call p%update_ghost_nodes
    call v%update_ghost_nodes

    ! Initialize max error to zero
    e_Fvx_max = 0.0_dp
    e_Fvy_max = 0.0_dp
    e_Fpx_max = 0.0_dp
    e_Fpy_max = 0.0_dp

    do s = 1,ns
        call C%create('mesh.txt', 'C')

        ! Make a copy of C
        C1 = C

        ! Initialise random number generator.
        call random_seed(size=r)               ! Get size of seed array.
        call random_seed(put=urandom_seed(r))  ! Put seed array into PRNG.

        ! Generate and output values.
        call random_number(shift)
        call mpi_allreduce(mpi_in_place,shift,1,mpi_real8,mpi_sum,mpi_comm_world,ierror)
        shift = shift/4

        ! Apply the shift
        do l = 1,C%number_of_mass_points
            C%mass_points(l)%X(1) = C%mass_points(l)%X(1) + shift
        end do
        call C%update_lagrangian_markers()
        call C%update_center_of_mass()
        !C%tra(1) = shift
        !call C%rigid_body_motion

        do l = 1,C%number_of_edges
            C%edges(l)%C%Fv = 0.0_dp
            C%edges(l)%C%Fp = 0.0_dp
        end do

        call compute_hydrodynamic_loads(C, v, p, mu, rho, [0.0_dp, 0.0_dp, 0.0_dp])

        open(newunit = fid1, file = 'local_forces.csv')
        write(fid1,'(A33)') 'theta,x,y,Fvx,Fvy,Fpx,Fpy,nx,ny,l'
        open(newunit = fid2, file = 'ref_forces.csv')
        write(fid2,'(A25)') 'theta,x,y,Fvx,Fvy,Fpx,Fpy'

        do l = 1,C%number_of_edges
            theta = atan2(C%edges(l)%C%X(2) - C%center_of_mass%X(2), C%edges(l)%C%X(1) - C%center_of_mass%X(1)) + pi
            if (comp_grid%rank == 0) write(fid1,'(*(E16.8,:,","))') theta, C%edges(l)%C%X(1:2), C%edges(l)%C%Fv(1:2), &
                                                                           C%edges(l)%C%Fp(1:2), C%edges(l)%n, C%edges(l)%l
            call reference_forces(C%edges(l), Fv, Fp)
            if (comp_grid%rank == 0) write(fid2,'(*(E16.8,:,","))') theta, C%edges(l)%C%X(1:2), Fv, Fp
            e_Fvx = abs(C%edges(l)%C%Fv(1) - Fv(1))
            e_Fvy = abs(C%edges(l)%C%Fv(2) - Fv(2))
            e_Fpx = abs(C%edges(l)%C%Fp(1) - Fp(1))
            e_Fpy = abs(C%edges(l)%C%Fp(2) - Fp(2))
            if (e_Fvx > e_Fvx_max(s)) e_Fvx_max(s) = e_Fvx
            if (e_Fvy > e_Fvy_max(s)) e_Fvy_max(s) = e_Fvy
            if (e_Fpx > e_Fpx_max(s)) e_Fpx_max(s) = e_Fpx
            if (e_Fpy > e_Fpy_max(s)) e_Fpy_max(s) = e_Fpy
        end do
    
        call mpi_allreduce(mpi_in_place,e_Fvx_max,1,mpi_real8,mpi_max,mpi_comm_world,ierror)
        call mpi_allreduce(mpi_in_place,e_Fvy_max,1,mpi_real8,mpi_max,mpi_comm_world,ierror)
        call mpi_allreduce(mpi_in_place,e_Fpx_max,1,mpi_real8,mpi_max,mpi_comm_world,ierror)
        call mpi_allreduce(mpi_in_place,e_Fpy_max,1,mpi_real8,mpi_max,mpi_comm_world,ierror)
   
        call C%destroy()
    end do
    if (comp_grid%rank == 0) then   
        print *, N, maxval(e_Fvx_max), maxval(e_Fvy_max), maxval(e_Fpx_max), maxval(e_Fpy_max) 
    endif

    call mpi_finalize(ierror)

contains

    !==============================================================================================
    subroutine init_fields

        integer  :: i, j, k, p_id, vx_id, vy_id
        real(dp) :: x, y

        open(newunit = p_id, file = 'p.dat')
        open(newunit = vx_id, file = 'vx.dat')
        open(newunit = vy_id, file = 'vy.dat')

        do k = comp_grid%lo(3),comp_grid%hi(3)
            do j = comp_grid%lo(2),comp_grid%hi(2)
                do i = comp_grid%lo(1),comp_grid%hi(1)
                    x = comp_grid%x(i) + stagger(1,0)*comp_grid%delta
                    y = comp_grid%y(j) + stagger(2,0)*comp_grid%delta
                    p%f(i,j,k)   = test_field(x, y)
                    write(p_id,*) x, y, p%f(i,j,k)
                    x = comp_grid%x(i) + stagger(1,1)*comp_grid%delta
                    y = comp_grid%y(j) + stagger(2,1)*comp_grid%delta
                    v%x%f(i,j,k) = test_field(x, y)
                    write(vx_id,*) x, y, v%x%f(i,j,k)
                    x = comp_grid%x(i) + stagger(1,2)*comp_grid%delta
                    y = comp_grid%y(j) + stagger(2,2)*comp_grid%delta
                    v%y%f(i,j,k) = test_field(x, y)
                    write(vy_id,*) x, y, v%y%f(i,j,k)
                end do
                write(p_id,*) ''
                write(vx_id,*) ''
                write(vy_id,*) ''
            end do
        end do

    end subroutine
    !==============================================================================================

    !==============================================================================================
    real(dp)  function test_field(x, y)

        real(dp), intent(in) :: x, y

        test_field = sin(x) + cos(y)

    end function
    !==============================================================================================

    !==============================================================================================
    real(dp)  function test_field_dx(x, y)

        real(dp), intent(in) :: x, y

        test_field_dx = cos(x)

    end function
    !==============================================================================================

    !==============================================================================================
    real(dp)  function test_field_dy(x, y)

        real(dp), intent(in) :: x, y

        test_field_dy = -sin(y)

    end function
    !==============================================================================================
        
    !==============================================================================================
    subroutine reference_forces(l_edge, Fvl, Fpl)

        use lagrangian_solid_mod , only : edge

        ! In/Out variables
        type(edge), intent(in   ) :: l_edge
        real(dp)  , intent(inout) :: Fvl(2), Fpl(2)

        ! Local variables
        real(dp) :: tau11, tau12, tau22, pl

        tau11 = 2.0_dp*test_field_dx(l_edge%C%X(1), l_edge%C%X(2))
        tau22 = 2.0_dp*test_field_dy(l_edge%C%X(1), l_edge%C%X(2))
        tau12 = test_field_dx(l_edge%C%X(1), l_edge%C%X(2)) + test_field_dy(l_edge%C%X(1), l_edge%C%X(2))
        Fvl(1) = (tau11*l_edge%n(1) + tau12*l_edge%n(2))*l_edge%l
        Fvl(2) = (tau22*l_edge%n(2) + tau12*l_edge%n(1))*l_edge%l

        pl = test_field(l_edge%C%X(1), l_edge%C%X(2))
        Fpl(1) = -pl*l_edge%n(1)*l_edge%l
        Fpl(2) = -pl*l_edge%n(2)*l_edge%l

    end subroutine
    !==============================================================================================
    
end program
