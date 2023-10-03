program projection

    ! Test the projection of non-divergence velocity field into a divergence one.

    use mpi
    use global_mod    , only : ierror, myrank
    use precision_mod , only : dp
    use grid_mod
    use scalar_mod
    use vector_mod
    use poisson_mod
    use fields_mod
    use utils_mod    , only : urandom_seed

    implicit none

    integer          :: Nx, Ny, Nz, s, i, j, k
    real(dp)         :: Lx, Ly, Lz, origin(3), r
    type(grid)       :: G
    type(scalar)     :: phi, div
    type(vector)     :: v, grad_phi
    type(bc_type)    :: bc(6)
    character(len=1) :: test_case

    ! Initialize MPI
    call mpi_init(ierror)
    call mpi_comm_rank(mpi_comm_world, myrank, ierror)
    
    ! The test accept as input argument the case
    ! case 1: test solver ppp
    ! case 2: test solver ppn
    call get_command_argument(1, test_case)
    bc(1)%s = 'Periodic'
    bc(2)%s = 'Periodic'
    bc(3)%s = 'Periodic'
    bc(4)%s = 'Periodic'
    select case(test_case)
    case('1')
        bc(5)%s = 'Periodic'
        bc(6)%s = 'Periodic'
        ! Open the output for the second test
        open(1, file = 'error_ppp')
    case('2')
        bc(5)%s = 'Wall'
        bc(6)%s = 'Wall'
        ! Open the output for the second test
        open(1, file = 'error_ppn')
    end select

    ! Create the grid
    Nx = 32
    Ny = 32
    Nz = 32
    Lx = 1.0_dp
    Ly = 1.0_dp        
    Lz = 1.0_dp
    origin = [0.0_dp, 0.0_dp, 0.0_dp]
    call G%setup(Nx, Ny, Nz, Lx, Ly, Lz, origin, 1, 1, bc)

    ! Create velocity vector field, divergence scalar field
    ! projection operator scalar field and its gradient 
    ! vector field 
    call v%allocate(G, 1)
    call div%allocate(G, 1)
    call phi%allocate(G, 1)
    call grad_phi%allocate(G, 1)

    ! Adapt field BC based on the testcase
    if (test_case == '2') then
        v%x%bc%type_front = 1
        v%y%bc%type_front = 1
        v%z%bc%type_front = 1
        phi%bc%type_front = 2
        v%x%bc%type_back = 1
        v%y%bc%type_back = 1
        v%z%bc%type_back = 1
        phi%bc%type_back = 2
    endif
    
    ! Initialize the Poisson Solver
    call init_poisson_solver(phi)

    ! Repeat the test several times
    do s = 1,100

        ! Initialise random number generator.
        call random_seed(size=Nx)               ! Get size of seed array.
        call random_seed(put=urandom_seed(Nx))  ! Put seed array into PRNG.

        ! Random initial velocity field
        do k = G%lo(3),G%hi(3)
            do j = G%lo(2),G%hi(2)
                do i = G%lo(1),G%hi(1)
                    call random_number(r)
                    v%x%f(i,j,k) = r
                    call random_number(r)
                    v%y%f(i,j,k) = r
                    call random_number(r)
                    v%z%f(i,j,k) = r
                end do
            end do
        end do
        call v%update_ghost_nodes()

        ! Evaluate divergence
        call divergence(v, div)

        ! Set RHS of poisson equation equal to - div
        phi%f = -div%f

        ! Solve the Poisson equation
        call solve_poisson(phi)
        call phi%update_ghost_nodes()

        ! Project velocity field 
        call gradient(phi, grad_phi)
        v%x%f = v%x%f + grad_phi%x%f
        v%y%f = v%y%f + grad_phi%y%f
        v%z%f = v%z%f + grad_phi%z%f
        call v%update_ghost_nodes()
        
        ! Compute final divergence
        call divergence(v, div)

        ! Evaluate maximum value
        print *, s, div%max_value()

    end do

    ! Free the memory
    call destroy_poisson_solver(phi)
    call phi%destroy()
    call G%destroy

    ! Finalize the simulation
    call MPI_FINALIZE(ierror)

end program projection
