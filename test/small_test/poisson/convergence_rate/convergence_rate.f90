program test_Poisson_FDS_2D

    ! Thest the Fast Direct Solver for the solution of the Poisson equation.

    use mpi
    use global_mod    , only : pi, ierror
    use precision_mod , only : dp
    use grid_mod      , only : grid, bc_type
    use scalar_mod
    use poisson_mod

    implicit none

    integer          :: prow, pcol, n, Nx, Ny, Nz, i, j, k
    real(dp)         :: Lx, Ly, Lz, origin(3), e, emax
    type(grid)       :: G
    type(scalar)     :: phi
#if DIM==3
    type(bc_type)    :: bc(6)
#else
    type(bc_type)    :: bc(4)
#endif
    character(len=1) :: test_case

    ! Initialize MPI
    call mpi_init(ierror)

    ! Use only 1 rank
    prow = 1
    pcol = 1

    call get_command_argument(1, test_case)

    ! The test accept as input argument the case
    ! 2D
    ! case 1: test solver nn
    ! case 2: test solver pp
    ! case 3: test solver pn
    ! 3D
    ! case 1: test solver ppp
    ! case 2: test solver ppn
    ! case 3: test solver npn
    ! case 4: test solver nnn
#if DIM==3
    bc(1)%s = 'Periodic'
    bc(2)%s = 'Periodic'
    bc(3)%s = 'Periodic'
    bc(4)%s = 'Periodic'
    select case(test_case)
    case('1')
        bc(5)%s = 'Periodic'
        bc(6)%s = 'Periodic'
        ! Open the output for the test
        open(1, file = 'error_ppp')
    case('2')
        bc(5)%s = 'Wall'
        bc(6)%s = 'Wall'
        ! Open the output for the test
        open(1, file = 'error_ppn')
    case('3')
        bc(1)%s = 'Wall'
        bc(2)%s = 'Wall'
        bc(3)%s = 'Periodic'
        bc(4)%s = 'Periodic'
        bc(5)%s = 'Wall'
        bc(6)%s = 'Wall'
        ! Open the output for the test
        open(1, file = 'error_npn')
    case('4')
        bc(1)%s = 'Wall'
        bc(2)%s = 'Wall'
        bc(3)%s = 'Wall'
        bc(4)%s = 'Wall'
        bc(5)%s = 'Wall'
        bc(6)%s = 'Wall'
        ! Open the output for the test
        open(1, file = 'error_nnn')
    end select
#else
    select case(test_case)
    case('1')
        bc(1)%s = 'Periodic'
        bc(2)%s = 'Periodic'
        bc(3)%s = 'Periodic'
        bc(4)%s = 'Periodic'
        open(1, file = 'error_pp')
    case('2')
        bc(1)%s = 'Periodic'
        bc(2)%s = 'Periodic'
        bc(3)%s = 'wall'
        bc(4)%s = 'wall'
        open(1, file = 'error_pn')
    case('3')
        bc(1)%s = 'wall'
        bc(2)%s = 'wall'
        bc(3)%s = 'wall'
        bc(4)%s = 'wall'
        open(1, file = 'error_nn')
    end select
#endif

    ! Solve Poisson equation for 5 levels of refinement
    do n = 1,5

        ! Create the grid
        Nx = 16*2**(n-1)
        Ny = 16*2**(n-1)
        Lx = 1.0_dp
        Ly = 1.0_dp
#if DIM==3
        Nz = Nx        
        Lz = Lx
#else
        Nz = 1
        Lz = Lx*float(Nz)/float(Nx)
#endif
        origin = [0.0_dp, 0.0_dp, 0.0_dp]
        call G%setup(Nx, Ny, Nz, Lx, Ly, Lz, origin, prow, pcol, bc)

        ! Create one scalar field for the RHS of the Poisson equation
        ! (note that ghos nodes are not needed for phi)
        call phi%allocate(G)

        ! Set Poisson RHS
        do k = G%lo(3), G%hi(3)
            do j = G%lo(2),G%hi(2)
                do i = G%lo(1),G%hi(1)
                    phi%f(i,j,k) = RHS(G%x(i), G%y(j), G%z(k))
                end do
            end do
        end do

        ! Initialize the Poisson Solver
        call init_poisson_solver(phi)

        ! Solve the Poisson equation
        call solve_poisson(phi)

        ! Compute the error
        emax = 0.0_dp
        do k = G%lo(3),G%hi(3)
            do j = G%lo(2),G%hi(2)
                do i = G%lo(1),G%hi(1)
                    e = abs(phi%f(i,j,k) - solution(G%x(i), G%y(j), G%z(k)))
                    if (e > emax) emax = e
                end do
            end do
        end do

        ! Write the error
        write(1,*) Nx, emax

        ! Free the memory
        call destroy_poisson_solver(phi)
        call phi%destroy()
        call G%destroy

    end do
    close(1)

    ! Finalize the simulation
    call MPI_FINALIZE(ierror)

contains 

    !===============================================================================================
    real(dp) function RHS(x, y, z)

        real(dp), intent(in) :: x, y, z

#if DIM==3
        select case(test_case)
        case('1')
            RHS = -12.0_dp*pi*pi*sin(2.0_dp*pi*x)* &
                                 cos(2.0_dp*pi*y)* &
                                 sin(2.0_dp*pi*z)
        case('2')
            RHS = -4.0_dp*pi*pi*(sin(2.0_dp*pi*x) + sin(2.0_dp*pi*y) + cos(2.0_dp*pi*z))
        case('3')
            RHS = -12.0_dp*pi*pi*cos(2.0_dp*pi*x)*cos(2.0_dp*pi*y)*cos(2.0_dp*pi*z)
        case('4')
            RHS = -12.0_dp*pi*pi*cos(2.0_dp*pi*x)*cos(2.0_dp*pi*y)*cos(2.0_dp*pi*z)
        end select
#else
        select case(test_case)
        case('1')
            RHS = -8.0_dp*pi*pi*sin(2.0_dp*pi*x)*cos(2.0_dp*pi*y)
        case('2')
            RHS = -4.0_dp*pi*pi*(sin(2.0_dp*pi*x) + cos(2.0_dp*pi*y))
        case('3')
            RHS = -8.0_dp*pi*pi*cos(2.0_dp*pi*x)*cos(2.0_dp*pi*y)
        end select
#endif
    end function
    !===============================================================================================

    !===============================================================================================
    real(dp) function solution(x, y, z)

        real(dp), intent(in) :: x, y, z

#if DIM==3
        select case(test_case)
        case('1')
            solution = sin(2.0_dp*pi*x)*cos(2.0_dp*pi*y)*sin(2.0_dp*pi*z)
        case('2')
            solution = sin(2.0_dp*pi*x) + sin(2.0_dp*pi*y) + cos(2.0_dp*pi*z)
        case('3')
            solution = cos(2.0_dp*pi*x)*cos(2.0_dp*pi*y)*cos(2.0_dp*pi*z)
        case('4')
            solution = cos(2.0_dp*pi*x)*cos(2.0_dp*pi*y)*cos(2.0_dp*pi*z)
        end select
#else
        select case(test_case)
        case('1')
            solution = sin(2.0_dp*pi*x)*cos(2.0_dp*pi*y)
        case('2')
            solution = sin(2.0_dp*pi*x) + cos(2.0_dp*pi*y)            
        case('3')
            solution = cos(2.0_dp*pi*x)*cos(2.0_dp*pi*y)
        end select
#endif
    end function
    !===============================================================================================

end program test_Poisson_FDS_2D
