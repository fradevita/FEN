module hypre_solver_mod

    ! This module contains the procedures for the solution of the Poisson equation using the
    ! MultiGrid solver of the library HYPRE.

    use precision_mod, only : dp

    implicit none

    ! Define the grid, the matrix A, the solution array b, the constant array d
    ! the solver and the stencil
    integer(8) :: HYPRE_grid, A, x, b, solver, stencil

    ! Define some properties of the solver
    integer  :: verbose = 0, periodic(2) = [0, 0]
    real(dp) :: tolerance = 1.0e-12_dp

contains

    !========================================================================================
    subroutine init_hypre_solver(comp_grid)

        ! This subroutin setup the HYPRE solver for the solution of the Poisson equation
        ! with iterative methods.

        use mpi
        use grid_mod

        ! In/out variables
        type(grid), intent(in) :: comp_grid

        ! Local variables
        integer :: ierr,  entry, i, j, nentries, nvalues, nx, ny
        integer, dimension(2) :: ilower, iupper, offset
        integer, dimension(5,2) :: offsets
        integer, dimension(:), allocatable :: stencil_indices
        real(dp) :: delta
        real(dp), dimension(:), allocatable :: values

        ! Local size of the grid
        nx = comp_grid%hi(1) - comp_grid%lo(1) + 1
        ny = comp_grid%hi(2) - comp_grid%lo(2) + 1

        ! **** 1. Create the HYPRE grid object ****
        call HYPRE_StructGridCreate(mpi_comm_world, 2, HYPRE_grid, ierr)

        ! Check for periodic directions
        if (comp_grid%boundary_conditions(1)%s == 'Periodic' .and. &
            comp_grid%boundary_conditions(2)%s == 'Periodic') then
            periodic(1) = comp_grid%nx
        endif
        if (comp_grid%boundary_conditions(3)%s == 'Periodic' .and. &
            comp_grid%boundary_conditions(4)%s == 'Periodic') then
            periodic(2) = comp_grid%ny
        endif
        call HYPRE_StructGridSetPeriodic(HYPRE_grid, periodic, ierr)

        ! Each rank own one box, add the box to the grid
        ilower = [comp_grid%lo(1), comp_grid%lo(2)]
        iupper = [comp_grid%hi(1), comp_grid%hi(2)]
        call HYPRE_StructGridSetExtents(HYPRE_grid, ilower, iupper, ierr)

        ! Assemble the HYPRE grid
        call HYPRE_StructGridAssemble(HYPRE_grid, ierr)

        ! **** 2. Define the discretization stencil ****

        ! Create a 2D 5-pt stencil object
        call HYPRE_StructStencilCreate(2, 5, stencil, ierr)

        ! Define the geometry of the stencil
        offsets(1,:) = [ 0, 0]
        offsets(2,:) = [-1, 0]
        offsets(3,:) = [ 1, 0]
        offsets(4,:) = [ 0,-1]
        offsets(5,:) = [ 0, 1]
        do entry = 1,5
            offset = offsets(entry,:)
            call HYPRE_StructStencilSetElement(stencil, entry - 1, offset, ierr)
        end do

        ! **** 3. Set up a Struct Matrix ****

        ! Create an empty matrix object
        call HYPRE_StructMatrixCreate(mpi_comm_world, HYPRE_grid, stencil, A, ierr)

        ! Indicate that the matrix coefficients are ready to be set
        call HYPRE_StructMatrixInitialize(A, ierr)

        ! Set the matrix cefficients. We first set the same stencil entries
        ! for each grid point. Then we make modifications to grid points
        ! near the boundary.

        ! local nx*ny grid points each with 5 stencil entries
        ilower = [comp_grid%lo(1), comp_grid%lo(2)]
        iupper = [comp_grid%hi(1), comp_grid%hi(2)]
        nentries = 5
        nvalues = nx*ny*nentries
        allocate(values(nvalues))
        allocate(stencil_indices(nentries))
        stencil_indices = [0, 1, 2, 3, 4]
        delta = comp_grid%delta
        do i = 1,nvalues,nentries
            values(i) = 4.0_dp/(delta**2)
            do j = 1,nentries-1
                values(i+j) = -1.0_dp/(delta**2)
            end do
        end do
        call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, nentries, stencil_indices, values, ierr)
        deallocate(values, stencil_indices)

        ! Set the coefficients outside of the domain
        if (periodic(1) == 0) then
            ! Values on the left of the box
            ilower = [comp_grid%lo(1), comp_grid%lo(2)]
            iupper = [comp_grid%lo(1), comp_grid%hi(2)]
            allocate(values(ny))
            allocate(stencil_indices(1))
            do i = 1,ny
                values(i) = 0.0
            end do
            stencil_indices(1) = 1
            call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)

            ! Dirichlet
            if (comp_grid%boundary_conditions(1)%s == 'Outflow') then
                do i = 1,ny
                if (i == 1 .or. i == ny) then
                    values(i) = 6.0_dp/(delta**2)
                else
                    values(i) = 5.0_dp/(delta**2)
                end if
                end do
            end if

            ! Neumann
            if (comp_grid%boundary_conditions(1)%s == 'Wall') then
                do i = 1,ny
                if (i == 1 .or. i == ny) then
                    values(i) = 2.0_dp/(delta**2)
                else
                    values(i) = 3.0_dp/(delta**2)
                end if
                end do
            end if

            stencil_indices(1) = 0
            call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)

            ! Values on the right of the box
            ilower = [comp_grid%hi(1), comp_grid%lo(2)]
            iupper = [comp_grid%hi(1), comp_grid%hi(2)]
            do i = 1,ny
                values(i) = 0.0_dp
            end do
            stencil_indices(1) = 2
            call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)

            ! Dirichlet
            if (comp_grid%boundary_conditions(2)%s == 'Outflow') then
                do i = 1,ny
                if (i == 1 .or. i == ny) then
                    values(i) = 6.0/(delta**2)
                else
                    values(i) = 5.0/(delta**2)
                end if
                end do
            end if

            ! Neumann
            if (comp_grid%boundary_conditions(2)%s == 'Wall') then
                do i = 1,ny
                if (i == 1 .or. i == ny) then
                    values(i) = 2.0_dp/(delta**2)
                else
                    values(i) = 3.0_dp/(delta**2)
                end if
                end do
            end if

            stencil_indices(1) = 0
            call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)
            deallocate(values, stencil_indices)
        end if

        if (periodic(2) == 0) then
            ! Values below the box
            ilower = [comp_grid%lo(1), comp_grid%lo(2)]
            iupper = [comp_grid%hi(1), comp_grid%lo(2)]
            allocate(values(nx))
            allocate(stencil_indices(1))
            stencil_indices(1) = 3
            do i = 1,nx
                values(i) = 0.0_dp
            end do
            call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)

            ! Dirichlet
            if (comp_grid%boundary_conditions(3)%s == 'Outflow') then
                do i = 1,nx
                if (i == 1 .or. i == nx) then
                    values(i) = 6.0_dp/(delta**2)
                else
                    values(i) = 5.0_dp/(delta**2)
                end if
                end do
            endif

            ! Neumann
            if (comp_grid%boundary_conditions(3)%s == 'Wall') then
                do i = 1,nx
                if (i == 1 .or. i == nx) then
                    values(i) = 2.0_dp/(delta**2)
                else
                    values(i) = 3.0_dp/(delta**2)
                end if
                end do
            end if

            stencil_indices(1) = 0
            call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)

            ! Values above the box
            ilower = [comp_grid%lo(1), comp_grid%hi(2)]
            iupper = [comp_grid%hi(1), comp_grid%hi(2)]
            stencil_indices(1) = 4
            do i = 1,nx
                values(i) = 0.0_dp
            end do
            call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)

            ! Dirichlet
            if (comp_grid%boundary_conditions(4)%s == 'Outflow') then
                do i = 1,nx
                if (i == 1 .or. i == nx) then
                    values(i) = 6.0_dp/(delta**2)
                else
                    values(i) = 5.0_dp/(delta**2)
                end if
                end do
            end if

            ! Neumann
            if (comp_grid%boundary_conditions(4)%s == 'Wall') then
                do i = 1,nx
                if (i == 1 .or. i == nx) then
                    values(i) = 2.0_dp/(delta**2)
                else
                    values(i) = 3.0_dp/(delta**2)
                end if
                end do
            end if

            stencil_indices(1) = 0
            call HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 1, stencil_indices, values, ierr)
            deallocate(values, stencil_indices)
        end if

        ! This is a collective call finalizing the matrix assembly.
        ! The matrix is now ``ready to be used''
        call HYPRE_StructMatrixAssemble(A, ierr)

        ! Create an empty vector object
        call HYPRE_StructVectorCreate(mpi_comm_world, HYPRE_grid, b, ierr)
        call HYPRE_StructVectorCreate(mpi_comm_world, HYPRE_grid, x, ierr)

        ! Indicate that the vector coefficients are ready to be set
        call HYPRE_StructVectorInitialize(x, ierr)
        call HYPRE_StructVectorInitialize(b, ierr)

        ! Create a PCG solver
        call HYPRE_StructPCGCreate(mpi_comm_world, solver, ierr)

        ! Set some parameters
        call HYPRE_StructPCGSetTol(solver, tolerance, ierr)      ! convergence tolerance
        call HYPRE_StructPCGSetPrintLevel(solver, verbose, ierr) ! amount of info. printed
        !call HYPRE_StructPCGSetMaxIter(solver, 50, ierr)         ! maximum number of iteration

    end subroutine init_hypre_solver
    !=========================================================================================

    !=========================================================================================
    subroutine solve_poisson(phi)

        ! This subroutine solve the poisson equation
        ! nabal^2 f = rhs
        ! with right-hand side rhs.
        ! RHS is an array of dimension nx*ny and at the end of the subroutine
        ! it will contain the solution array f.
        ! We solve the equation in the form Ax = b
        ! so d is equal to rhs and b is the solution array f.

        use grid_mod
        use scalar_mod


        implicit none

        ! In/Out variables
        type(scalar), intent(inout) :: phi

        ! Loval variables
        integer                             :: i, j, ierr, ilower(2), iupper(2), nx, ny
        type(grid)                          :: comp_grid
        real(dp), dimension(:), allocatable :: rhs
        
        comp_grid = phi%G

        nx = comp_grid%nx
        ny = comp_grid%ny

        ! The size of RHS is the size of the computational grid
        allocate(rhs(nx*ny))

        ! Get the RHS array for phi
        do j = comp_grid%lo(2),comp_grid%hi(2)
            do i = comp_grid%lo(1),comp_grid%hi(1)
                rhs(i + (j-1)*Ny) = phi%f(i,j,1)
            end do
        end do

        ! Set the array d equal to the rhs
        ilower = [1, 1]
        iupper = [comp_grid%nx, comp_grid%ny]
        call HYPRE_StructVectorSetBoxValues(b, ilower, iupper, rhs, ierr)
        Call HYPRE_StructVectorAssemble(b, ierr)

        ! Set the solution array b to zero
        rhs = 0.0_dp
        call HYPRE_StructVectorSetBoxValues(x, ilower, iupper, rhs, ierr)
        Call HYPRE_StructVectorAssemble(x, ierr)

        ! Solve the poisson equation
        call HYPRE_StructPCGSetup(solver, A, b, x, ierr)
        call HYPRE_StructPCGSolve(solver, A, b, x, ierr)

        ! Save the solution of the Poisson equation inside the rhs array
        call HYPRE_StructVectorGetBoxValues(x, ilower, iupper, rhs, ierr)

        ! Get the RHS array for phi
        do j = comp_grid%lo(2),comp_grid%hi(2)
            do i = comp_grid%lo(1),comp_grid%hi(1)
                phi%f(i,j,:) = -rhs(i + (j-1)*ny)
            end do
        end do

        deallocate(rhs)

    end subroutine solve_poisson
    !=========================================================================================

    !=========================================================================================
    subroutine destroy_hypre_solver

        ! Free the memory allocate by the HYPRE Poisson solver.

        ! Local variables
        integer :: ierr

        call HYPRE_StructGridDestroy(HYPRE_grid, ierr)
        call HYPRE_StructStencilDestroy(stencil, ierr)
        call HYPRE_StructMatrixDestroy(A, ierr)
        call HYPRE_StructVectorDestroy(b, ierr)
        call HYPRE_StructVectorDestroy(x, ierr)
        call HYPRE_StructPCGDestroy(solver, ierr)

    end subroutine destroy_hypre_solver
    !=========================================================================================

end module hypre_solver_mod
