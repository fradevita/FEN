program main

    use mpi
    use precision_mod 
    use global_mod              , only : ierror, myrank, pi
    use grid_mod
    use euclidean_mod           , only : distance
    use lagrangian_solid_2D_mod
    use ibm_mod
    use lagrangian_solid_mod    , only : expImp
    use navier_stokes_mod       , only : set_timestep, density, viscosity, v
    use solver_mod
    use IO_mod                  , only : stdout

    implicit none

    !**** Simulation Parameters ********************************************************************
    ! Flow parameters
    real(dp), parameter :: U = 1.0_dp    ! inflow velocity
    real(dp), parameter :: rhof = 1.0_dp ! fluid density
    real(dp), parameter :: Re = 200.0_dp ! Reynolds number = rhof*U*L/muf
 
    ! Flag parameters
    real(dp), parameter :: L = 1.0_dp       ! flag edge
    real(dp), parameter :: h = 0.01_dp*L    ! flag thicnkess
    real(dp), parameter :: gamma = 100.0_dp ! density ratio = rhos/rhof
    real(dp), parameter :: eps = 2500.0_dp  ! non-dimensional elastic modulus = Eh/(rhof*U*U*L)
    real(dp), parameter :: beta = 1.0e-4_dp ! non-dimensional bending stiffness = B/(rhof*U*U*L*L*L)

    ! Derived parameters
    real(dp), parameter :: E = eps*rhof*U*U*L/h    ! Young modulus     (from eps)
    real(dp), parameter :: B = beta*rhof*U*U*L*L*L ! Bending stiffness (from beta)
    real(dp), parameter :: muf = rhof*U*L/Re       ! Fluid viscosity   (from Re)

    ! Grid parameters
    real(dp), parameter :: Lx = 8.0_dp*L
    real(dp), parameter :: Ly = 2.0_dp*L
    real(dp), parameter :: Lz = 8.0_dp*L
    integer , parameter :: Npl = 50
    integer , parameter :: Nx = Npl*8
    integer , parameter :: Ny = Npl*2
    integer , parameter :: Nz = Npl*8
    !***********************************************************************************************

    ! Variables
    integer                           :: i, step, out_id, nbc, iABC(3), for_id
    integer, allocatable              :: ibc(:)
    real(dp)                          :: dt, time, lmean, x0, y0, z0
    type(grid)                        :: comp_grid
    type(bc_type)                     :: bc(6)
    real(dp), allocatable             :: Xbc(:,:)
    type(lagrangian_solid_2D), target :: S
    character(len= 7)                 :: sstep
    character(len=18)                 :: filename

    ! Setup MPI
    call mpi_init(ierror)
    call mpi_comm_rank(mpi_comm_world, myrank, ierror)

    !**** Define the computational grid ************************************************************
    bc(1)%s = 'Wall'
    bc(2)%s = 'Wall'
    bc(3)%s = 'Periodic'
    bc(4)%s = 'Periodic'
    bc(5)%s = 'Inflow'
    bc(6)%s = 'Outflow'
    call comp_grid%setup(Nx, Ny, Nz, Lx, Ly, Lz, [-4.0_dp*L, -1.0_dp*L, -1.0_dp*L], 2, 4, bc)
    call comp_grid%print_json()
    !***********************************************************************************************

    !**** Define the solid sphere ******************************************************************
    S%is_deformable = .true.
    S%is_open = .true.
    call S%loadFromFile('mesh.txt')
    allocate(solid_list(1))
    solid_list(1)%pS => S
    allocate(S%IM(6))

    ! Set material properties
    S%density = gamma*rhof
    S%Vol = L*L*h
    S%mass = S%density*S%Vol
    call S%setMass(S%mass)
    S%IM(1:3) = S%mass

    ! Literature relations
    call S%setStretchingConstant(E*h, 'V')
    call S%setBendingConstant(B*2.0_dp/sqrt(3.0_dp), 'cos')

    ! Select deformation solver and number of substeps
    S%advance_deformation => expImp
    S%nsubsteps = 50

    ! Set constrain
    S%apply_constraints => test_constrain

    ! Evaluate average edge size
    lmean = 0.0_dp
    do i = 1,S%number_of_edges
        lmean = lmean + S%edges(i)%l
    end do
    lmean = lmean/S%number_of_edges
    if (myrank == 0) print *, 'Average edge size: ', lmean, 'Grid spacing: ', comp_grid%delta, &
                                'Lagrangian/Eulerian grid spacing ratio: ', lmean/comp_grid%delta

    ! Find boudnary points
    nbc = 0
    do i = 1,S%number_of_mass_points
        if ( abs(S%mass_points(i)%X(1) - 0.0_dp) < 1.0e-11 .and. &
             abs(S%mass_points(i)%X(3) - 0.0_dp) < 1.0e-11) then
            nbc = nbc + 1
        end if
        
        ! Search point A
        x0 = sin(0.1_dp*pi)
        y0 = -0.5_dp 
        z0 = cos(0.1_dp*pi)
        if (distance(S%mass_points(i)%X(1:3), [x0, y0, z0]) < 1.0e-6_dp) iABC(1) = i
        
        ! Search point B
        y0 = 0.0_dp
        if (distance(S%mass_points(i)%X(1:3), [x0, y0, z0]) < 1.0e-6_dp) iABC(2) = i
        
        ! Search point C
        y0 = 0.5_dp
        if (distance(S%mass_points(i)%X(1:3), [x0, y0, z0]) < 1.0e-6_dp) iABC(3) = i
    end do

    allocate(ibc(nbc))
    allocate(Xbc(nbc,3))
    nbc = 1
    do i = 1,S%number_of_mass_points
        if ( abs(S%mass_points(i)%X(1) - 0.0_dp) < 1.0e-11 .and. &
             abs(S%mass_points(i)%X(3) - 0.0_dp) < 1.0e-11) then
            ibc(nbc) = i
            Xbc(nbc,:) = S%mass_points(i)%X(1:3)
            nbc = nbc + 1
        end if
    end do

    ! Open output file for ABC data
    open(newunit = out_id, file = 'data/ABC.csv')
    if (myrank == 0) then
        write(out_id,'(A28)') 't,XA,YA,ZA,XB,YB,ZB,XC,YC,ZC'
        write(out_id,'(*(E16.8,:,","))') 0.0_dp, S%mass_points(iABC(1))%X(1:3), &
                                                 S%mass_points(iABC(2))%X(1:3), &
                                                 S%mass_points(iABC(3))%X(1:3)
        flush(out_id)
    endif
    open(newunit = for_id, file = 'data/Forces.csv')
    if (myrank == 0) write(for_id,'(A10)') 't,Fx,Fy,Fz'
    !***********************************************************************************************

    !**** Initilize the solver *********************************************************************
    density = rhof
    viscosity = muf
    call init_solver(comp_grid)
    step = 0
    time = 0.0_dp
    call set_timestep(comp_grid, dt, U)
    dt = 0.4*dt
    call save_fields(step)
    call S%writeSTL('data/S_0000000.stl')

    ! Set BCs
    v%z%bc%front = U
    v%z%bc%left = U
    v%z%bc%right = U
    !***********************************************************************************************

    !**** Start Time loop **************************************************************************
    time_loop: do while(time < 20.0_dp)
    
        step = step + 1
        time = time + dt

        ! Advance in time the solution
        call advance_solution(comp_grid, step, dt)

        ! Advance solver status to log file
        call print_solver_status(stdout, step, time, dt)

        if (mod(step,10) == 0) then
            call S%integrate_hydrodynamic_forces()
            if (myrank == 0) then
                write(out_id,'(*(E16.8,:,","))') time, S%mass_points(iABC(1))%X(1:3), &
                                                       S%mass_points(iABC(2))%X(1:3), &
                                                       S%mass_points(iABC(3))%X(1:3)
                flush(out_id)
                write(for_id,'(*(E16.8,:,","))') time, S%center_of_mass%Fh(1:3)
                flush(for_id)
            endif
            
        endif

        if (mod(step,100) == 0) then
            call save_fields(step)
            write(sstep,'(I0.7)') step
            filename = 'data/S_'//sstep//'.stl'
            call S%writeSTL(filename)
        endif

    end do time_loop
    !***********************************************************************************************

    ! free memory
    call destroy_solver
    call comp_grid%destroy()

    call mpi_finalize(ierror)

contains

    !===============================================================================================
    subroutine test_constrain(self)

        use lagrangian_solid_mod

        class(lagrangian_solid), intent(inout), target :: self

        do i = 1,size(ibc)
            S%mass_points(ibc(i))%X(1:3) = Xbc(i,:)
            S%mass_points(ibc(i))%V(1:3) = 0.0_dp
            S%mass_points(ibc(i))%A(1:3) = 0.0_dp
        end do

    end subroutine
    !===============================================================================================

end program main
