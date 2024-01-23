program laminar_cylinder

    ! Flow around a fixed cylinder at Re = 20. Testace from
    ! http://www.mathematik.tu-dortmund.de/~featflow/en/benchmarks/cfdbenchmarking/flow/dfg_benchmark1_re20.html

    use mpi
    use global_mod          , only : ierror, myrank, pi
    use precision_mod       , only : dp
    use IO_mod              , only : stdout
    use grid_mod           
    use eulerian_circle_mod , only : Circle
    use ibm_mod                
    use eulerian_ibm_mod
    use navier_stokes_mod        , only : v, set_timestep, viscosity, p, mu, rho, g
    use solver_mod

    implicit none

    ! Parameters
    real(dp), parameter :: radius = 0.05_dp  !< cylinder radius
    real(dp), parameter :: Umax = 0.3_dp     !< maximum inflow velocity
    real(dp), parameter :: Umean = 0.2_dp    !< mean inflow velocity
    real(dp), parameter :: L = 2.0_dp*radius !< reference lenght
    real(dp), parameter :: Re = 20.0_dp      !< Reynolds number: Umean*L/nu

    ! Variables
    integer              :: Nx, Ny, Nz, step, out_id, i
    real(dp)             :: Lx, Ly, Lz, time, dt, origin(3)
    type(grid), target   :: comp_grid
    type(bc_type)        :: bc(4)
    type(Circle), target :: C
    character(len=3)     :: arg
    character(len=3)     :: sn
    character(len=11)    :: outfile
    character(len=12)    :: mesh_file

    ! Initialize MPI
    call mpi_init(ierror)
    call mpi_comm_rank(mpi_comm_world, myrank, ierror)
    
    call get_command_argument(1, arg)

    ! Set the resolution to about 24 points per radius
    read(arg,'(I03)') Nx!Nx = 64
    Ny = Nx*5
    
    ! The domain has size [0.41, 2.2]
    Lx = 0.41_dp
    Ly = Lx*5.0_dp

    ! Since 2D
    Nz = 1
    Lz = Lx*float(Nz)/float(Nx)

    ! Create the grid
    bc(1)%s = 'Wall'
    bc(2)%s = 'Wall'
    bc(3)%s = 'Inflow'
    bc(4)%s = 'Outflow'
    origin = [0.0_dp, 0.0_dp, 0.0_dp]
    call comp_grid%setup(Nx, Ny, Nz, Lx, Ly, Lz, origin, 4, 1, bc)

    ! Init the solid object
    C%R = radius
    C%name = 'C'
    C%G => comp_grid
    call C%setup()

    allocate(solid_list(1))
    solid_list(1)%pS => C
    write(sn,'(I0.3)') Nx
    mesh_file = 'mesh_'//sn//'.txt'
    call C%load_surface_points(mesh_file)

    C%center_of_mass%X(1:2) = [0.2_dp, 0.2_dp]

    ! Set the viscoisty
    viscosity = Umean*L/Re
    
    ! Initialize the solver
    call init_solver(comp_grid)
    step = 0
    time = 0.0_dp
    call set_timestep(comp_grid, dt, Umax)

    ! Set the inflow boundary condition
    do i = comp_grid%lo(1),comp_grid%hi(1)
        v%y%bc%bottom(i,:) = 4.0_dp*Umax*comp_grid%x(i)*(Lx - comp_grid%x(i))/Lx**2
    end do

    ! Open output file
    outfile = 'out_'//sn//'.txt'
    open(newunit = out_id, file = outfile)
    
    !==== Start Time loop ===================================================================
    time_loop: do while (time < 3.0_dp)

        step = step + 1
        time = time + dt

        ! Advance in time the solution
        call advance_solution(comp_grid, step, dt)

        ! Advance solver status to log file
        call print_solver_status(stdout, step, time, dt)

        ! Output force coefficients per unit lenght in z direction
        block 
            real(dp) :: Cd, Cl
            if (mod(step,10) == 0) then
                ! Evaluate forces as integral of eulerain forcing Fe
                Cl = 2.0_dp*Fe%x%integral()/Umean**2/L/Lz
                Cd = 2.0_dp*Fe%y%integral()/Umean**2/L/Lz
                ! Evalute forces with probes
                call compute_hydrodynamic_loads(C, v, p, mu, rho, g)
                if (myrank == 0) write(out_id,*) time, Cd, Cl, 2.0_dp*C%center_of_mass%Fh(2)/Umean**2/L, &
                                                               2.0_dp*C%center_of_mass%Fh(1)/Umean**2/L
            endif
        end block
   
    end do time_loop

    ! free memory
    call destroy_solver

    ! Finalize the simulation
    call MPI_FINALIZE(ierror)

end program laminar_cylinder
