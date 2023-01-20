program laminar_cylinder

    ! Flow around a fixed cylinder at Re = 20. Testace from
    ! http://www.mathematik.tu-dortmund.de/~featflow/en/benchmarks/cfdbenchmarking/flow/dfg_benchmark1_re20.html

    use mpi
    use constants            , only : pi
    use precision            , only : dp
    use io                   , only : stdout
    use class_Grid           , only : base_grid, bc_type
    use class_eulerian_circle, only : Circle
    use ibm                
    use eulerian_ibm
    use navier_stokes        , only : v, set_timestep, viscosity, p, mu, rho, g
    use solver
    use json

    implicit none

    ! Parameters
    real(dp), parameter :: radius = 0.05_dp  !< cylinder radius
    real(dp), parameter :: Umax = 0.3_dp     !< maximum inflow velocity
    real(dp), parameter :: Umean = 0.2_dp    !< mean inflow velocity
    real(dp), parameter :: L = 2.0_dp*radius !< reference lenght
    real(dp), parameter :: Re = 20.0_dp      !< Reynolds number: Umean*L/nu

    ! Variables
    integer              :: ierror, Nx, Ny, Nz, step, out_id, i
    real(dp)             :: Lx, Ly, Lz, time, dt, origin(3)
    type(bc_type)        :: bc(4)
    type(Circle), target :: C
    character(len=3)     :: sn
    character(len=11)    :: outfile
    character(len=12)    :: mesh_file

    ! Initialize MPI
    call mpi_init(ierror)
    
    ! Set the resolution to about 24 points per radius
    Nx = 64
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
    call base_grid%setup(Nx, Ny, Nz, Lx, Ly, Lz, origin, 8, 1, bc)

    ! Init the solid object
    C = circle(X = [0.2_dp, 0.2_dp, 0.0_dp], R = radius, name = 'C')
    allocate(Eulerian_Solid_list(1))
    Eulerian_solid_list(1)%pS => C
    write(sn,'(I0.3)') Nx
    mesh_file = 'mesh_'//sn//'.txt'
    call C%load_surface_points(mesh_file)

    ! Set the viscoisty
    viscosity = Umean*L/Re
    
    ! Initialize the solver
    call init_solver
    step = 0
    time = 0.0_dp
    call set_timestep(dt, Umax)

    ! Set the inflow boundary condition
    do i = base_grid%lo(1),base_grid%hi(1)
        v%y%bc%b(i,:) = 4.0_dp*Umax*base_grid%x(i)*(Lx - base_grid%x(i))/Lx**2
    end do

    ! Open output file
    outfile = 'out_'//sn//'.txt'
    open(newunit = out_id, file = outfile)

    ! Print setup json file
    call print_setup_json(dt)
    
    !==== Start Time loop ===================================================================
    time_loop: do while (time < 3.0_dp)

        step = step + 1
        time = time + dt

        ! Advance in time the solution
        call advance_solution(step, dt)

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
                if (base_grid%rank == 0) write(out_id,*) time, Cd, Cl, 2.0_dp*C%hF(2)/Umean**2/L, 2.0_dp*C%hF(1)/Umean**2/L
            endif
        end block
   
    end do time_loop

    ! free memory
    call destroy_solver

    ! Finalize the simulation
    call MPI_FINALIZE(ierror)

end program laminar_cylinder
