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
    use navier_stokes        , only : v, set_timestep, viscosity
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
    C = circle(X = [0.2_dp, 0.2_dp, 0.0_dp], R = radius)
    allocate(Eulerian_Solid_list(1))
    Eulerian_solid_list(1)%pS => C
    ! Set the viscoisty
    viscosity = Umean*L/Re
    
    ! Initialize the solver
    call init_solver
    call init_ibm
    call init_eulerian_ibm(Eulerian_solid_list)
    step = 0
    time = 0.0_dp
    call set_timestep(dt, Umax)

    ! Set the inflow boundary condition
    do i = base_grid%lo(1),base_grid%hi(1)
        v%y%bc%b(i,:) = 4.0_dp*Umax*base_grid%x(i)*(Lx - base_grid%x(i))/Lx**2
    end do

    ! Open output file
    write(sn,'(I0.3)') Nx
    outfile = 'out_'//sn//'.txt'
    open(newunit = out_id, file = outfile)

    ! Print setup json file
    call print_setup_json(dt)
    
    ! Print initial fields
    call save_fields(0)
    
    print *, v%x%bc%type_r
    !==== Start Time loop ===================================================================
    time_loop: do while (time < 3.0_dp)

        step = step + 1
        time = time + dt

        ! Advance in time the solution
        call advance_solution(step, dt)

        ! Advance solver status to log file
        call print_solver_status(stdout, step, time, dt)

        ! Output fields
        if (mod(step,1) == 0) call save_fields(step)

        ! Output force coefficients per unit lenght in z direction
        block 
            real(dp) :: Cd, Cl
            if (mod(step,10) == 0) then
                Cl = 2.0_dp*Fe%x%integral()/Umean**2/L/Lz
                Cd = 2.0_dp*Fe%y%integral()/Umean**2/L/Lz
                if (base_grid%rank == 0) write(out_id,*) time, Cd, Cl 
            endif
        end block
   
    end do time_loop

    ! free memory
    call destroy_solver

    ! Finalize the simulation
    call MPI_FINALIZE(ierror)

end program laminar_cylinder
