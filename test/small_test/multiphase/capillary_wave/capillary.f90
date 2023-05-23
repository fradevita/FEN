program main

    use mpi
    use precision_mod       , only : dp
    use global_mod          , only : pi
    use grid_mod
    use multiphase_mod      , only : rho_0, rho_1, mu_0, mu_1, sigma
    use volume_of_fluid_mod , only : distance, vof
    use navier_stokes_mod   , only : set_timestep
    use solver_mod          , only : init_solver, advance_solution, destroy_solver, print_solver_status
    use IO_mod
    
    implicit none

    ! Parameters
    real(dp), parameter :: a = 0.01_dp
    real(dp), parameter :: lambda = 1.0_dp
    real(dp), parameter :: rho_w = 1.0_dp
    real(dp), parameter :: wn = 2.0_dp*pi/lambda
    
    ! Variables
    integer           :: ierror, Nx, Ny, Nz, step, res, nprint
    real(dp)          :: Lx, Ly, Lz, Tmax, time, dt, origin(3)
    type(grid)        :: comp_grid
    type(bc_type)     :: bc(4)
    character(len=2)  :: sx 
    character(len=7)  :: ss
    character(len=19) :: filename

    ! Initialize the MPI library
    call mpi_init(ierror)

    open(1, file = 'timestep.dat')

    resolution_loop: do res = 1,4
  
        ! Create the grid
        Nx = 2**(2 + res)
        Ny = Nx*3
        Nz = 1
        Lx = lambda
        Ly = 3.0_dp*lambda
        Lz = Lx*float(Nz)/float(Nx)
        origin = [0.0_dp, -Ly/2.0_dp, 0.0_dp]
        bc(1)%s = 'Periodic'
        bc(2)%s = 'Periodic'
        bc(3)%s = 'Wall'
        bc(4)%s = 'Wall'
        call comp_grid%setup(Nx, Ny, Nz, Lx, Ly, Lz, origin, 1, 1, bc)
        comp_grid%name = 'grid'
        
        ! Setup the multiphase parameters
        rho_0 = rho_w
        rho_1 = rho_w
        mu_0 = 0.0182571749236_dp
        mu_1 = mu_0
        sigma = 1.0_dp

        ! Set the initial interface
        distance => wave

        ! Maximum time
        Tmax = 25.0_dp/11.1366559937_dp

        ! Initialize the solver
        call init_solver(comp_grid)
        step = 0
        time = 0.0_dp

        ! Set the timestep
        call set_timestep(comp_grid, dt, 1.0_dp)
        
        ! Write the initial vof
        write(sx,'(I0.2)') Nx
        write(ss,'(I0.7)') step
        filename = 'data_'//sx//'/vof_'//ss
        call vof%write(filename)

        ! Need to save dt for postprocessing reasons
        write(1,*) dt

        ! Set the output frequency in order to have about 80 files per resolution
        nprint = int(Tmax/80/dt)

        ! Print the json setup
        call case_setup()

        !==== Start Time loop ===================================================================
        do while(time < Tmax)

            step = step + 1
            time = time + dt
            
            ! Advance in time the solution
            call advance_solution(comp_grid, step, dt)
            
            ! Print solver status on the log
            call print_solver_status(stdout, step, time, dt)
            
            if (mod(step,nprint) == 0) then
            write(ss,'(I0.7)') step
            filename = 'data_'//sx//'/vof_'//ss
            call vof%write(filename)
            endif

        end do
        
        ! Free the memory allocated by the solver
        call destroy_solver
        call comp_grid%destroy

    end do resolution_loop

    ! Finalize the simulation
    call MPI_FINALIZE(ierror)

contains

    !=======================================================================================
    function wave(x, y) result (d)

        implicit none

        ! In/Out variables
        real(dp), intent(in) :: x, y

        ! Local variables
        real(dp) :: x1, y1, x2, y2, d
    
        x1 = x - comp_grid%delta/2.0_dp
        y1 = a*cos(wn*x1)

        x2 = x + comp_grid%delta/2.0_dp
        y2 = a*cos(wn*x2)
        
        d = -((x2 - x1)*(y1 - y) - (x1 - x)*(y2 - y1))/sqrt((x2 - x1)**2 + (y2 - y1)**2)

    end function
    !=======================================================================================

    !=======================================================================================
    subroutine case_setup()

        integer :: json_case_id

        if (comp_grid%rank == 0) then
            open(newunit = json_case_id, file = 'case.json')
            write(json_case_id,'(A1)') '{'
            write(json_case_id,'(4x,A9)') '"Case": {'
            write(json_case_id,'(8x,A10,1x,E16.8,A1)') '"lambda": ', lambda, ','
            write(json_case_id,'(8x,A8,1x,E16.8,A1)') '"rho": ', rho_w, ','
            write(json_case_id,'(8x,A10,1x,E16.8)') '"sigma": ', sigma
            write(json_case_id,'(4x,A3)') '}'
            write(json_case_id,'(A1)') '}'
            flush(json_case_id)
            close(json_case_id)
        endif

    end subroutine case_setup
    !=======================================================================================
    
end program main