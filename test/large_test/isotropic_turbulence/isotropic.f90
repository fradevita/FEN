program isotropic

    ! Test case for the forces isotropic turbulance, comparison with data from basilisk
    ! and hit3d.

    use mpi
    use precision_mod    , only : dp
    use global_mod       , only : ierror, pi
    use grid_mod
    use solver_mod
    use navier_stokes_mod , only : viscosity, set_timestep, v, constant_CFL, CFL, S
    use IO_mod

    implicit none

    ! Variables
    integer           :: Nx, Ny, Nz, step, oid
    real(dp)          :: Lx, Ly, Lz, time, dt
    type(grid)        :: comp_grid
    type(bc_type)     :: bc(6)
    
    ! Initialize MPI
    call mpi_init(ierror)

    ! The domain is a squared box of size 2 pi
    Lx = 2.0_dp*pi
    Ly = 2.0_dp*pi
    Lz = 2.0_dp*pi

    ! Set the resolution
    Nx = 128
    Ny = 128
    Nz = 128

    ! Set boundary conditions
    bc(1)%s = 'Periodic'
    bc(2)%s = 'Periodic'
    bc(3)%s = 'Periodic'
    bc(4)%s = 'Periodic'
    bc(5)%s = 'Periodic'
    bc(6)%s = 'Periodic'
    
    ! Create the grid
    call comp_grid%setup(Nx, Ny, Nz, Lx, Ly, Lz, [0.0_dp, 0.0_dp, 0.0_dp], 2, 4, bc)
    
    ! Set the viscosity
    viscosity = 1.0e-2_dp

    ! Initialize the solver
    call init_solver(comp_grid)
    step = 0
    time = 0.0_dp

    ! Compute the timestep
    call set_timestep(comp_grid, dt, 1.0_dp)
    constant_CFL = .true.
    CFL = 0.9_dp
    
    ! Initial condition
    call init_fields

    call save_fields(step)

    ! Open output file
    open(newunit = oid, file = 'output.txt')

    !==== Start Time loop ===================================================================
    time_loop: do while (time <= 300.0_dp)

        step = step + 1
        time = time + dt

        ! Apply the linear forcing
        call forcing()

        ! Advance in time the solution
        call advance_solution(comp_grid, step, dt)

        ! Advance solver status to log file
        call print_solver_status(stdout, step, time, dt)

        ! Compute output
        call output
        if (mod(step,1000)==0) call save_fields(step)

    end do time_loop

contains

    !==============================================================================================
    subroutine init_fields

        integer :: i, j, k

        ! Set the ABC flow
        do k = comp_grid%lo(3),comp_grid%hi(3)
        do j = comp_grid%lo(2),comp_grid%hi(2)
            do i = comp_grid%lo(1),comp_grid%hi(1)
                v%x%f(i,j,k) = cos(comp_grid%y(j)) + sin(comp_grid%z(k)) + rand()*1.0e-2_dp
                v%y%f(i,j,k) = sin(comp_grid%x(i)) + cos(comp_grid%z(k)) + rand()*1.0e-2_dp
                v%z%f(i,j,k) = cos(comp_grid%x(i)) + sin(comp_grid%y(j)) + rand()*1.0e-2_dp
            end do
        end do
        end do

        ! Update halos and bc
        call v%update_ghost_nodes()

    end subroutine init_fields
    !==============================================================================================

    !==============================================================================================
    subroutine forcing

        integer  :: i, j, k
        real(dp) :: ubar, vbar, wbar 

        ! Compute average velocity in each direction
        ubar = 0.0_dp
        vbar = 0.0_dp
        wbar = 0.0_dp
        do k = comp_grid%lo(3),comp_grid%hi(3)
        do j = comp_grid%lo(2),comp_grid%hi(2)
            do i = comp_grid%lo(1),comp_grid%hi(1)
                ubar = ubar + v%x%f(i,j,k)
                vbar = vbar + v%y%f(i,j,k)
                wbar = wbar + v%z%f(i,j,k)
            end do
        end do
        end do

        call mpi_allreduce(mpi_in_place,ubar,1,mpi_real8,mpi_sum,mpi_comm_world,ierror)
        call mpi_allreduce(mpi_in_place,vbar,1,mpi_real8,mpi_sum,mpi_comm_world,ierror)
        call mpi_allreduce(mpi_in_place,wbar,1,mpi_real8,mpi_sum,mpi_comm_world,ierror)
        
        ubar = ubar/(comp_grid%nx*comp_grid%ny*comp_grid%nz)
        vbar = vbar/(comp_grid%nx*comp_grid%ny*comp_grid%nz)
        wbar = wbar/(comp_grid%nx*comp_grid%ny*comp_grid%nz)
    
        ! Force the velocity field
        do k = comp_grid%lo(3),comp_grid%hi(3)
        do j = comp_grid%lo(2),comp_grid%hi(2)
            do i = comp_grid%lo(1),comp_grid%hi(1)
                S%x%f(i,j,k) = 0.1_dp*(v%x%f(i,j,k) - ubar)
                S%y%f(i,j,k) = 0.1_dp*(v%y%f(i,j,k) - vbar)
                S%z%f(i,j,k) = 0.1_dp*(v%z%f(i,j,k) - wbar)
            end do
        end do
        end do

    end subroutine forcing
    !==============================================================================================
        
    !==============================================================================================
    subroutine output

        ! Local variables
        integer  :: i, j, k
        real(dp) :: ubar, vbar, wbar, ke
        
        ! Compute average velocity in each direction
        ubar = 0.0_dp
        vbar = 0.0_dp
        wbar = 0.0_dp
        do k = comp_grid%lo(3),comp_grid%hi(3)
        do j = comp_grid%lo(2),comp_grid%hi(2)
            do i = comp_grid%lo(1),comp_grid%hi(1)
                ubar = ubar + v%x%f(i,j,k)
                vbar = vbar + v%y%f(i,j,k)
                wbar = wbar + v%z%f(i,j,k)
            end do
        end do
        end do

        call mpi_allreduce(mpi_in_place,ubar,1,mpi_real8,mpi_sum,mpi_comm_world,ierror)
        call mpi_allreduce(mpi_in_place,vbar,1,mpi_real8,mpi_sum,mpi_comm_world,ierror)
        call mpi_allreduce(mpi_in_place,wbar,1,mpi_real8,mpi_sum,mpi_comm_world,ierror)
        
        ubar = ubar/(comp_grid%nx*comp_grid%ny*comp_grid%nz)
        vbar = vbar/(comp_grid%nx*comp_grid%ny*comp_grid%nz)
        wbar = wbar/(comp_grid%nx*comp_grid%ny*comp_grid%nz)

        ! Compute kinetic energy and viscous dissipation
        ke = 0.0_dp
        do k = comp_grid%lo(3),comp_grid%hi(3)
        do j = comp_grid%lo(2),comp_grid%hi(2)
            do i = comp_grid%lo(1),comp_grid%hi(1)
                ke = ke + (v%x%f(i,j,k) - ubar)**2 + (v%y%f(i,j,k) - vbar)**2 + &
                        (v%z%f(i,j,k) - wbar)**2
            end do
        end do
        end do

        call mpi_allreduce(mpi_in_place,ke,1,mpi_real8,mpi_sum,mpi_comm_world,ierror)

        if (comp_grid%rank == 0) then
            write(oid,*) time, 0.5_dp*ke/float(nx*ny*nz) 
            flush(oid)
        endif
        
    end subroutine output
    !==============================================================================================
  
end program isotropic
