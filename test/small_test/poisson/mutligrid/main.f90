program PoissonMG

  ! Test case for the solution of the Poisson equation using the HYPRE library
  ! with four periodic boundaries.
  
  use mpi
  use precision_mod    , only : dp
  use global_mod       , only : pi
  use grid_mod         , only : grid, bc_type
  use hypre_solver_mod , only : init_hypre_solver, solve_poisson, destroy_hypre_solver
  use scalar_mod

  implicit none

  integer       :: ierror, Nx, Ny, Nz, i, j, res
  real(dp)      :: Lx, Ly, Lz, origin(3), e, emax
  type(grid)    :: comp_grid
  type(bc_type) :: bc(4)
  type(scalar)  :: phi
  
  call mpi_init(ierror)

  open(1, file = 'error')

  resolution_loop: do res = 1,5
  
     Nx = 2**(3+res)
     Ny = 2**(3+res)
     Nz = 1
     Lx = 1.0_dp
     Ly = 1.0_dp
     Lz = Nz*Lx/Nx
     origin = [0.0_dp, 0.0_dp, 0.0_dp]
     bc(1)%s = 'Periodic'
     bc(2)%s = 'Periodic'
     bc(3)%s = 'Periodic'
     bc(4)%s = 'Periodic'

     call comp_grid%setup(Nx, Ny, Nz, Lx, Ly, Lz, origin, 4, 1, bc)

     ! Initialize the Poisson solver
     call init_hypre_solver(comp_grid)

     ! Set RHS
     call phi%allocate(comp_grid)

     do j = comp_grid%lo(2),comp_grid%hi(2)
        do i = comp_grid%lo(1),comp_grid%hi(1)
           phi%f(i,j,:) = -8.0_dp*pi*pi*sin(2.0_dp*pi*comp_grid%x(i))*sin(2.0_dp*pi*comp_grid%y(j))
        end do
     end do

     ! Solve the poisson equation
     call solve_poisson(phi)

     ! Get error
     emax = 0.0_dp
     do j = comp_grid%lo(2),comp_grid%hi(2)
        do i = comp_grid%lo(1),comp_grid%hi(1)
           e = abs(phi%f(i,j,1) - sin(2.0_dp*pi*comp_grid%x(i))*sin(2.0_dp*pi*comp_grid%y(j)))
           if (e > emax) emax = e
        end do
     end do
     call mpi_allreduce(mpi_in_place,emax,1,mpi_real8,mpi_max,mpi_comm_world,ierror)
     if (comp_grid%rank == 0) write(1,*) Nx, emax
     
     ! Free the solver memory
     call comp_grid%destroy()
     call phi%destroy()
     call destroy_hypre_solver

  end do resolution_loop
  
  call mpi_finalize(ierror)
  
end program PoissonMG
