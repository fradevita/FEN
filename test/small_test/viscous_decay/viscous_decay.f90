program main

  use mpi
  use precision      , only : dp
  use class_Grid     , only : base_grid, bc_type
  use constants      , only : gravity
  use multiphase     , only : rho_0, rho_1, mu_0, mu_1
  use volume_of_fluid, only : distance, vof
  use solver         , only : init_solver, advance_solution, print_solver_status, destroy_solver
  use navier_stokes  , only : set_timestep, g, v
  use io

  implicit none
  
  integer           :: ierror, Nx, Ny, Nz, step, out_id
  real(dp)          :: Lx, Ly, Lz, time, dt, Tmax, origin(3)
  type(bc_type)     :: bc(4)
  character(len=7)  :: ss
  character(len=16) :: filename

  ! Initialize the MPI library
  call mpi_init(ierror)

  ! Create the grid
  Nx = 128
  Ny = 256
  Nz = 1
  Lx = 1.0_dp
  Ly = 2.0_dp
  Lz = Lx*float(Nz)/float(Nx)
  origin = [0.0_dp, 0.0_dp, 0.0_dp]
  bc(1)%s = 'Periodic'
  bc(2)%s = 'Periodic'
  bc(3)%s = 'Wall'
  bc(4)%s = 'Wall'
  call base_grid%setup(Nx, Ny, Nz, Lx, Ly, Lz, origin, 4, 1, bc)

  ! Setup multiphase parameter
  rho_0 = 1000.0_dp
  rho_1 = rho_0/850.0_dp
  mu_0 = rho_0*Lx*sqrt(gravity*Lx)/1.0d+4
  mu_1 = mu_0*1.9d-2

  ! Set gravity
  g(2) = -gravity
  
  ! Set the initial vof from a distance function
  distance => wave
  
  ! Initialize the solver
  call init_solver
  step = 0
  time = 0.0d0
  call set_timestep(dt, 1.0_dp)
  dt = 0.1_dp*dt

  ! Initialize the velocity field
  call init_velocity

  ! Set maximum time
  Tmax = 5.0_dp
  
  ! Open file for output
  open(newunit = out_id, file = 'energy.dat')

  !==== Start Time loop ===================================================================
  do while(time < Tmax)

     step = step + 1
     time = time + dt

     ! Advance in time the solution
     call advance_solution(step, dt)

     call print_solver_status(stdout, step, time, dt)

     if (mod(step,100) == 0) then
        write(ss,'(I0.7)') step
        filename = 'data/v_x_'//ss
        call v%x%write(filename)
        filename = 'data/vof_'//ss
        call vof%write(filename)
     endif
        
     ! Compute mechanical energy of the wave
     call wave_energy()
  end do

  close(out_id)
  
  ! Free the memory allocated by the solver
  call destroy_solver()
  
  ! Finalize the simulation
  call MPI_FINALIZE(ierror)

contains

  !==============================================================================
  subroutine init_velocity

    ! Initialize the velocity field

    use constants      , only : pi
    use navier_stokes  , only : v
    use volume_of_fluid, only : vof
    
    integer  :: i, j, k
    real(dp) :: wn, omega, x, y, f

    wn = 2.0_dp*pi/Lx
    omega = sqrt(abs(g(2))*wn)

    do k = base_grid%lo(3),base_grid%hi(3)
      do j = base_grid%lo(2),base_grid%hi(2)
        do i = base_grid%lo(1),base_grid%hi(1)
          x = i*base_grid%delta
          y = (j - 0.5_dp)*base_grid%delta - Ly/2.0_dp
          f = (vof%f(i+1,j,k) + vof%f(i,j,k))*0.5_dp
          v%x%f(i,j,k) = (1.0_dp - f)*0.005*omega*exp( wn*y)*cos(wn*x) - &
                                    f*0.005*omega*exp(-wn*y)*cos(wn*x)

          x = (i - 0.5_dp)*base_grid%delta
          y = j*base_grid%delta - Ly/2.0_dp
          f = (vof%f(i,j+1,k) + vof%f(i,j,k))*0.5_dp
          v%y%f(i,j,k) = (1.0_dp - f)*0.005*omega*exp( wn*y)*sin(wn*x) + &
                               f*0.005*omega*exp(-wn*y)*sin(wn*x)
        end do
      end do
    end do

    ! Update boundary conditions
    call v%apply_bc()

  end subroutine
  !==============================================================================

  !==============================================================================
  function wave(x, y) result (d)

    use constants, only : pi

    implicit none

    ! In/Out variables
    real(dp), intent(in) :: x, y

    ! Local variables
    real(dp) :: d

    d = (y - 0.005_dp*cos(2.0_dp*pi*x/Lx) - Ly/2.0_dp)

  end function wave
  !==============================================================================

  !==============================================================================
  subroutine wave_energy()

    ! Compute mechanical wave energy

    use navier_stokes  , only : v
    use volume_of_fluid, only : vof
    
    integer  :: i, j, k
    real(dp) :: uc, vc, Ep, Ek, y

    Ek = 0.0_dp
    Ep = 0.0_dp

    do k = base_grid%lo(3),base_grid%hi(3)
      do j = base_grid%lo(2),base_grid%hi(2)
        y = (j - 0.5_dp)*base_grid%delta - Ly/2.0_dp
        do i = base_grid%lo(1),base_grid%hi(1)
          uc = 0.5_dp*(v%x%f(i,j,k) + v%x%f(i-1,j,k))
          vc = 0.5_dp*(v%y%f(i,j,k) + v%y%f(i,j-1,k))
          Ek = Ek + 0.5_dp*rho_0*(uc**2 + vc**2)*(1.0_dp - vof%f(i,j,k))
          Ep = Ep + rho_0*abs(g(2))*y*(1.0_dp - vof%f(i,j,k))
        end do
      end do
    end do

    call mpi_allreduce(mpi_in_place,Ek,1,mpi_real8,mpi_sum,mpi_comm_world,ierror)
    call mpi_allreduce(mpi_in_place,Ep,1,mpi_real8,mpi_sum,mpi_comm_world,ierror)
    
    if(base_grid%rank == 0) then
       write(out_id,*) time, Ek*base_grid%delta*base_grid%delta, Ep*base_grid%delta*base_grid%delta
    endif
    
  end subroutine wave_energy
  !==============================================================================

end program main
