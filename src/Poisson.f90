module Poisson

  use precision , only : dp
  use class_Grid, only : base_grid
  use class_Scalar

  implicit none

  include 'fftw3.f'

  interface
     subroutine Poisson_Solver(f)
       use class_Scalar
       type(scalar), intent(inout) :: f
     end subroutine Poisson_Solver
  end interface

  procedure(Poisson_Solver), pointer :: solve_Poisson => Null()

  ! Modified wavenumber
  real(dp), dimension(:), allocatable :: mwn_x
  real(dp), dimension(:), allocatable :: mwn_y
#if DIM==3
  real(dp), dimension(:), allocatable :: mwn_z
#endif

  ! Tridiagonal solver coefficients
  real(dp), dimension(:), allocatable :: a, b, c

  ! FFTW3 variables
  integer(8)                                 :: pf_x, pb_x, pf_y, pb_y
  real(dp)   , dimension(:,:,:), allocatable :: outr_x, outr_x_y, c1, d1r
  complex(dp), dimension(:,:,:), allocatable :: outc_x, outc_y, outc_x_y, d1c
#if DIM==3
  integer(8)                                 :: pf_z, pb_z
  complex(dp), dimension(:,:,:), allocatable :: outc_y_z, outc_z
#endif

  private
  public init_Poisson_Solver, solve_Poisson, destroy_Poisson_solver

contains

  !========================================================================================
  subroutine init_Poisson_Solver(phi)

    ! This subroutine setup the Poisson Solver

    use constants, only : pi
    use io       , only : print_error_message

    ! In/Out variables
    type(scalar), intent(in) :: phi

    ! Local variables
    integer  :: i, j, nx, ny, lo(3), hi(3), lo_y(3), hi_y(3)
    real(dp) :: delta
    logical  :: periodic_bc(3)
#if DIM==3
    integer :: k, nz, lo_z(3), hi_z(3)
#endif
    
    lo = base_grid%lo
    hi = base_grid%hi
    nx = base_grid%nx
    ny = base_grid%ny
    delta = base_grid%delta
    lo_y = base_grid%lo_y
    hi_y = base_grid%hi_y
    periodic_bc = base_grid%periodic_bc

#if DIM==3
    lo_z = base_grid%lo_z
    hi_z = base_grid%hi_z
    nz = base_grid%nz
#endif
    
    ! Setup transformation in x
    if (periodic_bc(1) .eqv. .true.) then
       ! FFT in the x direction

       ! Modified wave number in x
       allocate(mwn_x(nx))
       do i = 1,nx
          mwn_x(i) = 2.0_dp*(dcos(2.0_dp*pi*(i - 1.0_dp)/float(nx)) - 1.0_dp)/delta**2
       end do

       ! Allocate memory for the complex output of the FFT in the x direction
       allocate(outc_x(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))

       ! Define FFT plans
       call dfftw_plan_dft_r2c_1d(pf_x, nx, phi%f(lo(1):hi(1),lo(2),lo(3)), &
                                  outc_x(lo(1):hi(1),lo(2),lo(3)), FFTW_ESTIMATE)
       call dfftw_plan_dft_c2r_1d(pb_x, nx, outc_x(lo(1):hi(1),lo(2),lo(3)), &
                                  phi%f(lo(1):hi(1),lo(2),lo(3)), FFTW_ESTIMATE)

       ! Allocate memory for the transoposition in y of the complex output in x
       allocate(outc_x_y(lo_y(1):hi_y(1),lo_y(2):hi_y(2),lo_y(3):hi_y(3)))

    elseif (periodic_bc(1) .eqv. .false.) then

       ! DCT in the x direction

       ! Modified wave number in x
       allocate(mwn_x(nx))
       do i = 1,nx
          mwn_x(i) = 2.0_dp*(dcos(1.0_dp*pi*(i - 1.0_dp)/float(nx)) - 1.0_dp)/delta**2
       end do

       ! Allocate memory for the real output of the DCT in the x direction
       allocate(outr_x(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))

       ! Define DCT plans
       call dfftw_plan_r2r_1d(pf_x, nx, phi%f(lo(1):hi(1),lo(2),lo(3)), &
                              outr_x(lo(1):hi(1),lo(2),lo(3)), FFTW_REDFT10, FFTW_ESTIMATE)
       call dfftw_plan_r2r_1d(pb_x, nx, outr_x(lo(1):hi(1),lo(2),lo(3)), &
                              phi%f(lo(1):hi(1),lo(2),lo(3)), FFTW_REDFT01, FFTW_ESTIMATE)

       ! Allocate memory for the transoposition in y of the real output in x
       allocate(outr_x_y(lo_y(1):hi_y(1),lo_y(2):hi_y(2),lo_y(3):hi_y(3)))
    else
       call print_error_message('WARNING: wrong bc in x direction')
    endif

    ! Setup transformation in y
    if (periodic_bc(2) .eqv. .true.) then
       ! FFT in the y direction

       ! Modified wave number in y
       allocate(mwn_y(ny))
       do j = 1,ny
          mwn_y(j) = 2.0_dp*(dcos(2.0_dp*pi*(j - 1)/float(ny)) - 1.0_dp)/delta**2
       end do
       
       ! Allocate memory for the complex output of the FFT in the y direction
       allocate(outc_y(lo_y(1):hi_y(1),lo_y(2):hi_y(2),lo_y(3):hi_y(3)))

       ! Define FFT plans
       call dfftw_plan_dft_1d(pf_y, ny, outc_x_y(lo_y(1),lo_y(2):hi_y(2),lo_y(3)), &
            outc_y(lo_y(1),lo_y(2):hi_y(2),lo_y(3)), FFTW_FORWARD, FFTW_ESTIMATE)
       call dfftw_plan_dft_1d(pb_y, ny, outc_y(lo_y(1),lo_y(2):hi_y(2),lo_y(3)), &
            outc_y(lo_y(1),lo_y(2):hi_y(2),lo_y(3)), FFTW_BACKWARD, FFTW_ESTIMATE)

    elseif (periodic_bc(2) .eqv. .false.) then

       ! Tridiagonal solver in y direction
       allocate(a(ny))
       allocate(b(ny))
       allocate(c(ny))
       do j = 1,ny
          a(j) =  1.0_dp/delta**2
          b(j) = -2.0_dp/delta**2
          c(j) =  1.0_dp/delta**2
       end do
       b(1)  = b(1) + a(1)
       if (base_grid%boundary_conditions(3)%s == 'Inflow' .and. base_grid%boundary_conditions(4)%s == 'Outflow') then
          b(ny) = b(ny) - c(ny)
       else
          b(ny) = b(ny) + c(ny)
       endif
       a(1)  = 0.0_dp
       c(ny) = 0.0_dp

       ! Based on the type of transform in x, arrays are real or complex
       if (periodic_bc(1) .eqv. .true.) then
          allocate(c1(lo_y(1):hi_y(1),lo_y(2):hi_y(2),lo_y(3):hi_y(3)))
          allocate(d1c(lo_y(1):hi_y(1),lo_y(2):hi_y(2),lo_y(3):hi_y(3)))
       else
          allocate(c1(lo_y(1):hi_y(1),lo_y(2):hi_y(2),lo_y(3):hi_y(3)))
          allocate(d1r(lo_y(1):hi_y(1),lo_y(2):hi_y(2),lo_y(3):hi_y(3)))
       endif
    else
       call print_error_message('WARNING: wrong bc in y direction')
    endif

    ! Select the Poisson solver
    if ((periodic_bc(1) .eqv. .true.) .and. (periodic_bc(2) .eqv. .true.)) then
       solve_Poisson => Poisson_Solver_PP
    elseif ((periodic_bc(1) .eqv. .true.) .and. (periodic_bc(2) .eqv. .false.)) then
       solve_Poisson => Poisson_Solver_PN
    elseif ((periodic_bc(1) .eqv. .false.) .and. (periodic_bc(2) .eqv. .false.)) then
       solve_Poisson => Poisson_Solver_NN
    endif

#if DIM==3
    ! Modified wave number in z
    allocate(mwn_z(nz))
    do k = 1,nz
       mwn_z(k) = 2.0_dp*(dcos(2.0_dp*pi*(k - 1.0_dp)/float(nz)) - 1.0_dp)/delta**2
    end do

    ! Allocate memory for the transoposition in z of the complex output in y
    allocate(outc_y_z(lo_z(1):hi_z(1),lo_z(2):hi_z(2),lo_z(3):hi_z(3)))

    allocate(outc_z(lo_z(1):hi_z(1),lo_z(2):hi_z(2),lo_z(3):hi_z(3)))

    ! Define FFT plans
    call dfftw_plan_dft_1d(pf_z, nz, outc_y_z(lo_z(1),lo_z(2),lo_z(3):hi_z(3)), &
                           outc_z(lo_z(1),lo_z(2),lo_z(3):hi_z(3)), FFTW_FORWARD, FFTW_ESTIMATE)
    call dfftw_plan_dft_1d(pb_z, nz, outc_z(lo_z(1),lo_z(2),lo_z(3):hi_z(3)), &
                           outc_y_z(lo_z(1),lo_z(2),lo_z(3):hi_z(3)), FFTW_BACKWARD, FFTW_ESTIMATE)

    
    solve_Poisson => Poisson_Solver_PPP
#endif

  end subroutine init_Poisson_Solver
  !========================================================================================

  !========================================================================================
  subroutine Poisson_Solver_PN(phi)

    use decomp_2d, only : transpose_x_to_y, transpose_y_to_x
    
    ! Poisson solver with periodic bc in x and neumann in y

    use mpi

    implicit none

    ! In/Out variable
    type(scalar), intent(inout) :: phi

    ! Local variables
    integer :: i, j, k, error, lo(3), hi(3), lo_y(3), hi_y(3), nx, ny, nz
    real(dp) :: frac, mean_phi

    lo = base_grid%lo
    hi = base_grid%hi
    lo_y = base_grid%lo_y
    hi_y = base_grid%hi_y
    nx = base_grid%nx
    ny = base_grid%ny
    nz = base_grid%nz
    
    ! Perform fft in the x direction of the RHS of Poisson equation
    outc_x = 0.0_dp
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          call dfftw_execute_dft_r2c(pf_x, phi%f(lo(1):hi(1),j,k), outc_x(lo(1):hi(1),j,k))
       end do
    end do

    ! Normalize output
    outc_x = outc_x / float(nx)

    ! Transpose x -> y pencil
    call transpose_x_to_y(outc_x, outc_x_y)

    ! Solve Tridiagonal system
    ! Forward step: compute c1 and d1
    do k = lo_y(3),hi_y(3)
       do i = lo_y(1),hi_y(1)
          c1(i,lo_y(2),k) = c(lo_y(2))/(b(lo_y(2)) + mwn_x(i))
          d1c(i,lo_y(2),k) = outc_x_y(i,lo_y(2),k)/(b(lo_y(2)) + mwn_x(i))
       end do
    end do
    do k = lo_y(3),hi_y(3)
       do j = lo_y(2)+1,hi_y(2)-1
          do i = lo_y(1),hi_y(1)
             c1(i,j,k) = c(j)/(b(j) - a(j)*c1(i,j-1,k) + mwn_x(i))
             d1c(i,j,k) = (outc_x_y(i,j,k) - a(j)*d1c(i,j-1,k))/(b(j) + mwn_x(i) - a(j)*c1(i,j-1,k))
          end do
       end do
    end do
    do k = lo_y(3),hi_y(3)
       do i = lo_y(1),hi_y(1)
          frac = (b(hi_y(2)) + mwn_x(i) - a(hi_y(2))*c1(i,hi_y(2)-1,k))
          if (frac /= 0.0d0) then
             d1c(i,hi_y(2),k) = (outc_x_y(i,hi_y(2),k) - a(hi_y(2))*d1c(i,hi_y(2)-1,k))/frac
          else
             d1c(i,hi_y(2),k) = 0.0_dp
          end if
       end do
    end do

    ! Backward step: solve for x
    do k = lo_y(3),hi_y(3)
       do i = lo_y(1),hi_y(1)
          outc_x_y(i,hi_y(2),k) = d1c(i,hi_y(2),k)
       end do
    end do
    do k = lo_y(3),hi_y(3)
       do j = hi_y(2)-1,lo_y(2),-1
          do i = lo_y(1),hi_y(1)
             outc_x_y(i,j,k) = d1c(i,j,k) - c1(i,j,k)*outc_x_y(i,j+1,k)
          end do
       end do
    end do

    ! Transpose y -> x pencil
    call transpose_y_to_x(outc_x_y, outc_x)

    ! Perform inverse fft in the x direction of the RHS of Poisson equation
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          call dfftw_execute_dft_c2r(pb_x, outc_x(lo(1):hi(1),j,k), phi%f(lo(1):hi(1),j,k))
       end do
    end do

    ! Normalize to zero mean value
    mean_phi = 0.d0
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             mean_phi = mean_phi + phi%f(i,j,k)
          end do
       end do
    end do

    call mpi_allreduce(mpi_in_place,mean_phi,1,mpi_real8,mpi_sum,mpi_comm_world,error)
    phi%f(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = phi%f(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) &
         - mean_phi/float(nx*ny*nz)

  end subroutine Poisson_Solver_PN
  !========================================================================================

  !========================================================================================
  subroutine Poisson_Solver_PP(phi)

    ! Poisson solver with periodic bc in x and y
    use mpi
    use decomp_2d, only : transpose_x_to_y, transpose_y_to_x

    ! In/Out variable
    type(scalar), intent(inout) :: phi

    ! Local variables
    integer  :: i, j, k, lo(3), hi(3), lo_y(3), hi_y(3), nx, ny, nz, ierror
    real(dp) :: mean_phi

    lo = base_grid%lo
    hi = base_grid%hi
    lo_y = base_grid%lo_y
    hi_y = base_grid%hi_y
    nx = base_grid%nx
    ny = base_grid%ny
    nz = base_grid%nz
    
    ! Perform fft in the x direction of the RHS of Poisson equation
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          call dfftw_execute_dft_r2c(pf_x, phi%f(lo(1):hi(1),j,k), outc_x(lo(1):hi(1),j,k))
       end do
    end do

    ! Transpose x -> y pencil
    call transpose_x_to_y(outc_x, outc_x_y)

    ! Perform fft in the y direction of the RHS of Poisson equation
    do k = lo_y(3),hi_y(3)
       do i = lo_y(1),hi_y(1)
          call dfftw_execute_dft(pf_y, outc_x_y(i,lo_y(2):hi_y(2),k), outc_y(i,lo_y(2):hi_y(2),k))
       end do
    end do

    ! Normalize output
    outc_y = outc_y/float(nx*ny)

    ! Solve Poisson equation
    do k = lo_y(3),hi_y(3)
       do j = lo_y(2),hi_y(2)
          do i = lo_y(1),hi_y(1)
             if (mwn_x(i) + mwn_y(j) == 0.0_dp) then
                outc_y(i,j,k) = 0.0_dp
             else
                outc_y(i,j,k) = outc_y(i,j,k)/(mwn_x(i) + mwn_y(j))
             endif
          end do
       end do
    end do

    ! Perform inverse fft in the y direction of the RHS of Poisson equation
    do k = lo_y(3),hi_y(3)
       do i = lo_y(1),hi_y(1)
          call dfftw_execute_dft(pb_y, outc_y(i,lo_y(2):hi_y(2),k), outc_x_y(i,lo_y(2):hi_y(2),k))
       end do
    end do

    ! Transpose y -> x pencil
    call transpose_y_to_x(outc_x_y, outc_x)

    ! Perform inverse fft in the x direction of the RHS of Poisson equation
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          call dfftw_execute_dft_c2r(pb_x, outc_x(lo(1):hi(1),j,k), phi%f(lo(1):hi(1),j,k))
       end do
    end do

    ! Normalize to zero mean value
    mean_phi = 0.0_dp
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             mean_phi = mean_phi + phi%f(i,j,k)
          end do
       end do
    end do

    call mpi_allreduce(mpi_in_place,mean_phi,1,mpi_real8,mpi_sum,mpi_comm_world,ierror)
    phi%f = phi%f - mean_phi/float(nx*ny*nz)

  end subroutine Poisson_Solver_PP
  !========================================================================================

  !========================================================================================
  subroutine Poisson_Solver_NN(phi)

    ! Poisson solver with neumann bc in x and tridiagonal solver in y.
    
    use decomp_2d, only : transpose_x_to_y, transpose_y_to_x
    use mpi

    ! In/Out variable
    type(scalar), intent(inout) :: phi

    ! Local variables
    integer  :: i, j, k, lo(3), hi(3), nx, ny, nz, lo_y(3), hi_y(3)
    real(dp) :: frac

    lo = base_grid%lo
    hi = base_grid%hi
    lo_y = base_grid%lo_y
    hi_y = base_grid%hi_y
    nx = base_grid%nx
    ny = base_grid%ny
    nz = base_grid%nz
    
    ! Perform dct in the x direction of the RHS of Poisson equation
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          call dfftw_execute_r2r(pf_x, phi%f(lo(1):hi(1),j,k), outr_x(lo(1):hi(1),j,k))
       end do
    end do

    ! Transpose x -> y pencil
    call transpose_x_to_y(outr_x, outr_x_y)

    ! Solve Tridiagonal system
    ! Forward step: compute c1 and d1
    do k = lo_y(3),hi_y(3)
       do i = lo_y(1),hi_y(1)
          c1(i,lo_y(2),k) = c(lo_y(2))/(b(lo_y(2)) + mwn_x(i))
          d1r(i,lo_y(2),k) = outr_x_y(i,lo_y(2),k)/(b(lo_y(2)) + mwn_x(i))
       end do
    end do
    do k = lo_y(3),hi_y(3)
       do j = lo_y(2)+1,hi_y(2)-1
          do i = lo_y(1),hi_y(1)
             c1(i,j,k) = c(j)/(b(j) - a(j)*c1(i,j-1,k) + mwn_x(i))
             d1r(i,j,k) = (outr_x_y(i,j,k) - a(j)*d1r(i,j-1,k))/(b(j) + mwn_x(i) - a(j)*c1(i,j-1,k))
          end do
       end do
    end do
    do k = lo_y(3),hi_y(3)
       do i = lo_y(1),hi_y(1)
          frac = (b(hi_y(2)) + mwn_x(i) - a(hi_y(2))*c1(i,hi_y(2)-1,k))
          if (frac /= 0.0d0) then
             d1r(i,hi_y(2),k) = (outr_x_y(i,hi_y(2),k) - a(hi_y(2))*d1r(i,hi_y(2)-1,k))/frac
          else
             d1r(i,hi_y(2),k) = 0.0_dp
          end if
       end do
    end do

    ! Backward step: solve for x
    do k = lo_y(3),hi_y(3)
       do i = lo_y(1),hi_y(1)
          outr_x_y(i,hi_y(2),k) = d1r(i,hi_y(2),k)
       end do
    end do
    do k = lo_y(3),hi_y(3)
       do j = hi_y(2)-1,lo_y(2),-1
          do i = lo_y(1),hi_y(1)
             outr_x_y(i,j,k) = d1r(i,j,k) - c1(i,j,k)*outr_x_y(i,j+1,k)
          end do
       end do
    end do

    ! Transpose y -> x pencil
    call transpose_y_to_x(outr_x_y, outr_x)

    ! Perform inverse fft in the x direction of the RHS of Poisson equation
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          call dfftw_execute_r2r(pb_x, outr_x(lo(1):hi(1),j,k), phi%f(lo(1):hi(1),j,k))
       end do
    end do

    ! Normalize output
    phi%f = phi%f / float(nx*2)

  end subroutine Poisson_Solver_NN
  !========================================================================================

#if DIM==3
  !========================================================================================
  subroutine Poisson_Solver_PPP(phi)

    ! 3D Poisson solver with periodic bc in x, y and z

    use decomp_2d, only : transpose_x_to_y, transpose_y_to_x
    use decomp_2d, only : transpose_y_to_z, transpose_z_to_y

    ! In/Out variable
    type(scalar), intent(inout) :: phi

    ! Local variables
    integer :: i, j, k, lo(3), hi(3), lo_y(3), hi_y(3), lo_z(3), hi_z(3), nx, ny, nz

    lo = base_grid%lo
    hi = base_grid%hi
    lo_y = base_grid%lo_y
    hi_y = base_grid%hi_y
    lo_z = base_grid%lo_z
    hi_z = base_grid%hi_z
    nx = base_grid%nx
    ny = base_grid%ny
    nz = base_grid%nz
    
    ! Perform fft in the x direction of the RHS of Poisson equation
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          call dfftw_execute_dft_r2c(pf_x, phi%f(lo(1):hi(1),j,k), outc_x(lo(1):hi(1),j,k))
       end do
    end do

    ! Transpose x -> y pencil
    call transpose_x_to_y(outc_x, outc_x_y)

    ! Perform fft in the y direction of the RHS of Poisson equation
    do k = lo_y(3),hi_y(3)
       do i = lo_y(1),hi_y(1)
          call dfftw_execute_dft(pf_y, outc_x_y(i,lo_y(2):hi_y(2),k), outc_y(i,lo_y(2):hi_y(2),k))
       end do
    end do

    ! Transpose y -> z pencil
    call transpose_y_to_z(outc_y, outc_y_z)

    ! Perform fft in the z direction of the RHS of Poisson equation
    do j = lo_z(2),hi_z(2)
       do i = lo_z(1),hi_z(1)
          call dfftw_execute_dft(pf_z, outc_y_z(i,j,lo_z(3):hi_z(3)), outc_z(i,j,lo_z(3):hi_z(3)))
       end do
    end do

    ! Normalize output
    outc_z = outc_z/float(nx*ny*nz)

    ! Solve Poisson equation
    do k = lo_z(3),hi_z(3)
       do j = lo_z(2),hi_z(2)
          do i = lo_z(1),hi_z(1)
             if (mwn_x(i) + mwn_y(j) + mwn_z(k) == 0.0_dp) then
                outc_z(i,j,k) = 0.0_dp
             else
                outc_z(i,j,k) = outc_z(i,j,k)/(mwn_x(i) + mwn_y(j) + mwn_z(k))
             endif
          end do
       end do
    end do

    ! Perform inverse fft in the z direction of the RHS of Poisson equation
    do j = lo_z(2),hi_z(2)
       do i = lo_z(1),hi_z(1)
          call dfftw_execute_dft(pb_z, outc_z(i,j,lo_z(3):hi_z(3)), outc_y_z(i,j,lo_z(3):hi_z(3)))
       end do
    end do

    ! Transpose z -> y pencil
    call transpose_z_to_y(outc_y_z, outc_y)

    ! Perform inverse fft in the y direction of the RHS of Poisson equation
    do k = lo_y(3),hi_y(3)
       do i = lo_y(1),hi_y(1)
          call dfftw_execute_dft(pb_y, outc_y(i,lo_y(2):hi_y(2),k), outc_x_y(i,lo_y(2):hi_y(2),k))
       end do
    end do

    ! Transpose y -> x pencil
    call transpose_y_to_x(outc_x_y, outc_x)

    ! Perform inverse fft in the x direction of the RHS of Poisson equation
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          call dfftw_execute_dft_c2r(pb_x, outc_x(lo(1):hi(1),j,k), phi%f(lo(1):hi(1),j,k))
       end do
    end do

  end subroutine Poisson_Solver_PPP
  !========================================================================================
#endif  
  !========================================================================================
  subroutine destroy_Poisson_solver

    deallocate(mwn_x)
    if (base_grid%periodic_bc(1) .eqv. .true. ) deallocate(outc_x, outc_x_y)
    if (base_grid%periodic_bc(1) .eqv. .false.) deallocate(outr_x, outr_x_y)
    if (base_grid%periodic_bc(2) .eqv. .true. ) deallocate(outc_y, mwn_y)
    if (base_grid%periodic_bc(2) .eqv. .false.) then
       deallocate(a,b,c,c1)
       if (base_grid%periodic_bc(1) .eqv. .true. ) deallocate(d1c)
       if (base_grid%periodic_bc(1) .eqv. .false.) deallocate(d1r)
    endif
#if DIM==3
    deallocate(mwn_z, outc_y_z, outc_z)
#endif

  end subroutine destroy_Poisson_solver
  !========================================================================================

end module Poisson
