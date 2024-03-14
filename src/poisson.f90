module poisson_mod

    use precision_mod, only : dp
    use scalar_mod

    implicit none

    include 'fftw3.f'

    interface
        subroutine poisson_Solver(f)
            use scalar_mod
            type(scalar), intent(inout) :: f
        end subroutine poisson_Solver
    end interface

    ! Pointer to specific solver, based on BCs.
    procedure(poisson_Solver), pointer :: solve_poisson => Null()

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
    real(dp)   , dimension(:,:,:), allocatable :: outr_y, outr_y_z
    complex(dp), dimension(:,:,:), allocatable :: outc_y_z, outc_z
#endif

    private
    public init_Poisson_Solver, solve_Poisson, destroy_Poisson_solver

contains

    !===============================================================================================
    subroutine init_poisson_solver(phi)

        ! This subroutine setup the Poisson Solver

        use global_mod, only : pi
        use IO_mod    , only : print_error_message

        ! In/Out variables
        type(scalar), intent(in) :: phi

        ! Local variables
        integer  :: i, j, nx, ny, lo(3), hi(3), lo_y(3), hi_y(3)
        real(dp) :: delta
        logical  :: periodic_bc(3)
#if DIM==3
        integer :: k, nz, lo_z(3), hi_z(3)
#endif
    
        lo = phi%G%lo
        hi = phi%G%hi
        nx = phi%G%nx
        ny = phi%G%ny
        delta = phi%G%delta
        lo_y = phi%G%lo_y
        hi_y = phi%G%hi_y
        periodic_bc = phi%G%periodic_bc

#if DIM==3
        lo_z = phi%G%lo_z
        hi_z = phi%G%hi_z
        nz = phi%G%nz
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
                mwn_y(j) = 2.0_dp*(dcos(2.0_dp*pi*(j - 1.0_dp)/float(ny)) - 1.0_dp)/delta**2
            end do
            
            ! Allocate memory for the complex output of the FFT in the y direction
            allocate(outc_y(lo_y(1):hi_y(1),lo_y(2):hi_y(2),lo_y(3):hi_y(3)))

            ! Define FFT plans
            call dfftw_plan_dft_1d(pf_y, ny, outc_x_y(lo_y(1),lo_y(2):hi_y(2),lo_y(3)), &
                    outc_y(lo_y(1),lo_y(2):hi_y(2),lo_y(3)), FFTW_FORWARD, FFTW_ESTIMATE)
            call dfftw_plan_dft_1d(pb_y, ny, outc_y(lo_y(1),lo_y(2):hi_y(2),lo_y(3)), &
                    outc_y(lo_y(1),lo_y(2):hi_y(2),lo_y(3)), FFTW_BACKWARD, FFTW_ESTIMATE)

        elseif (periodic_bc(2) .eqv. .false.) then
#if DIM==3
            ! DCT in the y direction

            ! Modified wave number in y
            allocate(mwn_y(ny))
            do j = 1,ny
                mwn_y(j) = 2.0_dp*(dcos(1.0_dp*pi*(j - 1.0_dp)/float(ny)) - 1.0_dp)/delta**2
            end do

            ! Allocate memory for the real output of the DCT in the y direction
            allocate(outr_y(lo_y(1):hi_y(1),lo_y(2):hi_y(2),lo_y(3):hi_y(3)))

            ! Define FFT plans
            call dfftw_plan_r2r_1d(pf_y, ny, outr_x_y(lo_y(1),lo_y(2):hi_y(2),lo_y(3)), &
                    outr_y(lo_y(1),lo_y(2):hi_y(2),lo_y(3)), FFTW_REDFT10, FFTW_ESTIMATE)
            call dfftw_plan_r2r_1d(pb_y, ny, outr_y(lo_y(1),lo_y(2):hi_y(2),lo_y(3)), &
                    outr_x_y(lo_y(1),lo_y(2):hi_y(2),lo_y(3)), FFTW_REDFT01, FFTW_ESTIMATE)
#else
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
            if (phi%G%boundary_conditions(3)%s == 'Inflow' .and. &
                phi%G%boundary_conditions(4)%s == 'Outflow') then
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
#endif
        else
            call print_error_message('WARNING: wrong bc in y direction')
        endif

        ! Select the Poisson solver
        if ((periodic_bc(1) .eqv. .true.) .and. (periodic_bc(2) .eqv. .true.)) then
            solve_poisson => poisson_solver_pp
        elseif ((periodic_bc(1) .eqv. .true.) .and. (periodic_bc(2) .eqv. .false.)) then
            solve_poisson => poisson_solver_pn
        elseif ((periodic_bc(1) .eqv. .false.) .and. (periodic_bc(2) .eqv. .false.)) then
            solve_poisson => poisson_solver_nn
        endif

#if DIM==3
        if (periodic_bc(2)) then
            ! Allocate memory for the transoposition in z of the complex output in y
            allocate(outc_y_z(lo_z(1):hi_z(1),lo_z(2):hi_z(2),lo_z(3):hi_z(3)))
        else
            ! Allocate memory for the transoposition in z of the real output in y
            allocate(outr_y_z(lo_z(1):hi_z(1),lo_z(2):hi_z(2),lo_z(3):hi_z(3)))
        endif

        if (periodic_bc(3)) then
            ! Modified wave number in z
            allocate(mwn_z(nz))
            do k = 1,nz
                mwn_z(k) = 2.0_dp*(dcos(2.0_dp*pi*(k - 1.0_dp)/float(nz)) - 1.0_dp)/delta**2
            end do

            ! Allocate memory for the complex output of the FFT in z
            allocate(outc_z(lo_z(1):hi_z(1),lo_z(2):hi_z(2),lo_z(3):hi_z(3)))

            ! Define FFT plans
            call dfftw_plan_dft_1d(pf_z, nz, outc_y_z(lo_z(1),lo_z(2),lo_z(3):hi_z(3)), &
                outc_z(lo_z(1),lo_z(2),lo_z(3):hi_z(3)), FFTW_FORWARD, FFTW_ESTIMATE)
            call dfftw_plan_dft_1d(pb_z, nz, outc_z(lo_z(1),lo_z(2),lo_z(3):hi_z(3)), &
                outc_y_z(lo_z(1),lo_z(2),lo_z(3):hi_z(3)), FFTW_BACKWARD, FFTW_ESTIMATE)

        elseif (periodic_bc(3) .eqv. .false.) then
            ! Tridiagonal solver in z direction
            allocate(a(nz))
            allocate(b(nz))
            allocate(c(nz))
            do k = 1,nz
                a(k) =  1.0_dp/delta**2
                b(k) = -2.0_dp/delta**2
                c(k) =  1.0_dp/delta**2
            end do
            b(1) = b(1) + a(1)
            if (phi%G%boundary_conditions(5)%s == 'Inflow' .and. &
                phi%G%boundary_conditions(6)%s == 'Outflow') then
                b(nz) = b(nz) - c(nz)
            else
                b(nz) = b(nz) + c(nz)
            endif
            a(1)  = 0.0_dp
            c(nz) = 0.0_dp

            if (periodic_bc(2)) then
                ! Allocate memory for the transoposition in z of the complex output in y
                allocate( c1(lo_z(1):hi_z(1),lo_z(2):hi_z(2),lo_z(3):hi_z(3)))
                allocate(d1c(lo_z(1):hi_z(1),lo_z(2):hi_z(2),lo_z(3):hi_z(3)))
            else
                ! Allocate memory for the transoposition in z of the real output in y
                allocate( c1(lo_z(1):hi_z(1),lo_z(2):hi_z(2),lo_z(3):hi_z(3)))
                allocate(d1r(lo_z(1):hi_z(1),lo_z(2):hi_z(2),lo_z(3):hi_z(3)))
            endif
        endif
      
        if ((periodic_bc(1) .eqv. .true.) .and. &
            (periodic_bc(2) .eqv. .true.) .and. &   
            (periodic_bc(3) .eqv. .true.)) then
            solve_poisson => poisson_solver_ppp
        elseif ((periodic_bc(1) .eqv. .true.) .and. &
                (periodic_bc(2) .eqv. .true.) .and. &   
                (periodic_bc(3) .eqv. .false.)) then
            solve_poisson => poisson_solver_ppn
        elseif ((periodic_bc(1) .eqv. .false.) .and. &
                (periodic_bc(2) .eqv. .false.) .and. &   
                (periodic_bc(3) .eqv. .false.)) then
            solve_poisson => poisson_solver_nnn
        endif
#endif

    end subroutine init_Poisson_Solver
    !===============================================================================================

    !===============================================================================================
    subroutine poisson_solver_pn(phi)
        ! Poisson solver with periodic bc in x and neumann in y

#ifdef MPI
        use mpi
        use global_mod, only : ierror
#endif
        use decomp_2d, only : transpose_x_to_y, transpose_y_to_x
    
        implicit none

        ! In/Out variable
        type(scalar), intent(inout) :: phi

        ! Local variables
        integer :: i, j, k, lo(3), hi(3), lo_y(3), hi_y(3), nx, ny, nz
        real(dp) :: frac, mean_phi

        lo = phi%G%lo
        hi = phi%G%hi
        lo_y = phi%G%lo_y
        hi_y = phi%G%hi_y
        nx = phi%G%nx
        ny = phi%G%ny
        nz = phi%G%nz
        
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
#ifdef MPI
        call mpi_allreduce(mpi_in_place,mean_phi,1,mpi_real8,mpi_sum,mpi_comm_world,ierror)
#endif
        phi%f(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = phi%f(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) &
            - mean_phi/float(nx*ny*nz)

    end subroutine poisson_solver_pn
    !===============================================================================================

    !===============================================================================================
    subroutine poisson_solver_pp(phi)

        ! Poisson solver with periodic bc in x and y
#ifdef MPI
        use mpi
        use global_mod, only : ierror
#endif
        use decomp_2d, only : transpose_x_to_y, transpose_y_to_x

        ! In/Out variable
        type(scalar), intent(inout) :: phi

        ! Local variables
        integer  :: i, j, k, lo(3), hi(3), lo_y(3), hi_y(3), nx, ny, nz
        real(dp) :: mean_phi

        lo = phi%G%lo
        hi = phi%G%hi
        lo_y = phi%G%lo_y
        hi_y = phi%G%hi_y
        nx = phi%G%nx
        ny = phi%G%ny
        nz = phi%G%nz
        
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

#ifdef MPI
        call mpi_allreduce(mpi_in_place,mean_phi,1,mpi_real8,mpi_sum,mpi_comm_world,ierror)
#endif
        phi%f = phi%f - mean_phi/float(nx*ny*nz)

    end subroutine poisson_solver_pp
    !===============================================================================================

    !===============================================================================================
    subroutine poisson_solver_nn(phi)

        ! Poisson solver with neumann bc in x and tridiagonal solver in y.
#ifdef MPI
        use mpi
#endif        
        use decomp_2d, only : transpose_x_to_y, transpose_y_to_x

        ! In/Out variable
        type(scalar), intent(inout) :: phi

        ! Local variables
        integer  :: i, j, k, lo(3), hi(3), nx, ny, nz, lo_y(3), hi_y(3)
        real(dp) :: frac

        lo = phi%G%lo
        hi = phi%G%hi
        lo_y = phi%G%lo_y
        hi_y = phi%G%hi_y
        nx = phi%G%nx
        ny = phi%G%ny
        nz = phi%G%nz
        
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
    !===============================================================================================

#if DIM==3
    !===============================================================================================
    subroutine poisson_solver_ppp(phi)

        ! 3D Poisson solver with periodic bc in x, y and z

        use decomp_2d, only : transpose_x_to_y, transpose_y_to_x
        use decomp_2d, only : transpose_y_to_z, transpose_z_to_y

        ! In/Out variable
        type(scalar), intent(inout) :: phi

        ! Local variables
        integer :: i, j, k, lo(3), hi(3), lo_y(3), hi_y(3), lo_z(3), hi_z(3), nx, ny, nz

        lo = phi%G%lo
        hi = phi%G%hi
        lo_y = phi%G%lo_y
        hi_y = phi%G%hi_y
        lo_z = phi%G%lo_z
        hi_z = phi%G%hi_z
        nx = phi%G%nx
        ny = phi%G%ny
        nz = phi%G%nz
        
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

    end subroutine
    !===============================================================================================

    !===============================================================================================
    subroutine poisson_solver_ppn(phi)

        ! 3D Poisson solver with periodic bc in x, y and z
        use mpi
        use global_mod, only : ierror
        use decomp_2d, only : transpose_x_to_y, transpose_y_to_x
        use decomp_2d, only : transpose_y_to_z, transpose_z_to_y

        ! In/Out variable
        type(scalar), intent(inout) :: phi

        ! Local variables
        integer  :: i, j, k, lo(3), hi(3), lo_y(3), hi_y(3), lo_z(3), hi_z(3), nx, ny, nz
        real(dp) :: mean_phi
        real(dp) :: factor

        ! To shorten notation
        lo = phi%G%lo
        hi = phi%G%hi
        lo_y = phi%G%lo_y
        hi_y = phi%G%hi_y
        lo_z = phi%G%lo_z
        hi_z = phi%G%hi_z
        nx = phi%G%nx
        ny = phi%G%ny
        nz = phi%G%nz

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

        ! Perform fft in the y direction of the RHS of Poisson equation
        do k = lo_y(3),hi_y(3)
            do i = lo_y(1),hi_y(1)
                call dfftw_execute_dft(pf_y, outc_x_y(i,lo_y(2):hi_y(2),k), outc_y(i,lo_y(2):hi_y(2),k))
            end do
        end do

        ! Normalize output
        outc_y = outc_y / float(ny)

        ! Transpose y -> z pencil
        call transpose_y_to_z(outc_y, outc_y_z)

        ! Solve Tridiagonal system
        ! Forward step: compute c1 and d1
        do j = lo_z(2),hi_z(2)
            do i = lo_z(1),hi_z(1)
                factor = 1.0_dp/(b(lo_z(3)) + mwn_x(i) + mwn_y(j))
                c1(i,j,lo_z(3)) = c(lo_z(3))*factor
                d1c(i,j,lo_z(3)) = outc_y_z(i,j,lo_z(3))*factor
            end do
        end do
        
        do k = lo_z(3)+1,hi_z(3)-1
            do j = lo_z(2),hi_z(2)
                do i = lo_z(1),hi_z(1)
                    factor = 1.0_dp/(b(k) + mwn_x(i) + mwn_y(j) - a(k)*c1(i,j,k-1))
                    c1(i,j,k) = c(k)*factor
                    d1c(i,j,k) = (outc_y_z(i,j,k) - a(k)*d1c(i,j,k-1))*factor
                end do
            end do
        end do

        do j = lo_z(2),hi_z(2)
            do i = lo_z(1),hi_z(1)
                factor = (b(hi_z(3)) + mwn_x(i) + mwn_y(j) - a(hi_z(3))*c1(i,j,hi_z(3)-1))
                if (factor /= 0.0d0) then
                    d1c(i,j,hi_z(3)) = (outc_y_z(i,j,hi_z(3)) - a(hi_z(3))*d1c(i,j,hi_z(3)-1))/factor
                else
                    d1c(i,j,hi_z(3)) = 0.0_dp
                end if
            end do
        end do

        ! Backward step: solve for x
        do j = lo_z(2),hi_z(2)
            do i = lo_z(1),hi_z(1)
                outc_y_z(i,j,hi_z(3)) = d1c(i,j,hi_z(3))
            end do
        end do
        do k = hi_z(3)-1,lo_z(3),-1
            do j = lo_z(2),hi_z(2)
                do i = lo_z(1),hi_z(1)
                    outc_y_z(i,j,k) = d1c(i,j,k) - c1(i,j,k)*outc_y_z(i,j,k+1)
                end do
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
        phi%f = 0.0_dp
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

#ifdef MPI
        call mpi_allreduce(mpi_in_place,mean_phi,1,mpi_real8,mpi_sum,mpi_comm_world,ierror)
#endif
        phi%f = phi%f - mean_phi/float(nx*ny*nz)

    end subroutine
    !===============================================================================================

    !===============================================================================================
    subroutine poisson_solver_nnn(phi)

        ! 3D Poisson solver with periodic bc in x, y and z
        use mpi
        use global_mod , only : ierror
        use decomp_2d  , only : transpose_x_to_y, transpose_y_to_x
        use decomp_2d  , only : transpose_y_to_z, transpose_z_to_y

        ! In/Out variable
        type(scalar), intent(inout) :: phi

        ! Local variables
        integer  :: i, j, k, lo(3), hi(3), lo_y(3), hi_y(3), lo_z(3), hi_z(3), nx, ny, nz
        real(dp) :: mean_phi
        real(dp) :: factor

        ! To shorten notation
        lo = phi%G%lo
        hi = phi%G%hi
        lo_y = phi%G%lo_y
        hi_y = phi%G%hi_y
        lo_z = phi%G%lo_z
        hi_z = phi%G%hi_z
        nx = phi%G%nx
        ny = phi%G%ny
        nz = phi%G%nz

        ! Perform DCT in the x direction of the RHS of Poisson equation
        outr_x = 0.0_dp
        do k = lo(3),hi(3)
            do j = lo(2),hi(2)
                call dfftw_execute_r2r(pf_x, phi%f(lo(1):hi(1),j,k), outr_x(lo(1):hi(1),j,k))
            end do
        end do

        ! Normalize output
        outr_x = outr_x / real(nx*2, dp)

        ! Transpose x -> y pencil
        call transpose_x_to_y(outr_x, outr_x_y)

        ! Perform DCT in the y direction of the RHS of Poisson equation
        do k = lo_y(3),hi_y(3)
            do i = lo_y(1),hi_y(1)
                call dfftw_execute_r2r(pf_y, outr_x_y(i,lo_y(2):hi_y(2),k), outr_y(i,lo_y(2):hi_y(2),k))
            end do
        end do

        ! Normalize output
        outr_y = outr_y / real(ny*2, dp)

        ! Transpose y -> z pencil
        call transpose_y_to_z(outr_y, outr_y_z)

        ! Solve Tridiagonal system
        ! Forward step: compute c1 and d1
        do j = lo_z(2),hi_z(2)
            do i = lo_z(1),hi_z(1)
                factor = 1.0_dp/(b(lo_z(3)) + mwn_x(i) + mwn_y(j))
                c1(i,j,lo_z(3)) = c(lo_z(3))*factor
                d1r(i,j,lo_z(3)) = outr_y_z(i,j,lo_z(3))*factor
            end do
        end do
        
        do k = lo_z(3)+1,hi_z(3)-1
            do j = lo_z(2),hi_z(2)
                do i = lo_z(1),hi_z(1)
                    factor = 1.0_dp/(b(k) + mwn_x(i) + mwn_y(j) - a(k)*c1(i,j,k-1))
                    c1(i,j,k) = c(k)*factor
                    d1r(i,j,k) = (outr_y_z(i,j,k) - a(k)*d1r(i,j,k-1))*factor
                end do
            end do
        end do

        do j = lo_z(2),hi_z(2)
            do i = lo_z(1),hi_z(1)
                factor = (b(hi_z(3)) + mwn_x(i) + mwn_y(j) - a(hi_z(3))*c1(i,j,hi_z(3)-1))
                if (factor /= 0.0d0) then
                    d1r(i,j,hi_z(3)) = (outr_y_z(i,j,hi_z(3)) - a(hi_z(3))*d1r(i,j,hi_z(3)-1))/factor
                else
                    d1r(i,j,hi_z(3)) = 0.0_dp
                end if
            end do
        end do

        ! Backward step: solve for x
        do j = lo_z(2),hi_z(2)
            do i = lo_z(1),hi_z(1)
                outr_y_z(i,j,hi_z(3)) = d1r(i,j,hi_z(3))
            end do
        end do
        do k = hi_z(3)-1,lo_z(3),-1
            do j = lo_z(2),hi_z(2)
                do i = lo_z(1),hi_z(1)
                    outr_y_z(i,j,k) = d1r(i,j,k) - c1(i,j,k)*outr_y_z(i,j,k+1)
                end do
            end do
        end do

        ! Transpose z -> y pencil
        call transpose_z_to_y(outr_y_z, outr_y)

        ! Perform inverse DCT in the y direction of the RHS of Poisson equation
        do k = lo_y(3),hi_y(3)
            do i = lo_y(1),hi_y(1)
                call dfftw_execute_r2r(pb_y, outr_y(i,lo_y(2):hi_y(2),k), outr_x_y(i,lo_y(2):hi_y(2),k))
            end do
        end do

        ! Transpose y -> x pencil
        call transpose_y_to_x(outr_x_y, outr_x)

        ! Perform inverse DCT in the x direction of the RHS of Poisson equation
        phi%f = 0.0_dp
        do k = lo(3),hi(3)
            do j = lo(2),hi(2)
                call dfftw_execute_r2r(pb_x, outr_x(lo(1):hi(1),j,k), phi%f(lo(1):hi(1),j,k))
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

#ifdef MPI
        call mpi_allreduce(mpi_in_place,mean_phi,1,mpi_real8,mpi_sum,mpi_comm_world,ierror)
#endif
        phi%f = phi%f - mean_phi/real(nx*ny*nz, dp)

    end subroutine
    !===============================================================================================

#endif  
    !===============================================================================================
    subroutine destroy_Poisson_solver(phi)

        type(scalar), intent(in) :: phi

        ! x dir
        if (phi%G%periodic_bc(1) .eqv. .true. ) deallocate(mwn_x, outc_x, outc_x_y)
        if (phi%G%periodic_bc(1) .eqv. .false.) deallocate(mwn_x, outr_x, outr_x_y)
        
        ! y dir
        if (phi%G%periodic_bc(2) .eqv. .true. ) deallocate(mwn_y, outc_y)
        if (phi%G%periodic_bc(2) .eqv. .false.) then
#if DIM==3
            deallocate(mwn_y, outr_y)
#else
            deallocate(a, b, c, c1)
            if (phi%G%periodic_bc(1)) then 
                deallocate(d1c)
            else
                deallocate(d1r)
            endif
#endif
        endif

        ! z dir
#if DIM==3
        if (phi%G%periodic_bc(2)) then
            deallocate(outc_y_z)
        else
            deallocate(outr_y_z)
        endif
        if (phi%G%periodic_bc(3)) then
            deallocate(mwn_z, outc_z)
        else
            deallocate(a, b, c)
            if (phi%G%periodic_bc(2)) then
                deallocate(c1, d1c)
            else
                deallocate(c1, d1r)                
            endif
        endif
#endif

    end subroutine destroy_Poisson_solver
    !===============================================================================================

end module
