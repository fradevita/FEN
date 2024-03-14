program main

    use precision_mod
    use lagrangian_solid_2D_mod
    use euclidean_mod

    implicit none

    ! Parameters
    real(dp), parameter :: l = 1.0_dp          ! Lenght of the plate [m]
    real(dp), parameter :: w = 0.1_dp          ! Width of the plate [m]
    real(dp), parameter :: t = 1.0e-3_dp       ! Thickness of the plate [m]
    real(dp), parameter :: V = l*w*t           ! Volume
    real(dp), parameter :: rho = 7850.0_dp     ! density
    real(dp), parameter :: M = rho*V           ! Plate mass [kg]
    real(dp), parameter :: E = 210.1e+10_dp    ! Young's Modulus [Pa]
    real(dp), parameter :: ni = 0.3_dp         ! Poisson's ratio []
    real(dp), parameter :: gamma = 0.0e-6_dp   ! Small damping
    real(dp), parameter :: dt = 6.0e-7_dp      ! timestep size [s]

    ! Variables
    integer                   :: n, step, nbcl1, nbcl2, nbcr, eid, tid, sstep
    integer     , allocatable :: ibcl1(:), ibcl2(:), ibcr(:)
    real(dp)                  :: deltal
    type(lagrangian_solid_2D) :: S
    character(len=8 )         :: char_step
    character(len=19)         :: outfile

    call S%loadFromFile('mesh.txt')
    call S%setStretchingConstant(1.704122401906595469e+09_dp)
    call S%setBendingConstant(1.512040229889113107e+04_dp) 
    S%is_deformable = .true.

    ! evaluate resolution
    deltal = S%mass_points(2)%X(1)

    ! Set material properties
    call S%setMass(M)

    ! Viscous damping
    S%gamma = gamma

    ! Find boundary points
    nbcl1 = 0
    nbcl2 = 0
    nbcr = 0
    do n = 1,S%number_of_mass_points
        if (abs(S%mass_points(n)%x(1) - 0.0_dp) < 1.0e-6_dp) nbcl1 = nbcl1 + 1
        if (abs(S%mass_points(n)%x(1) - deltal) < 1.0e-6_dp) nbcl2 = nbcl2 + 1
        if (abs(S%mass_points(n)%x(1) - l     ) < 1.0e-6_dp) nbcr = nbcr+ 1
    end do
    allocate(ibcl1(nbcl1))
    allocate(ibcl2(nbcl2))
    allocate(ibcr(nbcr))

    nbcl1 = 0
    nbcl2 = 0
    nbcr = 0
    do n = 1,S%number_of_mass_points
        if (abs(S%mass_points(n)%X(1) - 0.0_dp) < 1.0e-6_dp) then
            nbcl1 = nbcl1 + 1
            ibcl1(nbcl1) = n
        end if
        if (abs(S%mass_points(n)%X(1) - deltal) < 1.0e-6_dp) then
            nbcl2 = nbcl2 + 1
            ibcl2(nbcl2) = n
        endif
        if (abs(S%mass_points(n)%X(1) - l) < 1.0e-6_dp) then
            nbcr = nbcr + 1
            ibcr(nbcr) = n
        end if
    end do
    S%apply_constraints => test_constraints

    ! External force
    do n = 1,nbcr
        S%mass_points(ibcr(n))%Fs(3) = -22.5_dp/nbcr
    end do

    ! Print STL
    step = 0
    write(char_step,'(I0.8)') step
    outfile = 'data/P_'//char_step//'.stl'
    call S%writeSTL(outfile)
    outfile = 'data/P_'//char_step//'.gnu'    
    call S%writeGNU(outfile)

    open(newunit = eid, file = 'data/energy.dat')
    open(newunit = tid, file = 'data/tip.dat')

    S%nsubsteps = 1

    preload: do step = 1,2000000
        
        call S%advance(dt, [0.0_dp, 0.0_dp, 0.0_dp])
        
        if (mod(step,1000) == 0) then! Print STL
            write(char_step,'(I0.8)') step
            outfile = 'data/P_'//char_step//'.stl'
            call S%writeSTL(outfile)
            outfile = 'data/P_'//char_step//'.gnu'
            call S%writeGNU(outfile)
        end if
        
        ! Check sum of internal forces
        print *, 'PRELOAD: ', step, step*dt, S%internalForcesIntegral()
        if (mod(step,100) == 0) then
          write(eid,*) step*dt, S%getKineticEnergy(), S%getPotentialEnergy()
        endif

        ! Condition to stop preload
        do n = 1,S%number_of_mass_points
            if (S%mass_points(n)%X(3) <= -0.338_dp) exit preload
        end do

    end do preload

    ! Unload the solid
    do n = 1,S%number_of_mass_points
        S%mass_points(n)%Fs(3) = 0.0_dp
        S%mass_points(n)%v = ZEROV
        S%mass_points(n)%a = ZEROV
    end do

    ! Release the solid
    step = 0
    if (mod(step,1000) == 0) then! Print STL
        write(char_step,'(I0.8)') step
        outfile = 'data/S_'//char_step//'.stl'
        call S%writeSTL(outfile)
        outfile = 'data/S_'//char_step//'.gnu'
        call S%writeGNU(outfile)
    end if
    
    do while(step*dt <= 3.0_dp)
        
        call S%advance(dt, [0.0_dp, 0.0_dp, 0.0_dp])
        
        if (mod(step,5000) == 0) then! Print STL
            write(char_step,'(I0.8)') step
            outfile = 'data/S_'//char_step//'.stl'
            call S%writeSTL(outfile)
            outfile = 'data/S_'//char_step//'.gnu'
            call S%writeGNU(outfile)
        end if
    
        if (mod(step,100) == 0) then
            write(tid,'(E16.8,1x)',advance='no') (step-sstep)*dt
            do n = 1,nbcr
                write(tid,'(E16.8,1x)',advance='no') S%mass_points(ibcr(n))%x(3)
            end do
            write(tid,*) ''

            ! Output energy            
            write(eid,*) step*dt, S%getKineticEnergy(), S%getPotentialEnergy()
        endif
        print *, 'LOAD FREE: ', step*dt, S%internalForcesIntegral()
        step  = step + 1
    end do

    ! Free memory
    call S%destroy()

contains

    !===============================================================================================
    subroutine test_constraints(self)

        use lagrangian_solid_mod

        class(lagrangian_solid), intent(inout), target :: self

        ! At x = 0 set zero displacment, velocity and acceleration
        do n = 1,nbcl1
            self%mass_points(ibcl1(n))%X(1) = 0.0_dp
            self%mass_points(ibcl1(n))%X(3) = 0.0_dp
        end do

        ! At x = dl set zero vertical displacement, velocity and acceleration
        do n = 1,nbcl2
            self%mass_points(ibcl2(n))%X(3) = 0.0_dp
        end do

    end subroutine
    !===============================================================================================

end program main
