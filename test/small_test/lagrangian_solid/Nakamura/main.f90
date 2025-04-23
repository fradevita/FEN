program main

    use mpi
    use precision_mod 
    use global_mod              , only : ierror, myrank, pi
    use lagrangian_solid_2D_mod
    use lagrangian_solid_mod
    use IO_mod                  , only : stdout, print_error_message
    use euclidean_mod           , only : ZEROV

    implicit none

    ! Parameters
    real(dp), parameter :: a   = 3.3e-6_dp  ! sphere radius       [m]
    real(dp), parameter :: rho = 1050.0_dp  ! density [kg/m^3]
    
    real(dp), parameter :: Ks  = 15.0e-6_dp
    real(dp), parameter :: Kb  = 0.25*1.0e-10_dp
    real(dp), parameter :: Ka  = 3.0e-3_dp
    real(dp), parameter :: KTA = 2.0e-3_dp
    real(dp), parameter :: Kv  = 50.0_dp

    ! Variables
    integer                           :: i, step, out_id
    real(dp)                          :: dt, time, Ws, Wb, Wa, WTA, WV, W, Wo, Wdiff
    type(lagrangian_solid_2D), target :: S
    character(len= 7)                 :: sstep

    ! **** Setup MPI *****************************************************
    call mpi_init(ierror)
    call mpi_comm_rank(mpi_comm_world, myrank, ierror)
    !*********************************************************************

    !**** Define the solid sphere ******************************************************************
    S%is_deformable = .true.
    call S%loadFromFile('Sphere_3.txt')
    S%density = rho
    S%Vol = S%V0
    S%mass = S%density*S%Vol
    allocate(S%IM(6))
    S%IM(1:6) = S%mass

    ! Setup material properties
    call S%setMass(S%mass)
    call S%setStretchingConstant(Ks)
    call S%setBendingConstant(Kb)
    call S%setAreaConstant(Ka)
    call S%setTotalAreaConstant(Kb)
    call S%setVolumeConstant(Kv)

    ! Set resting theta
    do i = 1,S%number_of_edges
        S%edges(i)%theta0 = 0.0_dp
    end do

    ! Critically damped system
    S%gamma = sqrt(S%mass_points(1)%m*Ks)

    ! Volume penalization
    S%V0 = S%V0*0.6_dp

    ! Select time marching scheme
    S%advance_deformation => expImp
    !***********************************************************************************************

    !**** Start Time loop **************************************************************************
    time = 0.0_dp
    dt = 1.0e-7_dp
    open(newunit = out_id, file = 'energy.csv')
    write(out_id, '(A)') 't,Ws,Wb,Wa,WTA,WV,W'
    Wo = S%getPotentialEnergy()

    time_loop: do while(time < 1.2_dp)

        step = step + 1
        time = time + dt
        
        call S%advance(dt, [0.0_dp, 0.0_dp, 0.0_dp])
        W = S%getPotentialEnergy()
        Wdiff = abs(Wo - W)/abs(Wo + 1.0e-30_dp)
        Wo = W
        print *, 'step :', step, 'time: ', time, 'Wdiff: ', Wdiff 

        ! Evaluate energies
        Ws  = S%getStretchingEnergy()
        Wb  = S%getBendingEnergy()
        Wa  = S%getAreaEnergy()
        WTA = S%getTotalAreaEnergy()
        WV  = S%getVolumeEnergy()
        write(out_id,'(*(E16.8,:,","))') time, Ws, Wb, Wa, WTA, WV, S%getPotentialEnergy()

        if (mod(step,100) == 0) then
            write(sstep,'(I0.7)') step
            call S%writeSTL('data/S_'//sstep//'.stl')
        endif

        if (Wdiff < 1.0e-13_dp .and. step > 10000) exit time_loop

    end do time_loop
    !***********************************************************************************************

    ! free memory
    call S%destroy()
    call mpi_finalize(ierror)

end program main
