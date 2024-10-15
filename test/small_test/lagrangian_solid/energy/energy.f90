program energy

    ! Test for the structural spring-mass solver

    use precision_mod           , only : dp
    use global_mod              , only : pi
    use lagrangian_solid_1D_mod , only : lagrangian_solid_1D

    implicit none

    ! Parameters
    real(dp), parameter :: length = 1.0_dp                ! Lenght of the plate [m]
    real(dp), parameter :: width =  0.1_dp                ! Width of the plate [m]
    real(dp), parameter :: thickness = 1.0e-3_dp          ! Thickness of the plate [m]
    real(dp), parameter :: E = 210.1e+10_dp               ! Young's Modulus [Pa]
    real(dp), parameter :: IM = width*thickness**3/12.0d0 ! Second moment of area [m^4] 
    real(dp), parameter :: rhos = 7850.0_dp               ! Density of the plate [kg/m^3]
    real(dp), parameter :: F = 22.5_dp                    ! Point load on the free edge [N]
    real(dp), parameter :: nu = 0.3_dp                    ! Poisson modulus
    real(dp), parameter :: ymax = -0.338_dp               ! Maximum displacment in preload

    ! Variables
    integer                   :: step, eid
    real(dp)                  :: time, dt, Work, Ep, Ek
    type(lagrangian_solid_1D) :: S

    ! Set the mass of the solid body: rho*Volume
    S%vol = length*width*thickness
    S%density = rhos
    S%mass = S%Vol*S%density

    ! Turn on the deformable flag
    S%is_deformable = .true.

    ! Turn on the open flag
    S%is_open = .true.

    ! Create the lagrangian solid from the mesh file
    call S%create('mesh.txt')

    ! Half the last mass point mass
    S%mass_points(S%number_of_mass_points)%m = S%mass_points(S%number_of_mass_points)%m*0.5_dp

    ! Set the elastic in-plane constant of the springs
    ! Eq. 25 of de Tullio and Pascazio JCP 2016.
    ! ke = E*h*SumAi/l**2
    ! l is equal to length / number of mass points
    ! A = l*width
    S%edges%ks = E*thickness*width*real(S%number_of_edges, dp)/length

    ! Set the bending constant
    S%kb = E*IM*real(S%number_of_edges, dp)
    
    ! Set the constraints procedure
    S%apply_constraints => test_constraints
    
    ! Set timestep
    step = 0
    time = 0.0_dp
    dt = 5.0e-7_dp

    ! In the first run apply the load on the free tip
    S%mass_points(S%number_of_mass_points)%Fs(2) = -F

    open(eid, file = 'energy.dat')

    !**** Start Time loop
    time_loop: do while(time < 1.0_dp)

        step = step + 1
        time = time + dt
        write(*,*) 'setp: ', step, 'time: ', time, 'dt: ', dt
        
        ! Solve deformation due to internal foces
        call S%advance(dt, [0.0_dp, 0.0_dp, 0.0_dp])

        if (mod(step,100) == 0) then
            ! Evaluate work and mechanical energy
            Work = F*abs(S%mass_points(S%number_of_mass_points)%X(2))
            Ep = S%getPotentialEnergy()
            Ek = S%get_kinetic_energy()
            write(eid,*) time, Work, Ep, Ek
            flush(eid)
        endif

    end do time_loop

contains

    !===============================================================================================
    subroutine test_constraints(self)
        
        use lagrangian_solid_mod

        class(lagrangian_solid), intent(inout), target :: self

        self%mass_points(1)%X = 0.0_dp
        self%mass_points(1)%V = 0.0_dp
        self%mass_points(2)%X(2) = 0.0_dp
        self%mass_points(2)%V(2) = 0.0_dp
        self%mass_points(2)%A(2) = 0.0_dp
            
    end subroutine test_constraints
    !===============================================================================================
    
end program
