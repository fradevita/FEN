program cantilever

    ! Test for the structural spring-mass solver

    use precision_mod        , only : dp
    use global_mod           , only : pi
    use lagrangian_solid_mod , only : lagrangian_solid

    implicit none

    ! Parameters
    real(dp), parameter :: length = 1.0_dp                ! Lenght of the plate [m]
    real(dp), parameter :: width =  0.1_dp                ! Width of the plate [m]
    real(dp), parameter :: thickness = 1.0e-3_dp          ! Thickness of the plate [m]
    real(dp), parameter :: E = 210.1e+10_dp               ! Young's Modulus [Pa]
    real(dp), parameter :: IM = width*thickness**3/12.0d0 ! Second moment of area [m^4] 
    real(dp), parameter :: rhos = 7850                    ! Density of the plate [kg/m^3]
    real(dp), parameter :: F = 22.5_dp                    ! Point load on the free edge [N]
    real(dp), parameter :: nu = 0.3_dp                    ! Poisson modulus
    real(dp), parameter :: ymax = -0.338_dp               ! Maximum displacment in preload

    ! Variables
    integer                :: step, n, out_id, substep
    real(dp)               :: time, dt
    type(lagrangian_solid) :: test_solid

    ! Set the mass of the solid body: rho*Volume
    test_solid%M = rhos*length*width*thickness

    ! Turn on the deformable flag
    test_solid%is_deformable = .true.

    ! Turn on the open flag
    test_solid%is_open = .true.

    ! Create the lagrangian solid from the mesh file
    call test_solid%create('mesh.txt')

    ! Half the last mass point mass
    test_solid%mass_points(test_solid%number_of_mass_points)%M = &
        test_solid%mass_points(test_solid%number_of_mass_points)%M*0.5_dp

    ! Set the elastic in-plane constant of the springs
    ! Eq. 25 of de Tullio and Pascazio JCP 2016.
    ! ke = E*h*SumAi/l**2
    ! l is equal to length / number of mass points
    ! A = l*width
    test_solid%edges%ke = E*thickness*width*real(test_solid%number_of_edges, dp)/length
    test_solid%ke = E*thickness*width*real(test_solid%number_of_edges, dp)/length

    ! Set the bending constant
    test_solid%kb = E*IM*real(test_solid%number_of_edges, dp)
    
    ! Set the constraints procedure
    test_solid%apply_constraints => test_constraints
    
    ! Set timestep
    step = 0
    time = 0.0_dp
    dt = 1.0e-4_dp

    ! In the first run apply the load on the free tip
    test_solid%mass_points(test_solid%number_of_mass_points)%Fe(2) = -F

    !==== Start Time loop ===================================================================
    preload_time_loop: do while(time < 3.0_dp)

        step = step + 1
        time = time + dt
        write(*,*) 'setp: ', step, 'time: ', time, 'dt: ', dt
        
        ! Solve deformation due to internal foces nsubstep times per timestep
        do substep = 1,test_solid%nsubsteps
            ! Velocity Verlet
            call test_solid%velocity_verlet(dt/real(test_solid%nsubsteps, dp))
        end do

        ! Check for maximum deflection
        if (test_solid%mass_points(test_solid%number_of_mass_points)%X(2) <= ymax) exit preload_time_loop
    end do preload_time_loop
    
    ! Update the structure
    call test_solid%update_lagrangian_markers

    ! Open output file
    open(newunit = out_id, file = 'out.txt')

    ! Set external forces to zero
    test_solid%mass_points(test_solid%number_of_mass_points)%Fe(2) = 0.0_dp

    ! Set to zero velocity and accelerations
    do n = 1,test_solid%number_of_mass_points
        test_solid%mass_points(n)%V = 0.0_dp
        test_solid%mass_points(n)%A = 0.0_dp
    end do

    step = 0
    time = 0.0_dp
    !==== Start Time loop ===================================================================
    time_loop: do while(time < 3.0_dp)

        step = step + 1
        time = time + dt
        write(*,*) 'setp: ', step, 'time: ', time, 'dt: ', dt
        
        ! Solve deformation due to internal foces nsubstep times per timestep
        do substep = 1,test_solid%nsubsteps
            ! Velocity Verlet
            call test_solid%velocity_verlet(dt/real(test_solid%nsubsteps, dp))
        end do
        
        write(out_id,*) step*dt, test_solid%mass_points(test_solid%number_of_mass_points)%X(2)
        flush(out_id)

    end do time_loop

contains

    !========================================================================================
    subroutine test_constraints(self)

        class(lagrangian_solid), intent(inout) :: self

        self%mass_points(1)%X = 0.0_dp
        self%mass_points(1)%V = 0.0_dp
        self%mass_points(2)%X(2) = 0.0_dp
        self%mass_points(2)%V(2) = 0.0_dp
        self%mass_points(2)%A(2) = 0.0_dp
            
    end subroutine test_constraints
    !========================================================================================
  
end program cantilever
