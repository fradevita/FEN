program filament

    ! Test for the structural spring-mass solver

    use precision_mod        , only : dp
    use lagrangian_solid_mod

    implicit none

    ! Parameters
    real(dp), parameter :: length = 1.0_dp  ! Lenght of the filament
    real(dp), parameter :: Fr = 10.0_dp     ! Froud Number, Fr = gL/U^2

    ! Variables
    integer :: step, n, out_id, substep
    real(dp) :: time, dt
    type(lagrangian_solid) :: test_solid

    ! Set the mass of the solid body: rho*Volume
    test_solid%M = 1.00_dp

    ! Turn on the deformable flag
    test_solid%is_deformable = .true.

    ! Turn on the open flag
    test_solid%is_open = .true.

    ! Create the lagrangian solid from the mesh file
    call test_solid%create('mesh.txt')

    ! Half the last mass point mass
    ! test_solid%mass_points(test_solid%number_of_mass_points)%M = &
    !     test_solid%mass_points(test_solid%number_of_mass_points)%M*0.5_dp

    ! Set the elastic in-plane constant of the springs
    test_solid%ke = 1.0e+6_dp
    test_solid%edges%ke = 1.0e+6_dp

    ! Set the bending constant
    test_solid%kb = 0.0_dp

    ! Set the constraints procedure
    test_solid%apply_constraints => test_constraints

    ! Set timestep
    step = 0
    time = 0.0_dp
    dt = 2.0e-4_dp

    ! Add gravity
    do n = 1,test_solid%number_of_mass_points
        test_solid%mass_points(n)%Fe(2) = -Fr*test_solid%mass_points(n)%m
    end do

    ! Open output file
    open(newunit = out_id, file = 'out.txt')

    step = 0
    time = 0.0_dp

    call test_solid%print_configuration(step)

    !==== Start Time loop ===================================================================
    time_loop: do while(time < 20.0_dp)

        step = step + 1
        time = time + dt
        write(*,*) 'setp: ', step, 'time: ', time, 'dt: ', dt

        ! Solve deformation due to internal foces nsubstep times per timestep
        do substep = 1,test_solid%nsubsteps
            ! Velocity Verlet
            call test_solid%velocity_verlet(dt/real(test_solid%nsubsteps, dp))
        end do

        if (mod(step,100) == 0) then
            call test_solid%print_configuration(step)
        endif

        write(out_id,*) time, test_solid%mass_points(test_solid%number_of_mass_points)%X
        flush(out_id)

    end do time_loop

contains

    !========================================================================================
    subroutine test_constraints(self)

        class(lagrangian_solid), intent(inout) :: self

        self%mass_points(1)%X = 0.0_dp
        self%mass_points(1)%V = 0.0_dp

    end subroutine test_constraints
    !========================================================================================

end program filament
