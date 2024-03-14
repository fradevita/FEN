program filament

    ! Test for the structural spring-mass solver

    use precision_mod           , only : dp
    use lagrangian_solid_1D_mod

    implicit none

    ! Parameters
    real(dp), parameter :: length = 1.0_dp  ! Lenght of the filament
    real(dp), parameter :: Fr = 10.0_dp     ! Froud Number, Fr = gL/U^2

    ! Variables
    integer                   :: step, n, out_id
    real(dp)                  :: time, dt
    type(lagrangian_solid_1D) :: S

    ! Set the mass of the solid body: rho*Volume
    S%mass = 1.00_dp

    ! Turn on the deformable flag
    S%is_deformable = .true.

    ! Turn on the open flag
    S%is_open = .true.

    ! Create the lagrangian solid from the mesh file
    call S%create('mesh.txt')

    ! Half the last mass point mass
    ! S%mass_points(S%number_of_mass_points)%M = &
    !     S%mass_points(S%number_of_mass_points)%M*0.5_dp

    ! Set the elastic in-plane constant of the springs
    ! to an arbitrary value to enforce inextensibility
    S%edges%ks = 1.0e+6_dp

    ! Set the bending constant
    S%kb = 0.0_dp

    ! Set the constrain procedure
    S%apply_constraints => test_constraints

    ! Set timestep
    step = 0
    time = 0.0_dp
    dt = 5.0e-6_dp

    ! Open output file
    open(newunit = out_id, file = 'out.txt')

    step = 0
    time = 0.0_dp

    call S%print_configuration(step)

    !==== Start Time loop ===================================================================
    time_loop: do while(time < 20.0_dp)

        step = step + 1
        time = time + dt
        write(*,*) 'setp: ', step, 'time: ', time, 'dt: ', dt

        ! Solve deformation due to internal foces nsubstep times per timestep
        call S%advance(dt, [0.0_dp, -Fr])

        if (mod(step,200) == 0) then
            call S%print_configuration(step)
            write(out_id,*) time, S%mass_points(S%number_of_mass_points)%X
            flush(out_id)
        endif

    end do time_loop

contains

    !========================================================================================
    subroutine test_constraints(self)

        use lagrangian_solid_mod

        class(lagrangian_solid), intent(inout), target :: self

        self%mass_points(1)%X = 0.0_dp
        self%mass_points(1)%V = 0.0_dp

    end subroutine test_constraints
    !========================================================================================

end program filament
