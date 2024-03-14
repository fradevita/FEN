program Euler_Bernoulli

    ! Test for the structural spring-mass solver. The test compare the deflaction of a
    ! filament with the solution of the Euler-Bernoulli beam theory.

    use precision_mod           , only : dp
    use global_mod              , only : pi
    use marker_mod
    use lagrangian_solid_1D_mod , only : lagrangian_solid_1D

    implicit none

    ! Parameters
    real(dp), parameter :: length = 1.0_dp                  ! Lenght of the beam [m]
    real(dp), parameter :: width =  0.01_dp                 ! Width of the beam [m]
    real(dp), parameter :: thickness = 1.0e-3_dp            ! Thickness of the beam [m]
    real(dp), parameter :: E = 210.1e+10_dp                 ! Young's Modulus [Pa]
    real(dp), parameter :: IM = width*thickness**3/12.0d0   ! Second moment of area [m^4]
    real(dp), parameter :: D = E*IM                         ! Bending stiffness
    real(dp), parameter :: rhos = 7850.0_dp                 ! Density of the beam [kg/m^3]
    real(dp), parameter :: F = 1.0e-4_dp                    ! Point load on the free edge [N]

    ! Variables
    integer                       :: step, n, r, res
    real(dp)                      :: time, dt
    real(dp), allocatable         :: diff(:)
    type(lagrangian_solid_1D)     :: S
    type(marker), allocatable     :: oldX(:)
    character(len=3)              :: sres
    character(len=11)             :: outfile
    character(len=12)             :: meshfile
  
    resolution_cycle: do r = 1,4

        ! Current resolution
        res = 2**(4+r)

        ! Set the mass of the solid body: rho*Volume
        S%vol = length*width*thickness
        S%density = rhos
        S%mass = S%vol*S%density

        ! Turn on the deformable flag
        S%is_deformable = .true.

        ! Turn on the open flag
        S%is_open = .true.

        ! Create the lagrangian solid from the mesh file
        write(sres,'(I0.3)') res
        meshfile = 'mesh_'//sres//'.txt'
        print *, meshfile
        call S%create(meshfile)

        ! Half the last mass point mass
        !S%mass_points(S%number_of_mass_points)%M = &
        !     S%mass_points(S%number_of_mass_points)%M*0.5_dp

        ! Set the elastic in-plane constant of the springs
        ! Eq. 25 of de Tullio and Pascazio JCP 2016.
        ! ke = E*h*SumAi/l**2
        ! l is equal to length / number of mass points
        ! A = l*width
        S%edges%ks = E*thickness*width*real(S%number_of_edges, dp)/length

        ! Set the bending constant
        S%kb = D*real(S%number_of_edges, dp)
    
        ! Set the constraints procedure
        S%apply_constraints => test_constraints

        ! Set timestep
        step = 0
        time = 0.0_dp
        dt = 1.0e-4_dp
        if (r == 3) dt = dt / 2.0_dp
        if (r == 4) dt = dt / 4.0_dp
        
        ! Apply the load on the trailing edge
        S%mass_points(S%number_of_mass_points)%Fs(2) = -F
        
        ! Set a viscosity to dump oscillations at grid scale.
        S%c = 5.0e-6_dp

        ! Allocate variables for steady state check
        allocate(diff(S%number_of_mass_points))
        allocate(oldX(S%number_of_mass_points))
        do n = 1,S%number_of_mass_points
            oldX(n)%X = S%mass_points(n)%X
        end do
        diff = 1.0_dp
        
        !==== Start Time loop ===================================================================
        time_loop: do
            
            step = step + 1
            time = time + dt
            write(*,*) 'setp: ', step, 'time: ', time, 'dt: ', dt, 'dX: ', maxval(diff)
            
            ! Solve deformation due to internal forces nsubstep times per timestep
           call S%advance(dt, [0.0_dp, 0.0_dp, 0.0_dp])
            
            ! Check for steady state
            do n = 1,S%number_of_mass_points
                diff(n) = maxval(abs(oldX(n)%X - S%mass_points(n)%X))
            end do
            if (maxval(diff) < 1.0e-12_dp .and. step > 1) exit time_loop

            ! Update old position vector
            do n = 1,S%number_of_mass_points
                oldX(n)%X = S%mass_points(n)%X
            end do
            
        end do time_loop

        ! Print finale configuration
        outfile = 'out_'//sres//'.txt'
        open(unit = 1, file = outfile)
        do n = 1,S%number_of_mass_points
            write(1,*) S%mass_points(n)%X(1:2)
        end do
        close(1)

        call S%destroy()
        deallocate(diff)
        deallocate(oldX)
  end do resolution_cycle
  
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
  
end program Euler_Bernoulli
