module lagrangian_solid_mod

    ! This module contains the basic interface for simulating a solid object discretized on a 
    ! Lagrangian mesh. The description is based on mass points connected by edges 
    ! (and triangles in 3D).

    use precision_mod
    use global_mod
    use marker_mod
    use edge_mod
    use solid_mod

    implicit none
    private
    public :: lagrangian_solid, lagrangian_solid_pointer, forcing_element, velocity_verlet, expImp

    type forcing_element
        !< The forcing Element is a general interface to retrive information from elements
        !< where the MLS forcing is applied. 
        type(marker), pointer :: C    !< marker point
        real(dp)    , pointer :: n(:) !< normal vector to forcing element
        real(dp)    , pointer :: A    !< reference area
    end type forcing_element

    type, extends(solid), abstract :: lagrangian_solid
        !< Asbtract Lagrangian solid
        integer                            :: number_of_mass_points
        type(marker), allocatable          :: mass_points(:)
        integer                            :: number_of_edges
        type(edge), allocatable            :: edges(:)
        real(dp)                           :: traslation(tdof)     ! center of mass traslation
        real(dp)                           :: rotation(rdof)       ! center of mass rotation
        real(dp), pointer                  :: rotation_center(:)   ! Rotation center coordinates
        real(dp)                           :: gamma = 0.0_dp
        logical                            :: is_open = .false.
        logical                            :: is_deformable = .false.
        integer                            :: nsubsteps = 1
        type(forcing_element), allocatable :: forcing_elements(:)

        procedure(self_subroutine_interafce), pointer :: apply_constraints => Null()
        procedure(deformation_solver)       , pointer :: advance_deformation => velocity_verlet
        contains
            ! Object independent procedures
            procedure, pass(self) :: advance
            procedure, pass(self) :: rigid_body_motion
            procedure, pass(self) :: update_center_of_mass
            procedure, pass(self) :: integrate_hydrodynamic_forces

            ! Object dependent procedures
            procedure(self_subroutine_interafce), deferred :: update
            procedure(self_subroutine_interafce), deferred :: compute_internal_forces
            procedure(self_subroutine_interafce), deferred :: interpolate_from_forcing_element_to_mass_point
            procedure(self_subroutine_interafce), deferred :: update_forcing_elements
    end type lagrangian_solid

    interface
        subroutine self_subroutine_interafce(self)
            import lagrangian_solid
            class(lagrangian_solid), intent(inout), target :: self
        end subroutine self_subroutine_interafce
        subroutine deformation_solver(self, dt)
            use precision_mod, only : dp
            import lagrangian_solid
            class(lagrangian_solid), intent(inout), target :: self
            real(dp)               , intent(in   )         :: dt
        end subroutine deformation_solver
    end interface
    type lagrangian_solid_pointer
        class(lagrangian_solid), pointer :: pS => Null()
    end type lagrangian_solid_pointer

contains
    
    !===============================================================================================
    subroutine advance(self, dt, g, fluid_density)

        ! Objective: perform one timestep of the structural solver due to applied forces.

        ! In/Out varialbes
        class(lagrangian_solid), intent(inout)           :: self
        real(dp)               , intent(in   )           :: dt 
        real(dp)               , intent(in   )           :: g(:)
        real(dp)               , intent(in   ), optional :: fluid_density

        ! Local variables
        integer  :: n, substep
        real(dp) :: bf(dofs), Vnp1(dofs), Xnp1(dofs), forces(dofs)

        ! If the solid is deformable
        if (self%is_deformable) then

            ! Setup the external forces vector as sum of gravity forces, hydrodynamic forces
            ! and additional external sources
            do n = 1,self%number_of_mass_points
                self%mass_points(n)%Fe(1:tdof) = self%mass_points(n)%Fh(1:tdof)  + &
                                                 self%mass_points(n)%m*g(1:tdof) + &
                                                 self%mass_points(n)%Fs(1:tdof)
            end do

            if (present(fluid_density)) then
                if (self%is_open) then 
                    ! Must account directly for buoyancy since the volume is zero
                    do n = 1,self%number_of_mass_points
                        self%mass_points(n)%Fe(1:tdof) = self%mass_points(n)%Fe(1:tdof) - &
                                        fluid_density*self%Vol*g(1:tdof)/self%number_of_mass_points
                    end do
                endif
            endif

            ! Solve deformation due to internal forces nsubstesp per timestep
            do substep = 1,self%nsubsteps
                call self%advance_deformation(dt/real(self%nsubsteps, dp))
            end do

        else

            ! For rigid solid body solve Newton's second law for the center of mass

            ! First integrate the local forces to get center of mass forces
            call self%integrate_hydrodynamic_forces()

            ! Add body forces
            bf = 0.0_dp
            if (self%is_open) then
                ! For open solid buoyancy must be added explicitly
                bf(1:Ndim) = self%Vol*(self%density - fluid_density)*g(1:Ndim)
            else
                bf(1:Ndim) = self%IM(1:Ndim)*g(1:Ndim)
            endif

            ! Assemble all forces
            forces = self%center_of_mass%Fh + self%center_of_mass%Fs + bf

            ! Find new values
            Vnp1 = self%center_of_mass%V + dt*forces/self%IM
            Xnp1 = self%center_of_mass%X + dt*0.5_dp*(Vnp1 + self%center_of_mass%V)

            ! Save traslation and rotation
            self%traslation = Xnp1(1:tdof) - self%center_of_mass%X(1:tdof)
            self%rotation = Xnp1(tdof+1:tdof+rdof) - self%center_of_mass%X(tdof+1:tdof+rdof)

            ! Update solution
            self%center_of_mass%V = Vnp1
            self%center_of_mass%X = Xnp1
            self%center_of_mass%A = forces/self%IM

            ! Moves all lagrangian points accordingly
            call self%rigid_body_motion()
        endif

    end subroutine advance
    !===============================================================================================

    !===============================================================================================
    subroutine integrate_hydrodynamic_forces(self)

        ! Objective: compute net hydrodynamic force on center of mass as sum over all 
        ! forcing element centers.

        ! In/Out variables
        class(lagrangian_solid), intent(inout), target :: self

        ! Local variables
        integer :: n

        ! Hydrodynamic forces
        self%center_of_mass%Fh = 0.0_dp
        self%center_of_mass%Fv = 0.0_dp
        self%center_of_mass%Fp = 0.0_dp

        ! Cycle over lagrangian markers
        do n = 1,size(self%forcing_elements)
            self%center_of_mass%Fh = self%center_of_mass%Fh + self%forcing_elements(n)%C%Fh
            self%center_of_mass%Fv = self%center_of_mass%Fv + self%forcing_elements(n)%C%Fv
            self%center_of_mass%Fp = self%center_of_mass%Fp + self%forcing_elements(n)%C%Fp
        end do

    end subroutine
    !===============================================================================================

    !===============================================================================================
    subroutine update_center_of_mass(self)

        class(lagrangian_solid), intent(inout) :: self

        ! Local variables
        integer :: n

        do n = 1,tdof
            self%center_of_mass%X(n) = 0.0_dp
        end do
        
        do n = 1,self%number_of_mass_points
            self%center_of_mass%X(1:tdof) = self%center_of_mass%X(1:tdof) + self%mass_points(n)%X
        end do
        self%center_of_mass%X(1:tdof) = self%center_of_mass%X(1:tdof)/self%number_of_mass_points

    end subroutine
    !===============================================================================================

    !===============================================================================================
    subroutine rigid_body_motion(self)

        ! Objective: update mass points and forcing points according to center of mass traslation 
        !            and rotation as consequence of rigid body motion.

        use euclidean_mod, only : Rx, Ry, Rz, crossProduct

        ! In/Out variables
        class(lagrangian_solid), intent(inout), target :: self

        ! Local variables
        integer  :: n
        real(dp) :: R(Ndim, Ndim), d(Ndim), V(tdof), A(tdof), omega(rdof), alpha(rdof)
        real(dp) :: omega_x_r(3), omega_x_omega_x_r(3), alpha_x_r(3)

        ! **** First step: update mass points **** 

        ! Build rotation matrix
#if DIM==3
        R = matmul(matmul(Rz(self%rotation(3)), &
                          Ry(self%rotation(2))), &
                          Rx(self%rotation(1)))
#else
        ! Build rotation matrix
        R(1,1) =  cos(self%rotation(1))
        R(1,2) = -sin(self%rotation(1))
        R(2,1) =  sin(self%rotation(1))
        R(2,2) =  cos(self%rotation(1))
#endif

        ! Apply traslation and rotation to each mass point
        do n = 1,self%number_of_mass_points
            ! Apply traslation:
            self%mass_points(n)%X = self%mass_points(n)%X + self%traslation
            
            ! Distance from rotation center
            d = self%mass_points(n)%X(1:Ndim) - self%rotation_center(1:Ndim)
            
            ! Apply rotation:
            self%mass_points(n)%X = self%rotation_center + matmul(R,d)
        end do
       
        ! **** Second step: update forcing points ****
        call self%update()
        do n = 1,size(self%forcing_elements)

            ! Distance between forcing point and rotation center
            d = self%forcing_elements(n)%C%X(1:Ndim) - self%rotation_center(1:Ndim)

            ! Select center of mass velocity and acceleracion
            V = self%center_of_mass%V(1:tdof)
            omega = self%center_of_mass%V(tdof+1:dofs)
            A = self%center_of_mass%A(1:tdof)
            alpha = self%center_of_mass%A(tdof+1:dofs)

            ! Velocity of the forcing point n is given by
            ! v = v_cm + omega x r
#if DIM==3
            omega_x_r = crossProduct(omega, d)
#else
            omega_x_r = crossProduct([0.0_dp, 0.0_dp, omega], [d(1), d(2), 0.0_dp])
#endif
            self%forcing_elements(n)%C%V(1:tdof) = V(1:tdof) + omega_x_r(1:tdof)
            !self%forcing_elements(n)%C%V(1) = self%center_of_mass%V(1) - self%center_of_mass%V(3)*d(2)
            !self%forcing_elements(n)%C%V(2) = self%center_of_mass%V(2) + self%center_of_mass%V(3)*d(1)
            
            ! The angular velocity is the same of the center of mass
            self%forcing_elements(n)%C%V(tdof+1:dofs) = omega!self%center_of_mass%V(tdof+1:dofs)

            ! Acceleration of the lagrangian marker l
            ! a = a_cm + omega x (omega x r) + alpha x r
#if DIM==3
            omega_x_omega_x_r = crossProduct(omega, crossProduct(omega, d))
            alpha_x_r = crossProduct(alpha, d)
#else
            omega_x_omega_x_r = crossProduct([0.0_dp, 0.0_dp, omega], &
                                    crossProduct([0.0_dp, 0.0_dp, omega], [d(1), d(2), 0.0_dp]))
            alpha_x_r = crossProduct([0.0_dp, 0.0_dp, alpha], [d(1), d(2), 0.0_dp])
#endif
            self%forcing_elements(n)%C%A(1:tdof) = A(1:tdof) + omega_x_omega_x_r(1:tdof) + &
                                                        alpha_x_r(1:tdof)
            ! self%forcing_elements(n)%C%A(1) = self%center_of_mass%A(1) - d(1)*self%center_of_mass%V(3)**2 - &
            !                                     self%center_of_mass%A(3)*d(2)
            ! self%forcing_elements(n)%C%A(2) = self%center_of_mass%A(2) - d(2)*self%center_of_mass%V(3)**2 + &
            !                                     self%center_of_mass%A(3)*d(1)

            ! The angular acceleration is the same
            !self%forcing_elements(n)%C%A(3) = self%center_of_mass%A(3)
            self%forcing_elements(n)%C%A(tdof+1:dofs) = alpha
        end do

    end subroutine
    !===============================================================================================

    !===============================================================================================
    subroutine velocity_verlet(self, dt)

        ! In/Out variables
        class(lagrangian_solid), intent(inout), target :: self
        real(dp)               , intent(in   )         :: dt

        ! Local variables
        integer               :: i
        type(marker), pointer :: mp

        !***** First step of the velocity verlet

        ! Evaluate internal forces from potential
        call self%compute_internal_forces()
       
        do i = 1,self%number_of_mass_points
            ! Select current mass point
            mp => self%mass_points(i) 
            ! Acceleration at n
            mp%a = (mp%Fe + mp%Fi)/mp%m
            ! Velocity at n + 1/2
            mp%v = mp%v + 0.5_dp*dt*mp%a
            ! Position at n + 1
            mp%X = mp%X + dt*mp%v 
        end do

        ! Enforce BC
        if (associated(self%apply_constraints)) then
           call self%apply_constraints()
        endif

        ! Update the solid with the new mass points position
        call self%update()

        !**** Second step of the velocity verlet

        ! Get new internal forces
        call self%compute_internal_forces()

        do i = 1,self%number_of_mass_points
            ! Select current mass point
            mp => self%mass_points(i)
            ! Accleration at n + 1
            mp%a = (mp%Fe + mp%Fi)/mp%m
            ! Veloicyt at n + 1
            mp%v = mp%v*(1.0_dp - self%gamma) + 0.5_dp*dt*mp%a
        end do

        ! Enforce BC
        if (associated(self%apply_constraints)) then
            call self%apply_constraints()
        endif
 
        ! Update the solid with the new mass points position
        call self%update()

    end subroutine
    !===============================================================================================

    !===============================================================================================
    subroutine expImp(self, dt)

        ! In/Out variables
        class(lagrangian_solid), intent(inout), target :: self
        real(dp)               , intent(in   )         :: dt

        ! Local variables
        integer               :: i
        real(dp)              :: a, b, c
        type(marker), pointer :: mp

        ! **** First step is explicit
        ! Evaluate internal forces from potential at step n
        call self%compute_internal_forces()
       
        do i = 1,self%number_of_mass_points
            ! Select current mass point
            mp => self%mass_points(i) 

            ! Coefficients
            a = 0.5_dp*dt/mp%m
            b = a*self%gamma
            c = 1.0_dp/(1.0_dp + b)
    
            ! Velocity at n + 1/2
            mp%v = mp%v*(1.0_dp - b) + a*(mp%Fe + mp%Fi)
            
            ! Position at n + 1
            mp%X = mp%X + dt*mp%v
        end do

        ! Enforce BC on position at n + 1
        if (associated(self%apply_constraints)) then
           call self%apply_constraints()
        endif

        ! Update the solid with the new mass points position
        call self%update()

        ! **** Second step is implicit

        ! Evaluate internal forces from potential at step n+1
        call self%compute_internal_forces()

        do i = 1,self%number_of_mass_points
            ! Select current mass point
            mp => self%mass_points(i)
            
            ! Coefficients
            a = 0.5_dp*dt/mp%m
            b = a*self%gamma
            c = 1.0_dp/(1.0_dp + b)
    
            ! Veloicyt at n + 1
            mp%v = (mp%v + a*(mp%Fe + mp%Fi))*c
        end do

        ! Enforce BC on position at n + 1
        if (associated(self%apply_constraints)) then
            call self%apply_constraints()
        endif
 
        ! Update the solid with the new mass points position
        call self%update()

    end subroutine
    !===============================================================================================

end module lagrangian_solid_mod
