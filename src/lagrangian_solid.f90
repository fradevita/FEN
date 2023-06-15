!> ! Define the lagrangian solid objects for the Immersed Boundary Method
module lagrangian_solid_mod

    use precision_mod
    use global_mod
    use marker_mod
    use solid_mod

    implicit none
    private
    public :: lagrangian_solid, mass_point, lagrangian_solid_pointer, edge

    ! Define the mass point type
    type, extends(marker) ::  mass_point
        !< A mass point is a marker with a mass connected to one or more edges.
        real(dp)                           :: Fi(tdof) = 0.0_dp !< Internal forces vector
        integer                            :: number_of_edges   !< self exp.
        integer, dimension(:), allocatable :: edges_index       !< self exp.
    end type mass_point

    ! Define the edge type
    type edge
        integer                   :: mass_point_index(2)
        type(mass_point), pointer :: x1 !< pointer to first mass point
        type(mass_point), pointer :: x2 !< pointer to second mass point
        real(dp)                  :: l0 !< stress-free length
        real(dp)                  :: l  !< actual length
        real(dp), dimension(2)    :: n  !< Normal vector to the edge        
        type(marker)              :: C  !< Centroid of the edge
        real(dp)                  :: ke !< elastic constant for in-plane deformation
    contains
        procedure :: update_length                !< self expl.
        procedure :: update_norm                  !< self expl.
        procedure :: update_centroid              !< self expl.
        procedure :: get_interpolated_mass_values !< Get values from adjacent mass points
    end type edge

    ! Define a solid object
    type, extends(solid) :: lagrangian_solid
        integer                       :: number_of_mass_points   !< self exp.
        integer                       :: number_of_edges         !< self exp.
        type(mass_point), allocatable :: mass_points(:)          !< array of mass points
        type(edge), allocatable       :: edges(:)                !< array of edges
        real(dp)                      :: Vol                     !< volume of the solid
        real(dp)                      :: rho                     !< density of the solid
        real(dp)                      :: M(tdof+rdof) = 1.0_dp   !< inertial mass (applied to the center of mass)
        real(dp)                      :: tra(tdof)               !< rigid body traslation
        real(dp)                      :: rot(rdof)               !< rigid body rotation
        integer                       :: nsubsteps = 200         !< number of substep for the structural solver
        real(dp)                      :: ke                      !< elastic constant in-plane deformation
        real(dp)                      :: kb                      !< elastic constant for bending 
        real(dp)                      :: c = 0.0_dp              !< damping for the structural solver
        integer                       :: output_file_id          !< id number for the output file
        logical                       :: is_deformable = .false. !< flag for solving deformation
        logical                       :: is_open = .false.       !< flag for open structure
        logical                       :: is_out = .false.        !< flag for checking if the solid is outside the domain
        procedure(constraints), pass(self), pointer :: apply_constraints => Null()
    contains
        procedure :: create
        procedure :: get_center_of_mass                    
        procedure :: update_center_of_mass
        procedure :: rigid_body_motion
        procedure :: compute_internal_forces
        procedure :: integrate_hydrodynamic_forces
        procedure :: update_lagrangian_markers
        procedure :: print_configuration
        procedure :: get_total_length
        procedure :: get_potential_energy
        procedure :: get_kinetic_energy
        procedure :: velocity_verlet
        procedure :: destroy
    end type lagrangian_solid

    abstract interface
        subroutine constraints(self)
            import lagrangian_solid
            class(lagrangian_solid), intent(inout) :: self
        end subroutine constraints
    end interface
    type lagrangian_solid_pointer
        class(lagrangian_solid), pointer :: pS => Null()
    end type lagrangian_solid_pointer

contains

    !========================================================================================
    subroutine create(self, filename, name)

        ! This subroutine initialize the variables of the solid.
        ! The file filename must contains the location of the mass points.

        ! In/Out variables
        class(lagrangian_solid)    , intent(inout), target   :: self
        character(len=*), intent(in   )           :: filename
        character(len=*), intent(in   ), optional :: name

        ! Local variables
        integer  :: fid, nl, io, n, l

        ! Open the mesh file
        open(newunit = fid, file = filename)

        ! Read the number of lines, which correspond to the number of mass points
        nl = 0
        do
            read(fid,*,iostat = io)
            if (io /= 0) exit
            nl = nl + 1
        end do
        rewind(fid)

        ! **** Setup mass points ****
        ! Set the number of mass points
        self%number_of_mass_points = nl

        ! Allocate the array of mass points
        allocate(self%mass_points(self%number_of_mass_points))

        ! Read the position of mass points
        do n = 1,self%number_of_mass_points
            
            ! Allocate the mass point variables
            call self%mass_points(n)%create(tdof) ! only traslational motion

            ! Read the position from file
            read(fid,*) self%mass_points(n)%X

            ! Set the discrete mass of the point
            self%mass_points(n)%m = self%M(1)/self%number_of_mass_points
        end do
        close(fid)

        ! **** Setup Center of Mass ****
        ! Compute center of mass of the solid body
        call self%center_of_mass%create(dofs)
        call self%update_center_of_mass()
        
        ! **** Setup edges ****
        ! Set the number of edges
        ! If the structure is closed the number of edges is equal to the number of mass points,
        ! if it is open the number of edges is equal to the number of mass points -1.
        if (self%is_open) then
            self%number_of_edges = self%number_of_mass_points - 1
        else
            self%number_of_edges = self%number_of_mass_points
        endif

        ! Allocate the array of edges
        allocate(self%edges(self%number_of_edges))

        ! Initialize all edges
        do l = 1,self%number_of_edges - 1
            ! Each edge connect two masses, select the indexes
            self%edges(l)%mass_point_index = [l, l + 1]
            ! Each edge connect two masses, select the indexes
            self%edges(l)%x1 => self%mass_points(l)
            self%edges(l)%x2 => self%mass_points(l+1)
        end do

        ! Fix the last edge
        if (self%is_open) then
            l = self%number_of_edges
            self%edges(l)%x1 => self%mass_points(l)
            self%edges(l)%x2 => self%mass_points(l+1)
            self%edges(l)%mass_point_index = [l, l + 1]
        else
            l = self%number_of_edges
            self%edges(l)%x1 => self%mass_points(l)
            self%edges(l)%mass_point_index(2) = 1
            self%edges(l)%x2 => self%mass_points(1)
        end if

        ! Setup edge variables
        do l = 1,self%number_of_edges
            ! Compute edge length
            call self%edges(l)%update_length()

            ! Set the initial equilibrium length
            self%edges(l)%l0 = self%edges(l)%l

            ! Compute the normal vector
            call self%edges(l)%update_norm()

            ! Allocate centroid variables
            call self%edges(l)%C%create(dofs)

            ! Compute the centroid
            call self%edges(l)%update_centroid()
        end do

        ! Set the index of the mass points. By default consider the structure closed
        do n = 1,self%number_of_mass_points
            self%mass_points(n)%number_of_edges = 2
            allocate(self%mass_points(n)%edges_index(2))
            self%mass_points(n)%edges_index = [n - 1, n]
        end do
        self%mass_points(1)%edges_index(1) = self%number_of_edges
        self%mass_points(self%number_of_mass_points)%edges_index(2) = 1

        ! If the solid is open, the first and last mass point have only one edge, fix it.
        if (self%is_open) then
            deallocate(self%mass_points(1)%edges_index)
            self%mass_points(1)%number_of_edges = 1
            allocate(self%mass_points(1)%edges_index(1))
            self%mass_points(1)%edges_index(1) = 1
            deallocate(self%mass_points(self%number_of_mass_points)%edges_index)
            self%mass_points(self%number_of_mass_points)%number_of_edges = 1
            allocate(self%mass_points(self%number_of_mass_points)%edges_index(1))
            self%mass_points(self%number_of_mass_points)%edges_index(1) = self%number_of_edges
        end if

        if (present(name)) self%name = name

        ! Open output file
        open(newunit = self%output_file_id, file = trim('data/'//self%name))

    end subroutine create
    !========================================================================================

    !========================================================================================
    subroutine update_center_of_mass(self)

        ! Find the center of mass location of the solid.

        ! In/Out variables
        class(lagrangian_solid), intent(inout) :: self

        ! Local variables
        integer  :: n
        real(dp) :: Xcm(tdof)

        Xcm = 0.0_dp
        do n = 1,self%number_of_mass_points
            Xcm = Xcm + self%mass_points(n)%X
        end do

        self%center_of_mass%X(1:tdof) = Xcm/real(self%number_of_mass_points, dp)

    end subroutine
    !========================================================================================

    !========================================================================================
    function get_center_of_mass(self) result(Xcm)

        ! Find the center of mass location of the solid.

        ! In/Out variables
        class(lagrangian_solid), intent(in) :: self

        ! Local variables
        integer  :: n
        real(dp) :: Xcm(tdof)

        Xcm = 0.0_dp
        do n = 1,self%number_of_mass_points
            Xcm = Xcm + self%mass_points(n)%X
        end do

        Xcm = Xcm/real(self%number_of_mass_points, dp)

    end function
    !========================================================================================

    !========================================================================================
    subroutine update_length(self)

        ! In/Out variables
        class(edge), intent(inout) :: self

        self%l = sqrt( (self%x1%X(1) - self%x2%X(1))**2 + (self%x1%X(2) - self%x2%X(2))**2 )

    end subroutine
    !========================================================================================

    !========================================================================================
    subroutine update_norm(self)

        ! In/Out variables
        class(edge)     , intent(inout) :: self
       
        ! Local variables
        real(dp) :: t(2), mod_t

#if DIM==3
#else
        ! Tangent vector
        t = self%x2%X - self%x1%X
        mod_t = sqrt(t(1)**2 + t(2)**2 + 1.0e-14_dp)
        t = t / mod_t

        ! Get norm vector
        self%n = [t(2), -t(1)]
#endif

    end subroutine
    !========================================================================================

    !========================================================================================
    subroutine update_centroid(self)

        ! In/Out variables
        class(edge), intent(inout) :: self
      
        ! Centroid location
        self%C%X = (self%x1%X*self%x1%m + self%x2%X*self%x2%m)/(self%x1%m + self%x2%m)

    end subroutine
    !========================================================================================

    !========================================================================================
    function get_total_length(self) result(L)

        ! In/Out variables
        class(lagrangian_solid), intent(in) :: self

        ! Local variables
        integer  :: n
        real(dp) :: L

        L = 0.0_dp
        do n = 1,self%number_of_edges
            L = L + self%edges(n)%l
        end do

    end function
    !========================================================================================

    !========================================================================================
    subroutine get_interpolated_mass_values(self, mass_points)

        ! In/Out variables
        class(edge)     , intent(inout) :: self
        type(mass_point), intent(in   ) :: mass_points(:)

        ! Local variables
        type(mass_point) :: mp1, mp2

        ! Select mass points
        mp1 = mass_points(self%mass_point_index(1))
        mp2 = mass_points(self%mass_point_index(2))

        ! Interpolate values on the lagrangian marker location
        self%C%V = (mp1%V*mp1%m + mp2%V*mp2%m)/(mp1%m + mp2%m)
        self%C%A = (mp1%A*mp1%m + mp2%A*mp2%m)/(mp1%m + mp2%m)

    end subroutine get_interpolated_mass_values
    !========================================================================================

    !========================================================================================
    subroutine compute_internal_forces(self)

        ! In/Out variables
        class(lagrangian_solid), intent(inout) :: self

        ! Local variables
        integer    :: n, l, np1, nm1
        real(dp)   :: r12(Ndim), r21(Ndim), d1(Ndim), d2(Ndim), fac, prefac, fac1, fac2, C(2,2)
        type(edge) :: l_edge

        ! First set internal forces to zero
        do n = 1,self%number_of_mass_points
            self%mass_points(n)%Fi = 0.0_dp
        end do

        ! Each edge exert an in-plane elastic force on the two mass points connected.
        do l = 1,self%number_of_edges
            
            ! Select the l-th edge (to shorten the notation)
            l_edge = self%edges(l)

            ! Distance vector
            r12 = l_edge%x1%X - l_edge%x2%X
            r21 = l_edge%x2%X - l_edge%x1%X

            ! Elastic in-place forces, eqs (26-27) of de Tullio and Pascazio JCP 2016.
            l_edge%x1%Fi = l_edge%x1%Fi - l_edge%ke*(l_edge%l - l_edge%l0)*r12/l_edge%l
            l_edge%x2%Fi = l_edge%x2%Fi - l_edge%ke*(l_edge%l - l_edge%l0)*r21/l_edge%l
        end do

      ! Each mass node connected to two edges has also blending rigidity
      do n = 1,self%number_of_mass_points
         ! Check if the mass is connected to two edges
         if (self%mass_points(n)%number_of_edges == 2) then

            ! **** Allen & Tildesley, Compute Simulation of Liquids, appendix C ****
            ! Compute distance between mass points
            nm1 = merge(n - 1, self%number_of_mass_points, n > 1)
            np1 = merge(n + 1, 1, n < self%number_of_mass_points)
            d1 = self%mass_points(n  )%X(1:Ndim) - self%mass_points(nm1)%X(1:Ndim)
            d2 = self%mass_points(np1)%X(1:Ndim) - self%mass_points(n  )%X(1:Ndim)

            ! Store the C coefficients in a matrix
            C(1,1) = dot_product(d1,d1)
            C(1,2) = dot_product(d1,d2)
            C(2,1) = dot_product(d2,d1)
            C(2,2) = dot_product(d2,d2)

            prefac = 1.0_dp/sqrt(C(2,2)*C(1,1))
            fac = C(2,1)

            ! The bending potential of mass point n influences mass point n-1, n and n + 1
            fac1 = fac / C(2,2)
            fac2 = fac / C(1,1)
            self%mass_points(np1)%Fi = self%mass_points(np1)%Fi - &
                                       self%kb*prefac*(fac1*d2 - d1)
            self%mass_points( n )%Fi = self%mass_points(n  )%Fi + &
                                       self%kb*prefac*(fac1*d2 - fac2*d1 + d2 - d1)
            self%mass_points(nm1)%Fi = self%mass_points(nm1)%Fi + &
                                       self%kb*prefac*(fac2*d1 - d2)

         endif
      end do

    end subroutine compute_internal_forces
    !========================================================================================

    !========================================================================================
    subroutine integrate_hydrodynamic_forces(self)

        ! Get integral forces on the solid.

        ! In/Out variables
        class(lagrangian_solid), intent(inout) :: self

        ! Local variables
        integer :: l

        ! Hydrodynamic forces
        self%center_of_mass%Fh = 0.0_dp
        self%center_of_mass%Fv = 0.0_dp
        self%center_of_mass%Fp = 0.0_dp

        ! Cycle over lagrangian markers
        do l = 1,self%number_of_edges
            self%center_of_mass%Fh = self%center_of_mass%Fh + self%edges(l)%C%Fh
            self%center_of_mass%Fv = self%center_of_mass%Fv + self%edges(l)%C%Fv
            self%center_of_mass%Fp = self%center_of_mass%Fp + self%edges(l)%C%Fp
        end do

    end subroutine
    !========================================================================================

    !========================================================================================
    subroutine rigid_body_motion(self, central_axis)

        ! This subroutine moves the solid according to the value of the traslation and 
        ! rotation array.

        use IO_mod     , only : print_error_message
        use global_mod , only : Ndim, tdof, rdof

        ! In/Out variables
        class(lagrangian_solid), intent(inout)           :: self
        real(dp)   , intent(in   ), optional :: central_axis(tdof)

        ! Local variables
        integer  :: n, l
        real(dp) :: R(Ndim,Ndim), rot_cen(tdof), d(tdof)

#if DIM==3
#else
        ! Build rotation matrix
        R(1,1) =  cos(self%rot(1))
        R(1,2) = -sin(self%rot(1))
        R(2,1) =  sin(self%rot(1))
        R(2,2) =  cos(self%rot(1))

        ! Set center of rotation
        if (present(central_axis)) then
            rot_cen = central_axis
        else
            rot_cen = self%center_of_mass%X(1:tdof)
        endif

        ! First apply traslation and rotation to each mass point
        do n = 1,self%number_of_mass_points

            ! Apply traslation
            self%mass_points(n)%X = self%mass_points(n)%X + self%tra

            ! Distance from center of rotation
            d = self%mass_points(n)%X - rot_cen

            ! Apply rotation
            self%mass_points(n)%X = rot_cen + matmul(R,d)

        end do

        ! When center of rotation is equal to center of mass, check that the center of mass 
        ! location is preserved
        if (present(central_axis))then
            ! do nothing
        else 
            if (maxval(abs(self%center_of_mass%X - self%get_center_of_mass())) > 1.0e-13_dp) then
                call print_error_message('ERROR: rigid mody motion do not preserve center of mass location')
            endif
        endif

        ! Update lagrangian markers with new mass node position
        call self%update_lagrangian_markers

        ! Now compute velocity and acceleration on each lagrangian marker
        do l = 1,self%number_of_edges

            ! Distance between edge center and rotation center
            d = self%edges(l)%C%X - rot_cen

            ! Velocity of the lagrangian point n is given by
            ! v = v_cm + omega x r
            self%edges(l)%C%V(1) = self%center_of_mass%V(1) - self%center_of_mass%V(3)*d(2)
            self%edges(l)%C%V(2) = self%center_of_mass%V(2) + self%center_of_mass%V(3)*d(1)

            ! The angular velocity is the same
            self%edges(l)%C%V(3) = self%center_of_mass%V(3)

            ! Acceleration of the lagrangian marker l
            ! a = a_cm + omega x (omega x r) + alpha x r
            self%edges(l)%C%A(1) = self%center_of_mass%A(1) - d(1)*self%center_of_mass%V(3)**2 - &
                                    self%center_of_mass%A(3)*d(2)
            self%edges(l)%C%A(2) = self%center_of_mass%A(2) - d(2)*self%center_of_mass%V(3)**2 + &
                                    self%center_of_mass%A(3)*d(1)

            ! The angular acceleration is the same
            self%edges(l)%C%A(3) = self%center_of_mass%A(3)

        end do
#endif

    end subroutine rigid_body_motion
    !==============================================================================================

    !==============================================================================================
    subroutine update_lagrangian_markers(self)

        ! Every time the structures moves it is necessary to update the lagrangian markers
        ! data structure

        ! In/Out variables
        class(lagrangian_solid), intent(inout) :: self

        ! Local variables
        integer :: l

        if (associated(self%apply_constraints)) then
            call self%apply_constraints()
        endif

        do l = 1,self%number_of_edges
            ! Compute edge length
            call self%edges(l)%update_length()

            ! Compute the normal vector corresponding to obj edge
            call self%edges(l)%update_norm()

            ! Compute the centroid
            call self%edges(l)%update_centroid()

            ! Interpolate mass points values on the lagrangian marker
            if (self%is_deformable) call self%edges(l)%get_interpolated_mass_values(self%mass_points)
        end do

    end subroutine update_lagrangian_markers
    !==============================================================================================

    !==============================================================================================
    function get_potential_energy(self) result(Ep)

        use utils_mod, only : clamp

        ! In/Out variables
        class(lagrangian_solid), intent(in) :: self
        real(dp) :: Ep

        ! Local variables
        integer :: l, n, nm1, np1
        real(dp) :: d1(2), d2(2), C(2,2), prefac, fac, pot, theta

        Ep = 0.0_dp
        ! In-plane spring potential energy
        do l = 1,self%number_of_edges
            !Ep = Ep + 0.5_dp*self%edges(l)%ke*(self%edges(l)%l - self%edges(l)%l0)**2
            Ep = Ep + 0.5_dp*self%ke*(self%edges(l)%l - self%edges(l)%l0)**2
        end do

        ! Bending potential
        do n = 1,self%number_of_mass_points
            ! Check if the mass is connected to two edges
            if (self%mass_points(n)%number_of_edges == 2) then

                ! Compute distance between mass points
                nm1 = merge(n - 1, self%number_of_mass_points, n > 1)
                np1 = merge(n + 1, 1, n < self%number_of_mass_points)
                d1 = self%mass_points(n  )%X(1:2) - self%mass_points(nm1)%X(1:2)
                d2 = self%mass_points(np1)%X(1:2) - self%mass_points(n  )%X(1:2)

                ! Store the C coefficients in a matrix
                C(1,1) = dot_product(d1,d1)
                C(1,2) = dot_product(d1,d2)
                C(2,1) = dot_product(d2,d1)
                C(2,2) = dot_product(d2,d2)

                prefac = 1.0_dp/sqrt(C(2,2)*C(1,1))
                fac = C(2,1)
                !pot = -self%kb*prefac*fac ! This is -cos(theta)
                pot = -prefac*fac ! This is -cos(theta)
                theta = acos(clamp(-pot, -1.0_dp, 1.0_dp))
                Ep = Ep + self%kb*(1.0_dp - cos(theta))
            end if
        end do

    end function get_potential_energy
    !==============================================================================================

    !==============================================================================================
    function get_kinetic_energy(self) result(Ek)

        ! In/Out variables
        class(lagrangian_solid), intent(in) :: self
        real(dp) :: Ek

        ! Local variables
        integer :: n

        Ek = 0.0_dp
        do n = 1,self%number_of_mass_points
            Ek = Ek + 0.5_dp*(self%mass_points(n)%m*self%mass_points(n)%V(1)**2 + &
                self%mass_points(n)%m*self%mass_points(n)%V(2)**2)
        end do

    end function get_kinetic_energy
    !==============================================================================================

    !==============================================================================================
    subroutine velocity_verlet(self, dt)

        ! Advance the solution array with the velocity verlet scheme

        ! In/Out variables
        class(lagrangian_solid), intent(inout) :: self
        real(dp)    , intent(in   ) :: dt

        ! Local variables
        integer :: n

        ! Update internal forces
        call self%compute_internal_forces()

        do n = 1,self%number_of_mass_points
            self%mass_points(n)%V = self%mass_points(n)%V + &
                0.5_dp*dt*(self%mass_points(n)%Fe + self%mass_points(n)%Fi)/self%mass_points(n)%m
            self%mass_points(n)%X = self%mass_points(n)%X + dt*self%mass_points(n)%V
        end do

        call self%update_lagrangian_markers
        call self%compute_internal_forces()

        do n = 1,self%number_of_mass_points
            self%mass_points(n)%V = self%mass_points(n)%V*(1.0_dp - self%c) + &
                0.5_dp*dt*(self%mass_points(n)%Fi + self%mass_points(n)%Fe)/self%mass_points(n)%m
        end do

        call self%update_lagrangian_markers

    end subroutine velocity_verlet
    !==============================================================================================

    !==============================================================================================
    subroutine write_csv(self, time)

        ! In/Out variables
        class(lagrangian_solid), intent(in) :: self
        real(dp)    , intent(in) :: time

        ! If writing the first line write first the header
        if (time == 0.0_dp) then
            write(self%output_file_id,11) 't,x,y,theta,u,v,omega,ax,ay,alpha,Fsx,Fsy,Fpx,Fpy,Fx,Fy,Mz'
        endif
11      format(A58)

        write(self%output_file_id,12) time, ",", &
            self%center_of_mass%X(1) , ",", self%center_of_mass%X(2) , ",", self%center_of_mass%X(3), ",", &
            self%center_of_mass%V(1) , ",", self%center_of_mass%V(2) , ",", self%center_of_mass%V(3), ",", &
            self%center_of_mass%A(1) , ",", self%center_of_mass%A(2) , ",", self%center_of_mass%A(3), ",", &
            self%center_of_mass%Fv(1), ",", self%center_of_mass%Fv(2), ",", &
            self%center_of_mass%Fp(1), ",", self%center_of_mass%Fp(2), ",", &
            self%center_of_mass%Fh(1), ",", self%center_of_mass%Fh(2), ",", self%center_of_mass%Fh(3)
12      format(E16.8,16(A1,E16.8))

    end subroutine write_csv
    !==============================================================================================

    !==============================================================================================
    subroutine print_configuration(self, step)

        ! In/Out variables
        class(lagrangian_solid), intent(in) :: self
        integer     , intent(in) :: step

        ! Local vraiables
        integer :: out_id, l
        character(len=7) :: sn
        character(len=15) :: filename

        ! Open output file
        write(sn,'(I0.7)') step
        filename = 'data/sb_'//sn

        open(newunit = out_id, file = filename)
        do l = 1,self%number_of_edges
            write(out_id,*) self%edges(l)%x1%X
            write(out_id,*) self%edges(l)%x2%X
            write(out_id,*) ''
        end do
        close(out_id)

    end subroutine print_configuration
    !==============================================================================================

    !==============================================================================================
    subroutine destroy(self)

        class(lagrangian_solid), intent(inout) :: self

        deallocate(self%edges)
        deallocate(self%mass_points)
        call self%center_of_mass%destroy

    end subroutine destroy
    !==============================================================================================

end module
