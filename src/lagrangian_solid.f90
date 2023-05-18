!> ! Define the lagrangian solid objects for the Immersed Boundary Method
module lagrangian_solid

    use precision, only : dp

    implicit none

    type marker
        !< Base type to define a point
        real(dp), dimension(:), allocatable :: X  !< Position vector
        real(dp), dimension(:), allocatable :: V  !< Velocity vector
        real(dp), dimension(:), allocatable :: A  !< Acceleration vector
        real(dp), dimension(:), allocatable :: F  !< Generic  force vector
        real(dp), dimension(:), allocatable :: Fe !< External force vector
        real(dp), dimension(:), allocatable :: Fv !< Viscous  force vector
        real(dp), dimension(:), allocatable :: Fp !< Pressure force vector
        real(dp), dimension(:), allocatable :: Fi !< Internal force vector
    end type marker

    ! Define the mass point type
    type, extends(marker) ::  mass_point
        !< A mass point is a marker with a mass connected to one or more edges.                      
        real(dp)                           :: m               !< Discrete inertial mass
        integer                            :: number_of_edges !< self exp.
        integer, dimension(:), allocatable :: edges_index     !< self exp.
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
    type solid
        integer                       :: number_of_mass_points   !< self exp.
        integer                       :: number_of_edges         !< self exp.
        type(mass_point), allocatable :: mass_points(:)          !< array of mass points
        type(edge), allocatable       :: edges(:)                !< array of edges
        real(dp)                      :: Vol                     !< volume of the solid
        real(dp)                      :: rho                     !< density of the solid
        real(dp)                      :: M(3) = 1.0_dp           !< inertial mass (applied to the center of mass)
        type(marker)                  :: center_of_mass          !< self exp.
        real(dp)                      :: tra(2)                  !< rigid body traslation
        real(dp)                      :: rot(1)                  !< rigid body rotation
        integer                       :: nsubsteps = 200         !< number of substep for the structural solver
        real(dp)                      :: ke                      !< elastic constant in-plane deformation
        real(dp)                      :: kb                      !< elastic constant for bending 
        real(dp)                      :: c = 0.0_dp              !< damping for the structural solver
        integer                       :: outpit_file_id          !< id number for the output file
        logical                       :: is_deformable = .false. !< flag for solving deformation
        logical                       :: is_open = .false.       !< flag for open structure
        logical                       :: is_out = .false.        !< flag for checking if the solid is outside the domain  
    contains
        procedure :: create                       
        procedure :: get_center_of_mass
        procedure :: compute_internal_forces
        procedure :: integrate_hydrodynamic_forces
        procedure :: update_lagrangian_markers
        procedure :: print_configuration
        procedure :: get_total_length
        procedure :: get_potential_energy
        procedure :: get_kinetic_energy
        procedure :: velocity_verlet
        procedure :: write_csv
        procedure :: destroy
    end type solid

    interface
        subroutine constraints(self)
            import solid
            type(solid), intent(inout) :: self
        end subroutine constraints
    end interface
    procedure(constraints), pointer :: apply_constraints => Null()

contains

    !========================================================================================
    subroutine create(self, filename)

        ! This subroutine initialize the variables of the solid.
        ! The file filename must contains the location of the mass points.

        ! In/Out variables
        class(solid)    , intent(inout), target :: self
        character(len=*), intent(in   )         :: filename

        ! Local variables
        integer  :: fid, nl, io, n, l
        real(dp) :: Xcm(2)

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
            allocate(self%mass_points(n)%X(2))
            read(fid,*) self%mass_points(n)%X(1), self%mass_points(n)%X(2)

            ! Set other variables to zero
            allocate(self%mass_points(n)%V(2))
            self%mass_points(n)%V = 0.0_dp
            allocate(self%mass_points(n)%A(2))
            self%mass_points(n)%A = 0.0_dp
            allocate(self%mass_points(n)%Fe(2))
            self%mass_points(n)%Fe = 0.0_dp
            allocate(self%mass_points(n)%Fi(2))
            self%mass_points(n)%Fi = 0.0_dp

            ! Set the discrete mass of the point
            self%mass_points(n)%m = self%M(1)/self%number_of_mass_points
        end do
        close(fid)

        ! **** Setup Center of Mass ****
        ! Compute center of mass of the solid body
        allocate(self%center_of_mass%X(3))
        allocate(self%center_of_mass%V(3))
        allocate(self%center_of_mass%A(3))
        allocate(self%center_of_mass%F(3))
        allocate(self%center_of_mass%Fv(3))
        allocate(self%center_of_mass%Fp(3))
        Xcm = self%get_center_of_mass()
        self%center_of_mass%X(1:2) = [Xcm(1), Xcm(2)]
        
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
        do l = 1,self%number_of_edges
            ! Each edge connect two masses, select the indexes
            self%edges(l)%mass_point_index = [l, l + 1]
            ! Each edge connect two masses, select the indexes
            self%edges(l)%x1 => self%mass_points(l)
            self%edges(l)%x2 => self%mass_points(l+1)
        end do

        ! Fix the last edge if the structure is closed
        if (self%is_open .eqv. .false.) then
            l = self%number_of_edges
            self%edges(l)%mass_point_index(2) = 1
        end if

        ! Setup edge variables
        do l = 1,self%number_of_edges
            ! Compute edge length
            call self%edges(l)%update_length()

            ! Set the initial equilibrium length
            self%edges(l)%l0 = self%edges(l)%l

            ! Compute the normal vector
            call self%edges(l)%update_norm()

            allocate(self%edges(l)%C%X(3))
            allocate(self%edges(l)%C%V(3))
            allocate(self%edges(l)%C%A(3))
            allocate(self%edges(l)%C%Fv(3))
            allocate(self%edges(l)%C%Fp(3))
            allocate(self%edges(l)%C%F(3))

            ! Compute the centroid
            call self%edges(l)%update_centroid()

            ! Set velocity, acceleration and forces to zero
            self%edges(l)%C%V = 0.0_dp
            self%edges(l)%C%A = 0.0_dp
            self%edges(l)%C%F = 0.0_dp
            self%edges(l)%C%Fv = 0.0_dp
            self%edges(l)%C%Fp = 0.0_dp
        end do

        ! Set the index of the mass points. By default consider the structure closed
        do n = 1,self%number_of_mass_points
            self%mass_points(n)%number_of_edges = 2
            !allocate(self%mass_points(n)%edges_index(2))
            !self%mass_points(n)%edges_index = [n - 1, n]
        end do
        !self%mass_points(1)%edges_index(1) = self%number_of_edges
        !self%mass_points(self%number_of_mass_points)%edges_index(2) = 1

        ! If the solid is open, the first and last mass point have only one edge, fix it.
        if (self%is_open) then
            !deallocate(self%mass_points(1)%edges_index)
            self%mass_points(1)%number_of_edges = 1
            !allocate(self%mass_points(1)%edges_index(1))
            !self%mass_points(1)%edges_index(1) = 1
            !deallocate(self%mass_points(self%number_of_mass_points)%edges_index)
            self%mass_points(self%number_of_mass_points)%number_of_edges = 1
            !allocate(self%mass_points(self%number_of_mass_points)%edges_index(1))
            !self%mass_points(self%number_of_mass_points)%edges_index(1) = self%number_of_edges
        end if

    end subroutine create
    !========================================================================================

    !========================================================================================
    function get_center_of_mass(self) result(Xcm)

        ! Find the center of mass location of the solid.

        ! In/Out variables
        class(solid), intent(in) :: self

        ! Local variables
        integer  :: n
        real(dp) :: Xcm(2)

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
        mod_t = sqrt(t(1)**2 + t(2)**2 + 1.0e-16_dp)
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
        class(solid), intent(in) :: self

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
        class(solid), intent(inout) :: self

        ! Local variables
        integer    :: n, l, i_mp1, i_mp2, np1, nm1
        real(dp)   :: r12(2), r21(2), d1(2), d2(2), fac, prefac, fac1, fac2, C(2,2)
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
            l_edge%x1%Fi = l_edge%x1%Fi - self%ke*(l_edge%l - l_edge%l0)*r12/l_edge%l
            l_edge%x2%Fi = l_edge%x2%Fi - self%ke*(l_edge%l - l_edge%l0)*r21/l_edge%l
        end do

      ! Each mass node connected to two edges has also blending rigidity
      do n = 1,self%number_of_mass_points
         ! Check if the mass is connected to two edges
         if (self%mass_points(n)%number_of_edges == 2) then

            ! **** Allen & Tildesley, Compute Simulation of Liquids, appendix C ****
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
        class(solid), intent(inout) :: self

        ! Local variables
        integer :: l

        ! Hydrodynamic forces
        self%center_of_mass%F  = 0.0_dp
        self%center_of_mass%Fv = 0.0_dp
        self%center_of_mass%Fp = 0.0_dp

        ! Cycle over lagrangian markers
        do l = 1,self%number_of_edges
            self%center_of_mass%F  = self%center_of_mass%F  + self%edges(l)%C%F
            self%center_of_mass%Fv = self%center_of_mass%Fv + self%edges(l)%C%Fv
            self%center_of_mass%Fp = self%center_of_mass%Fp + self%edges(l)%C%Fp
        end do

    end subroutine
    !========================================================================================

    !========================================================================================
    subroutine update_lagrangian_markers(self)

        ! Every time the structures moves it is necessary to update the lagrangian markers
        ! data structure

        ! In/Out variables
        class(solid), intent(inout) :: self

        ! Local variables
        integer :: l

        if (associated(apply_constraints)) then
            call apply_constraints(self)
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
    !========================================================================================

    !========================================================================================
    function get_potential_energy(self) result(Ep)

        use utils, only : clamp

        ! In/Out variables
        class(solid), intent(in) :: self
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
    !========================================================================================

    !========================================================================================
    function get_kinetic_energy(self) result(Ek)

        ! In/Out variables
        class(solid), intent(in) :: self
        real(dp) :: Ek

        ! Local variables
        integer :: n

        Ek = 0.0_dp
        do n = 1,self%number_of_mass_points
            Ek = Ek + 0.5_dp*(self%mass_points(n)%m*self%mass_points(n)%V(1)**2 + &
                self%mass_points(n)%m*self%mass_points(n)%V(2)**2)
        end do

    end function get_kinetic_energy
    !========================================================================================

    !========================================================================================
    subroutine velocity_verlet(self, dt)

        ! Advance the solution array with the velocity verlet scheme

        ! In/Out variables
        class(solid), intent(inout) :: self
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
    !========================================================================================

    !========================================================================================
    subroutine write_csv(self, fid, time)

        ! In/Out variables
        class(solid), intent(in) :: self
        integer     , intent(in) :: fid
        real(dp)    , intent(in) :: time

        ! If writing the first line write first the header
        if (time == 0.0_dp) then
            write(fid,11) 't,x,y,theta,u,v,omega,ax,ay,alpha,Fsx,Fsy,Fpx,Fpy,Fx,Fy,Mz'
        endif
11      format(A58)

        write(fid,12) time, ",", &
            self%center_of_mass%X(1) , ",", self%center_of_mass%X(2) , ",", self%center_of_mass%X(3), ",", &
            self%center_of_mass%V(1) , ",", self%center_of_mass%V(2) , ",", self%center_of_mass%V(3), ",", &
            self%center_of_mass%A(1) , ",", self%center_of_mass%A(2) , ",", self%center_of_mass%A(3), ",", &
            self%center_of_mass%Fv(1), ",", self%center_of_mass%Fv(2), ",", &
            self%center_of_mass%Fp(1), ",", self%center_of_mass%Fp(2), ",", &
            self%center_of_mass%F(1) , ",", self%center_of_mass%F(2) , ",", self%center_of_mass%F(3)
12      format(E16.8,16(A1,E16.8))

    end subroutine write_csv
    !========================================================================================

    !========================================================================================
    subroutine print_configuration(self, step)

        ! In/Out variables
        class(solid), intent(in) :: self
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
            write(out_id,*) self%mass_points(self%edges(l)%mass_point_index(1))%X
            write(out_id,*) self%mass_points(self%edges(l)%mass_point_index(2))%X
            write(out_id,*) ''
        end do
        close(out_id)

    end subroutine print_configuration
    !========================================================================================

    !========================================================================================
    subroutine destroy(self)

        class(solid), intent(inout) :: self

        deallocate(self%edges)
        deallocate(self%mass_points)
        deallocate(self%center_of_mass%X)
        deallocate(self%center_of_mass%V)
        deallocate(self%center_of_mass%A)
        deallocate(self%center_of_mass%F)
        deallocate(self%center_of_mass%Fv)
        deallocate(self%center_of_mass%Fp)

    end subroutine destroy
    !========================================================================================

end module lagrangian_solid
