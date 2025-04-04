!> ! Define the lagrangian solid objects for the Immersed Boundary Method
module lagrangian_solid_1D_mod

    use precision_mod
    use global_mod
    use marker_mod
    use edge_mod
    use lagrangian_solid_mod, only : lagrangian_solid

    implicit none
    private
    public :: lagrangian_solid_1D

    ! Define a solid object
    type, extends(lagrangian_solid) :: lagrangian_solid_1D
        integer, allocatable :: mass_point_edges_index(:,:)
        integer, allocatable :: edge_mass_points_index(:,:)

        real(dp)                      :: central_axis(tdof)
        real(dp)                      :: kb                      !< elastic constant for bending 
        real(dp)                      :: c = 0.0_dp              !< damping for the structural solver
        
        integer                       :: output_file_id          !< id number for the output file
        logical                       :: is_out = .false.        !< flag for checking if the solid is outside the domain

    contains
        procedure :: create
        procedure :: update
        procedure :: get_center_of_mass

        procedure :: compute_internal_forces
        procedure :: interpolate_from_forcing_element_to_mass_point
        procedure :: update_forcing_elements
 
        procedure :: print_configuration
        procedure :: get_total_length
        procedure :: getPotentialEnergy => get_potential_energy
        procedure :: get_kinetic_energy
        procedure :: destroy

    end type lagrangian_solid_1D

    abstract interface
        subroutine constraints(self)
            import lagrangian_solid_1D
            class(lagrangian_solid_1D), intent(inout) :: self
        end subroutine constraints
    end interface

contains

    !===============================================================================================
    subroutine create(self, filename, name)

        ! This subroutine initialize the variables of the solid.
        ! The file filename must contains the location of the mass points.

        ! In/Out variables
        class(lagrangian_solid_1D), intent(inout), target   :: self
        character(len=*)          , intent(in   )           :: filename
        character(len=*)          , intent(in   ), optional :: name

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
            self%mass_points(n) = marker(tdof, n)

            ! Read the position from file
            read(fid,*) self%mass_points(n)%X

            ! Set the discrete mass of the point
            self%mass_points(n)%m = self%mass/self%number_of_mass_points
        end do
        close(fid)

        ! **** Setup Center of Mass ****
        ! Compute center of mass of the solid body
        self%center_of_mass = marker(dofs, 0)
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

        ! Allocate connection matrix
        allocate(self%mass_point_edges_index(self%number_of_mass_points,2))
        allocate(self%edge_mass_points_index(self%number_of_edges,2))

        ! Initialize all edges
        do l = 1,self%number_of_edges - 1
            ! Each edge connect two masses, select the indexes
            self%edge_mass_points_index(l,:) = [l, l + 1]
            
            ! Each edge connect two masses, select the indexes
            self%edges(l)%m1 => self%mass_points(l)
            self%edges(l)%m2 => self%mass_points(l+1)
        end do

        ! Fix the last edge
        l = self%number_of_edges
        if (self%is_open) then
            self%edges(l)%m1 => self%mass_points(l)
            self%edges(l)%m2 => self%mass_points(l+1)
            self%edge_mass_points_index(l,:) = [l, l + 1]
        else
            self%edges(l)%m1 => self%mass_points(l)
            self%edges(l)%m2 => self%mass_points(1)
            self%edge_mass_points_index(l,:) = [l, 1]
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
            self%mass_point_edges_index(n,:) = [n - 1, n]
        end do
        self%mass_point_edges_index(1,1) = self%number_of_edges
        self%mass_point_edges_index(self%number_of_mass_points,2) = 1

        ! If the solid is open, the first and last mass point have only one edge, fix it.
        if (self%is_open) then
            self%mass_point_edges_index(1,1) = -1
            self%mass_point_edges_index(self%number_of_mass_points,2) = -1
        end if

        ! Setup the forcing elements if using Immersed boundary method
#ifdef IBM
        allocate(self%forcing_elements(self%number_of_edges))
        do n = 1,self%number_of_edges
            self%forcing_elements(n)%C => self%edges(n)%C
            self%forcing_elements(n)%A => self%edges(n)%l
            self%forcing_elements(n)%n => self%edges(n)%n
        end do
        allocate(self%p_s(self%number_of_edges))
        allocate(self%tau_s(6,self%number_of_edges))
#endif

        if (present(name)) self%name = name

    end subroutine create
    !===============================================================================================

    !===============================================================================================
    subroutine update(self)

        ! Objective update the lagrangian markers data structure

        ! In/Out variables
        class(lagrangian_solid_1D), intent(inout), target :: self

        ! Local variables
        integer :: l

        do l = 1,self%number_of_edges
            ! Compute edge length
            call self%edges(l)%update_length()

            ! Compute the normal vector corresponding to obj edge
            call self%edges(l)%update_norm()

            ! Compute the centroid
            call self%edges(l)%update_centroid()
        end do

    end subroutine update
    !===============================================================================================

    !===============================================================================================
    function get_center_of_mass(self) result(Xcm)

        ! Find the center of mass location of the solid.

        ! In/Out variables
        class(lagrangian_solid_1D), intent(in) :: self

        ! Local variables
        integer  :: n
        real(dp) :: Xcm(tdof)

        Xcm = 0.0_dp
        do n = 1,self%number_of_mass_points
            Xcm = Xcm + self%mass_points(n)%X
        end do

        Xcm = Xcm/real(self%number_of_mass_points, dp)

    end function
    !===============================================================================================

    !===============================================================================================
    function get_total_length(self) result(L)

        ! In/Out variables
        class(lagrangian_solid_1D), intent(in) :: self

        ! Local variables
        integer  :: n
        real(dp) :: L

        L = 0.0_dp
        do n = 1,self%number_of_edges
            L = L + self%edges(n)%l
        end do

    end function
    !===============================================================================================
    
    !===============================================================================================
    subroutine interpolate_from_forcing_element_to_mass_point(self)

        ! Objective: distribute the forces from the Lagrangian markers to the mass points. The force
        !            on the mass point is given by the weigthed average between the froces on the
        !            lagragian markers using as weight the size of the edges.
    
        use global_mod, only : Ndim

        ! In/out variables
        class(lagrangian_solid_1D), intent(inout), target :: self

        ! Local variables
        integer  :: n

        do n = 1,self%number_of_mass_points
            self%mass_points(n)%Fh = 0.0_dp
        end do

        do n = 1,self%number_of_edges
            self%edges(n)%m1%Fh = self%edges(n)%m1%Fh + 0.5_dp*self%edges(n)%C%Fh(1:Ndim)
            self%edges(n)%m2%Fh = self%edges(n)%m2%Fh + 0.5_dp*self%edges(n)%C%Fh(1:Ndim)
        end do

    end subroutine
    !===============================================================================================

    !===============================================================================================
    subroutine update_forcing_elements(self)

        ! Objective: compute velocity and acceleration on the forcing element center.
        !            For deformable objects, distribute velocity and acceleration from mass points;
        !            for rigid objecst, get information from center of mass.

        use global_mod, only : Ndim

        ! In/out variables
        class(lagrangian_solid_1D), intent(inout), target :: self

        ! Local variables
        integer  :: n

        if (self%is_deformable) then
            ! Initialize to zero
            do n = 1,self%number_of_edges
                self%edges(n)%C%V(1:2) = 0.5_dp*(self%edges(n)%m1%V(1:2) + self%edges(n)%m2%V(1:2))
                self%edges(n)%C%A(1:2) = 0.5_dp*(self%edges(n)%m1%A(1:2) + self%edges(n)%m2%A(1:2))
            end do
        else

        endif

    end subroutine
    !===============================================================================================
    
    !===============================================================================================
    subroutine compute_internal_forces(self)

        ! In/Out variables
        class(lagrangian_solid_1D), intent(inout), target :: self

        ! Local variables
        integer    :: n, l, np1, nm1, n1, n2
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
            r12 = l_edge%m1%X - l_edge%m2%X
            r21 = l_edge%m2%X - l_edge%m1%X

            ! Elastic in-place forces, eqs (26-27) of de Tullio and Pascazio JCP 2016.
            l_edge%m1%Fi = l_edge%m1%Fi - l_edge%ks*(l_edge%l - l_edge%l0)*r12/l_edge%l
            l_edge%m2%Fi = l_edge%m2%Fi - l_edge%ks*(l_edge%l - l_edge%l0)*r21/l_edge%l
        end do

        ! Each mass node connected to two edges has also blending rigidity
        if (self%is_open) then
            n1 = 2
            n2 = self%number_of_mass_points - 1
        else
            n1 = 1
            n2 = self%number_of_mass_points
        endif
        do n = n1,n2
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
            self%mass_points(np1)%Fi = self%mass_points(np1)%Fi - self%kb*prefac*(fac1*d2 - d1)
            self%mass_points( n )%Fi = self%mass_points(n  )%Fi + self%kb*prefac*(fac1*d2 - &
                                            fac2*d1 + d2 - d1)
            self%mass_points(nm1)%Fi = self%mass_points(nm1)%Fi + self%kb*prefac*(fac2*d1 - d2)
        end do

    end subroutine compute_internal_forces
    !===============================================================================================

    !===============================================================================================
    function get_potential_energy(self) result(Ep)

        use utils_mod, only : clamp

        ! In/Out variables
        class(lagrangian_solid_1D), intent(in), target :: self
        real(dp) :: Ep

        ! Local variables
        integer :: l, n, nm1, np1, n1, n2
        real(dp) :: d1(2), d2(2), C(2,2), prefac, fac, pot, theta

        Ep = 0.0_dp
        ! In-plane spring potential energy
        do l = 1,self%number_of_edges
            Ep = Ep + 0.5_dp*self%edges(l)%ks*(self%edges(l)%l - self%edges(l)%l0)**2
            !Ep = Ep + 0.5_dp*self%ke*(self%edges(l)%l - self%edges(l)%l0)**2
        end do

        ! Each mass node connected to two edges has also blending rigidity
        if (self%is_open) then
            n1 = 2
            n2 = self%number_of_mass_points - 1
        else
            n1 = 1
            n2 = self%number_of_mass_points
        endif

        ! Bending potential
        do n = n1, n2
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
        end do

    end function get_potential_energy
    !===============================================================================================

    !===============================================================================================
    function get_kinetic_energy(self) result(Ek)

        ! In/Out variables
        class(lagrangian_solid_1D), intent(in) :: self
        real(dp) :: Ek

        ! Local variables
        integer :: n

        Ek = 0.0_dp
        do n = 1,self%number_of_mass_points
            Ek = Ek + 0.5_dp*(self%mass_points(n)%m*self%mass_points(n)%V(1)**2 + &
                              self%mass_points(n)%m*self%mass_points(n)%V(2)**2)
        end do

    end function get_kinetic_energy
    !===============================================================================================

    !===============================================================================================
    subroutine print_configuration(self, step)

        ! In/Out variables
        class(lagrangian_solid_1D), intent(in) :: self
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
            write(out_id,*) self%edges(l)%m1%X
            write(out_id,*) self%edges(l)%m2%X
            write(out_id,*) ''
        end do
        close(out_id)

    end subroutine print_configuration
    !===============================================================================================

    !===============================================================================================
    subroutine destroy(self)

        class(lagrangian_solid_1D), intent(inout) :: self

        deallocate(self%edges)
        deallocate(self%mass_points)
        deallocate(self%mass_point_edges_index)
        deallocate(self%edge_mass_points_index)
        call self%center_of_mass%destroy
#ifdef IBM
        deallocate(self%forcing_elements)
#endif

    end subroutine destroy
    !===============================================================================================

end module
