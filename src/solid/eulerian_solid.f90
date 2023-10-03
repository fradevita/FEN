!> ! Define the eulerian solid objects for the Immersed Boundary Method
module eulerian_solid_mod

    use precision_mod, only : dp
    use global_mod
    use grid_mod     , only : grid
    use solid_mod
    use marker_mod

    implicit none
    private
    public :: eulerian_solid
    public :: eulerian_solid_pointer

    !> Base object for eulerian IBM
    type, extends(solid), abstract :: eulerian_solid
        real(dp) :: mass = 1.0_dp                      !< mass of the solid
        real(dp) :: M(6) = 1.0_dp                      !< inertial mass
        real(dp) :: rotation(3) = 0.0_dp               !< center of mass rotation
        real(dp) :: traslation(3) = 0.0_dp             !< center of mass traslation
        real(dp) :: rot_center(3) = 0.0_dp             !< center of rotation
        type(grid), pointer :: G => Null()             !< grid where the solid is defined
        type(marker), allocatable :: surface_points(:) !< array of surface points
    contains
        ! List of object-dependent procedure
        procedure(distance_interface          ), pass(self), deferred :: distance           !< distance from the solid surface
        procedure(norm_interface              ), pass(self), deferred :: norm               !< local norm vector
        procedure(volume_interface            ), pass(self), deferred :: volume             !< volume (area in 2D) of the solid
        procedure(rotational_inertia_interface), pass(self), deferred :: rotational_inertia !< rotational inertia around center of mass
        ! List of object-independent procedure
        procedure, pass(self) :: setup
        procedure, pass(self) :: advance               !< Timestep advancing of the solid
        procedure, pass(self) :: load_surface_points   !< Initialize surface points for the probe
        procedure, pass(self) :: update_surface_points !< Update position of the surface points
        procedure, pass(self) :: velocity              !< Rigid body velocity
        procedure, pass(self) :: acceleration          !< Rigid body acceleration
        procedure, pass(self) :: check_periodicity
        !procedure, pass(self) :: print_csv             !< Output solid center of mass data
        procedure, pass(self) :: print_surface_points  !< Output surface points
    end type eulerian_solid

    abstract interface
        !< Abastract interface for eulerian solid object
        function distance_interface(self, x) result(f)
            !< Function to compute distance from solid surface at given location x.
            use precision_mod, only : dp
            import eulerian_solid
            class(eulerian_solid), intent(in) :: self    !< eulerian solid object
            real(dp)             , intent(in) :: x(3)    !< location where to compute distance
            real(dp)                          :: f       !< distance result
        end function distance_interface

        function norm_interface(self, x) result(n)
            !< Function to compute local norm with respect to solid surface at given location x.
            use precision_mod, only : dp
            import :: eulerian_solid
            class(eulerian_solid), intent(in) :: self !< eulerian solid object
            real(dp)             , intent(in) :: x(3) !< location where to compute normal vector
            real(dp)                          :: n(3) !< local norm vector
        end function norm_interface

        function volume_interface(self) result(V)
            !< Function to compute local norm with respect to solid surface at given location x.
            use precision_mod, only : dp
            import :: eulerian_solid
            class(eulerian_solid), intent(in) :: self !< eulerian solid object
            real(dp)                          :: V    !< volume of the solid
        end function volume_interface

        function rotational_inertia_interface(self) result(I)
            !< Function to compute local norm with respect to solid surface at given location x.
            use precision_mod, only : dp
            import :: eulerian_solid
            class(eulerian_solid), intent(in) :: self !< eulerian solid object
            real(dp)                          :: I(3) !< rotational inertia 
        end function rotational_inertia_interface
    end interface

    type eulerian_solid_pointer
        class(eulerian_solid), pointer :: pS => Null()
    end type eulerian_solid_pointer

contains

    !==============================================================================================
    subroutine setup(self)

        ! In/Out variables
        class(eulerian_solid), intent(inout) :: self

        ! Local variables
        character(len=99) :: outfile
    
        ! Evaluate solid mass
        self%mass = self%density*self%volume()

        ! Build inertial mass matrix
        self%M = [self%mass, self%mass, self%mass, [self%rotational_inertia()]]

        ! Create the center of mass marker
        call self%center_of_mass%create(6)

        ! Set the rotation center equal to the center of mass
        self%rot_center = self%center_of_mass%X(1:3)

        ! Open output file
        outfile = './data/'//trim(self%name)
        open(newunit = self%file_id, file = trim(outfile))
        write(self%file_id,11) 't,x,y,z,alpha,beta,gamma,u,v,w,omegax,omegay,omegaz,ax,&
                               &ay,az,alphax,alphay,alphaz,Fx,Fy,Fz,Mx,My,Mz'
11    format(A99)

    end subroutine setup
    !==============================================================================================
    
    !==============================================================================================
    subroutine advance(self, dt)

        ! In/Out variables
        class(eulerian_solid), intent(inout), target :: self
        real(dp)             , intent(in   )         :: dt

        ! Local variables
        real(dp), dimension(6) :: Vnp1, Xnp1
        type(marker), pointer  :: CoM

        CoM => self%center_of_mass

        ! Acceleration
        CoM%A = (CoM%Fh + CoM%Fe)/self%M
        
        ! New velocity and position
        Vnp1 = CoM%V + dt*CoM%A
        Xnp1 = CoM%X + dt*0.5_dp*(Vnp1 + CoM%V)

        ! Evaluate rotation and traslation
        self%traslation = Xnp1(1:3) - CoM%X(1:3)
        self%rotation   = Xnp1(4:6) - CoM%X(4:6)

        ! Update solution
        CoM%V = Vnp1
        CoM%X = Xnp1

        ! Update rotation center
        self%rot_center = CoM%X(1:3)

        ! Update surface points position
        call self%update_surface_points()

    end subroutine advance
    !==============================================================================================
   
    !==============================================================================================
    subroutine load_surface_points(self, mesh_file)

        ! When using the probe method to evaluate hydrodynamic forces, must load surface
        ! points on the solid.
  
        ! In/Out variables
        class(eulerian_solid), intent(inout) :: self
        character(len=*)     , intent(in   ) :: mesh_file

        ! Local variables
        integer  :: mesh_file_id, nl, status, n
        real(dp) :: X(3)

        open(newunit = mesh_file_id, file = mesh_file, action = 'read')
        nl = 0
        do
            read(mesh_file_id, *, iostat = status)
            if (status /= 0) exit
            nl = nl + 1
        end do
        rewind(mesh_file_id)

        allocate(self%surface_points(nl))
        do n = 1,nl
            call self%surface_points(n)%create(6)
            read(mesh_file_id,*) X(1), X(2), X(3)
            self%surface_points(n)%X = X
        end do
        close(mesh_file_id)

    end subroutine load_surface_points
    !==============================================================================================
  
    !==============================================================================================
    subroutine update_surface_points(self)
  
        use utils_mod, only : Rx, Ry, Rz

        ! In/Out variables
        class(eulerian_solid), intent(inout) :: self

        ! Local variables
        integer :: n
        real(dp) :: R(3,3), d(3)

        ! Build rotation matrix
        R = matmul(matmul(Rz(self%rotation(3)),Ry(self%rotation(2))),Rx(self%rotation(1)))
        
        ! Apply traslation and rotation to each mass point
        do n = 1,size(self%surface_points)
            ! Apply traslation:
            self%surface_points(n)%X = self%surface_points(n)%X + self%traslation

            ! Distance from rotation center
            d = self%surface_points(n)%X - self%rot_center

            ! Apply rotation:
            self%surface_points(n)%X = self%rot_center + matmul(R,d)
        end do
  
    end subroutine update_surface_points
    !==============================================================================================

    !==============================================================================================
    subroutine print_surface_points(self, step)
  
        ! In/Out variables
        class(eulerian_solid), intent(in) :: self
        integer              , intent(in) :: step

        ! Local variables
        integer           :: n, out_id
        character(len=7 ) :: sn
        character(len=99) :: outfile 

        ! Output file
        write(sn,'(I0.7)') step
        outfile = 'data/'//trim(self%name)//'_'//sn
        open(newunit = out_id, file = trim(outfile))
        do n = 1,size(self%surface_points)
            write(out_id, *) self%surface_points(n)%X
        end do
        close(out_id)

    end subroutine print_surface_points
    !==============================================================================================
    
    !==============================================================================================
    function velocity(self, X, dir) result(v)

        ! Return dir component of velocity due to rigid body motion in X
        ! Velocity is given by the sum of the center of mass velocity and the rotation
        ! contribution.

        use precision_mod, only : dp
        use utils_mod    , only : cross_product

        class(eulerian_solid), intent(in) :: self
        real(dp)             , intent(in) :: X(3)
        integer              , intent(in) :: dir
        real(dp)                          :: v

        ! Local variable
        integer  :: i
        real(dp) :: d(3), r(3), omega_vec_r(3)

        ! Position vector of the point X with respect to the rotation center
        d(1) = sqrt( (X(1) - self%rot_center(1)            )**2 + (X(2) - self%rot_center(2))**2)
        d(2) = sqrt( (X(1) - self%rot_center(1) + self%G%Lx)**2 + (X(2) - self%rot_center(2))**2)
        d(3) = sqrt( (X(1) - self%rot_center(1) - self%G%Lx)**2 + (X(2) - self%rot_center(2))**2)
        i = minloc(d, 1)
        if (i == 1) then
            r = X - self%rot_center
        elseif (i == 2) then
            r = X + [self%G%Lx, 0.0_dp, 0.0_dp] - self%rot_center
        elseif (i == 3) then
            r = X - [self%G%Lx, 0.0_dp, 0.0_dp] - self%rot_center
        endif

        ! Velocity due to rotation
        omega_vec_r = cross_product(self%center_of_mass%V(4:6), r)

        ! Total velocity
        v = self%center_of_mass%V(dir) + omega_vec_r(dir)

    end function velocity
    !==============================================================================================

    !==============================================================================================
    function acceleration(self, X, dir) result(a)

        ! Return dir component of acceleration due to rigid body motion in X
        ! Velocity is given by the sum of the center of mass velocity and the rotation
        ! contribution.

        use precision_mod, only : dp
        use utils_mod    , only : cross_product

        class(eulerian_solid), intent(in) :: self
        real(dp)             , intent(in) :: X(3)
        integer              , intent(in) :: dir
        real(dp)                          :: a

        ! Local variable
        integer  :: i
        real(dp) :: d(3), r(3), omega_vec_r(3), omega_vec_omega_vec_r(3), alpha_vec_r(3)

        ! Position vector of the point X with respect to the rotation center
        d(1) = sqrt( (X(1) - self%rot_center(1))**2 + (X(2) - self%rot_center(2))**2)
        d(2) = sqrt( (X(1) + self%G%Lx - self%rot_center(1))**2 + (X(2) - self%rot_center(2))**2)
        d(3) = sqrt( (X(1) - self%G%Lx - self%rot_center(1))**2 + (X(2) - self%rot_center(2))**2)
        i = minloc(d, 1)
        if (i == 1) then
            r = X - self%rot_center
        elseif (i == 2) then
            r = X + [self%G%Lx, 0.0_dp, 0.0_dp] - self%rot_center
        elseif (i == 3) then
            r = X - [self%G%Lx, 0.0_dp, 0.0_dp] - self%rot_center
        endif

        ! Velocity due to rotation
        omega_vec_r = cross_product(self%center_of_mass%V(4:6), r)

        ! Acceleration due to this velocity
        omega_vec_omega_vec_r = cross_product(self%center_of_mass%V(4:6), omega_vec_r)

        ! Angular acceleration contribution
        alpha_vec_r = cross_product(self%center_of_mass%A(4:6), r)

        ! Final acceleration
        a = self%center_of_mass%A(dir) + omega_vec_omega_vec_r(dir) + alpha_vec_r(dir)

    end function acceleration
    !==============================================================================================

    !==============================================================================================
    subroutine check_periodicity(self)

        use precision_mod     , only : dp

        ! In/Out variables
        class(eulerian_solid), intent(inout) :: self

        ! Local variables
        integer :: n

        if (self%G%periodic_bc(1) .eqv. .true.) then
            if (self%center_of_mass%X(1) > self%G%origin(1) + self%G%Lx) then
                self%center_of_mass%X(1) = self%center_of_mass%X(1) - self%G%Lx
                do n = 1,size(self%surface_points)
                    self%surface_points(n)%X = self%surface_points(n)%X - [self%G%Lx, 0.0_dp, 0.0_dp]
                end do
            elseif (self%center_of_mass%X(1) < self%G%origin(1)) then
                self%center_of_mass%X(1) = self%center_of_mass%X(1) + self%G%Lx
                do n = 1,size(self%surface_points)
                    self%surface_points(n)%X = self%surface_points(n)%X + [self%G%Lx, 0.0_dp, 0.0_dp]
                end do
            endif
        endif

        ! Update rotation center
        self%rot_center = self%center_of_mass%X(1:3)

    end subroutine check_periodicity
    !==============================================================================================

!     !==============================================================================================
!     subroutine print_csv(self, time)

!         ! Output center of mass variables

!         ! In/Out variables
!         class(eulerian_solid), intent(in) :: self
!         real(dp)             , intent(in) :: time
  
!         write(self%file_id,12) time, &
!            ",", self%center_of_mass%X(1) , ",", self%center_of_mass%X(2) , ",", self%center_of_mass%X(3) , &
!            ",", self%center_of_mass%X(4) , ",", self%center_of_mass%X(5) , ",", self%center_of_mass%X(6) , &
!            ",", self%center_of_mass%V(1) , ",", self%center_of_mass%V(2) , ",", self%center_of_mass%V(3) , &
!            ",", self%center_of_mass%V(4) , ",", self%center_of_mass%V(5) , ",", self%center_of_mass%V(6) , &
!            ",", self%center_of_mass%A(1) , ",", self%center_of_mass%A(2) , ",", self%center_of_mass%A(2) , &
!            ",", self%center_of_mass%A(4) , ",", self%center_of_mass%A(5) , ",", self%center_of_mass%A(6) , &
!            ",", self%center_of_mass%Fh(1), ",", self%center_of_mass%Fh(2), ",", self%center_of_mass%Fh(3), &
!            ",", self%center_of_mass%Fh(4), ",", self%center_of_mass%Fh(5), ",", self%center_of_mass%Fh(6)
  
!         flush(self%file_id)
! 12 format(E16.8,24(A1,E16.8))
  
!     end subroutine print_csv
!     !==============================================================================================
  
end module 

