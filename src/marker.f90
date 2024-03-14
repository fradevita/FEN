!< Module to define the marker class.
module marker_mod

    use precision_mod, only : dp

    implicit none

    private
    public :: marker, markersAreEqual

    type marker
        integer                             :: index      !< marker index
        real(dp)                            :: m = 1.0_dp !< Discrete mass (default 1)
        real(dp), dimension(:), allocatable :: x          !< Position vector
        real(dp), dimension(:), allocatable :: v          !< Velocity vector
        real(dp), dimension(:), allocatable :: a          !< Acceleration vector
        real(dp), dimension(:), allocatable :: Fv         !< Viscous forces vector
        real(dp), dimension(:), allocatable :: Fp         !< pressure forces vector
        real(dp), dimension(:), allocatable :: Fh         !< Hydrodynamic forces vector (Fv + Fp)
        real(dp), dimension(:), allocatable :: Fe         !< External forces vector
        real(dp), dimension(:), allocatable :: Fi         !< Internal forces vector
        real(dp), dimension(:), allocatable :: Fs         !< Source forces vector
    contains
        procedure, pass(self) :: create
        procedure, pass(self) :: destroy
    end type marker

    interface marker
        procedure constructor
    end interface

contains

    !===============================================================================================
    function constructor(degree_of_freedom, index, x)

        ! Objective: create a marker. It allocates memory and initialize to zero all
        ! arrays. Optionally, set initial value of x.

        ! In/Out variables
        integer , intent(in)           :: degree_of_freedom
        integer , intent(in)           :: index
        real(dp), intent(in), optional :: x(3)              !< marker position vector
        type(marker)                   :: constructor       !< output marker object

        ! Allocate memory
        allocate(constructor%x(degree_of_freedom))
        allocate(constructor%v(degree_of_freedom))
        allocate(constructor%a(degree_of_freedom))
        allocate(constructor%Fv(degree_of_freedom))
        allocate(constructor%Fp(degree_of_freedom))
        allocate(constructor%Fh(degree_of_freedom))
        allocate(constructor%Fe(degree_of_freedom))
        allocate(constructor%Fi(degree_of_freedom))
        allocate(constructor%Fs(degree_of_freedom))

        ! Initialize to zero
        constructor%x = 0.0_dp
        constructor%v = 0.0_dp
        constructor%a = 0.0_dp
        constructor%Fv = 0.0_dp
        constructor%Fp = 0.0_dp
        constructor%Fh = 0.0_dp
        constructor%Fe = 0.0_dp
        constructor%Fi = 0.0_dp
        constructor%Fs = 0.0_dp

        ! Set index and position
        constructor%index = index
        if (present(X)) constructor%X = X

    end function
    !===============================================================================================

    !===============================================================================================
    subroutine create(self, degree_of_freedom)

        !< Allocate memory and set to zero all variables

        class(marker), intent(inout) :: self              !< marker object
        integer      , intent(in   ) :: degree_of_freedom !< degree of freedom = array size

        ! Allocate memory
        allocate(self%X(degree_of_freedom))
        allocate(self%V(degree_of_freedom))
        allocate(self%A(degree_of_freedom))
        allocate(self%Fv(degree_of_freedom))
        allocate(self%Fp(degree_of_freedom))
        allocate(self%Fh(degree_of_freedom))
        allocate(self%Fe(degree_of_freedom))
        allocate(self%Fi(degree_of_freedom))
        allocate(self%Fs(degree_of_freedom))

        ! Initialize to zero
        self%X = 0.0_dp
        self%V = 0.0_dp
        self%A = 0.0_dp
        self%Fv = 0.0_dp
        self%Fp = 0.0_dp
        self%Fh = 0.0_dp
        self%Fe = 0.0_dp
        self%Fi = 0.0_dp
        self%Fs = 0.0_dp

    end subroutine
    !===============================================================================================

    !===============================================================================================
    logical function markersAreEqual(m1, m2)

        ! Objective: check if two points are equal when reading a set of data points.

        ! In/Out variables
        type(marker), intent(in) :: m1, m2

        if (m1%X(1) == m2%X(1) .and.   &
            m1%X(2) == m2%X(2) .and.   &
            m1%X(3) == m2%X(3)) then
            markersAreEqual = .true.
        else
            markersAreEqual = .false.
        endif

    end function
    !===============================================================================================

    !===============================================================================================
    subroutine destroy(self)

        ! Objective: Free the memory of the marker.

        class(marker), intent(inout) :: self !< marker object

        deallocate(self%x, self%v, self%a, self%Fv, self%Fp, self%Fh, self%Fe, self%Fi, self%Fs)

    end subroutine
    !===============================================================================================

end module