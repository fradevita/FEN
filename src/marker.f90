!< Module to define the marker derived type.
module marker_mod

    use precision_mod, only : dp

    implicit none

    private
    public :: marker

    type marker
        !< Base type to define a point
        real(dp)                            :: m = 1.0_dp !< Discrete mass
        real(dp), dimension(:), allocatable :: X          !< Position vector
        real(dp), dimension(:), allocatable :: V          !< Velocity vector
        real(dp), dimension(:), allocatable :: A          !< Acceleration vector
        real(dp), dimension(:), allocatable :: Fv         !< Viscous forces vector
        real(dp), dimension(:), allocatable :: Fp         !< pressure forces vector
        real(dp), dimension(:), allocatable :: Fh         !< Hydrodynamic forces vector
        real(dp), dimension(:), allocatable :: Fe         !< external forces vector
    contains
        procedure, pass(self) :: create  !< self exp.
        procedure, pass(self) :: destroy !< self exp.
    end type marker

contains

    !==============================================================================================
    subroutine create(self, degree_of_freedom)
        !< Allocate memory and set to zero all variables

        class(marker), intent(inout) :: self              !< marker object
        integer      , intent(in   ) :: degree_of_freedom !< number of degree of freedom of the marker

        ! Allocate memory
        allocate(self%X(degree_of_freedom))
        allocate(self%V(degree_of_freedom))
        allocate(self%A(degree_of_freedom))
        allocate(self%Fv(degree_of_freedom))
        allocate(self%Fp(degree_of_freedom))
        allocate(self%Fh(degree_of_freedom))
        allocate(self%Fe(degree_of_freedom))

        ! Initialize to zero
        self%X = 0.0_dp
        self%V = 0.0_dp
        self%A = 0.0_dp
        self%Fv = 0.0_dp
        self%Fp = 0.0_dp
        self%Fh = 0.0_dp
        self%Fe = 0.0_dp

    end subroutine
    !==============================================================================================

    !==============================================================================================
    subroutine destroy(self)
        ! Free the memory of the marker

        class(marker), intent(inout) :: self !< marker object

        deallocate(self%X, self%V, self%A, self%Fv, self%Fp, self%Fh, self%Fe)

    end subroutine
    !==============================================================================================

end module