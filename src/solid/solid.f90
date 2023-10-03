module solid_mod

    use precision_mod, only : dp
    use marker_mod   , only : marker

    implicit none

    private
    public :: solid, solid_pointer

    !< Base abstract object for a solid body
    type, abstract :: solid
        real(dp)          :: density = 1.0_dp !< density of the solid
        type(marker)      :: center_of_mass   !< center of mass of the solid
        integer           :: file_id = -1     !< index of the output file
        character(len=99) :: name = 'unset'   !< solid name
    contains
        ! List of object-independent procedure
        procedure, pass(self) :: write_csv 
    end type

    type solid_pointer
        class(solid), pointer :: pS => Null() !< pointer to solid
    end type solid_pointer

contains

    !==============================================================================================
    subroutine write_csv(self, time)

        class(solid), intent(in) :: self
        real(dp)    , intent(in) :: time

        write(self%file_id,'(*(E16.8,:,","))') time, self%center_of_mass%X , self%center_of_mass%V,  &
                                                     self%center_of_mass%A , self%center_of_mass%Fh

    end subroutine write_csv
    !==============================================================================================

end module