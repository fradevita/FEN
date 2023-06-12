module solid_mod

    use precision_mod
    use marker_mod

    implicit none

    private
    public :: solid, solid_pointer

    !< Base abstract object for a solid body
    type, abstract :: solid
        real(dp)          :: density = 1.0_dp !< density of the solid
        type(marker)      :: center_of_mass   !< center of mass of the solid
        integer           :: file_id = -1     !< index of the output file
        character(len=99) :: name = 'unset'   !< solid name
    end type

    type solid_pointer
        class(solid), pointer :: pS => Null() !< pointer to solid
    end type solid_pointer

end module