module solid_mod

    ! This module define the abstract solid object. It will be extended either in Eulerian or 
    ! Lagrangian formulation.

    use precision_mod, only : dp
    use marker_mod   , only : marker
    use global_mod   , only : dofs

    implicit none

    private
    public :: solid, solid_pointer

    !< Base abstract object for a solid body
    type, abstract :: solid
        real(dp)              :: density = 1.0_dp !< density of the solid
        real(dp)              :: vol = 1.0_dp     !< volume
        real(dp)              :: mass = 1.0_dp    !< mass
        real(dp), allocatable :: IM(:)            !< intertial mass vector
        type(marker)          :: center_of_mass   !< center of mass of the solid
        integer               :: file_id = -1     !< index of the output file
        character(len=99)     :: name = 'unset'   !< solid name
        logical               :: fsi = .true.     !< flag to activate fsi
    contains
        ! List of object-independent procedure
        procedure, pass(self) :: write_csv
    end type

    ! The solid pointer class is used to create a list of pointers to solid objects
    ! used in the simulation.
    type solid_pointer
        class(solid), pointer :: pS => Null() !< pointer to solid object
    end type solid_pointer

contains

    !===============================================================================================
    subroutine write_csv(self, time)

        ! Objective: output on file_id center of mass variables at given time. 

        use IO_mod

        ! In/Out variables
        class(solid), intent(in) :: self
        real(dp)    , intent(in) :: time
        
        ! Local variables
        logical :: is_opened

        inquire(unit = self%file_id, opened = is_opened)
        if (is_opened) then
            write(self%file_id,'(*(E16.8,:,","))') time, self%center_of_mass%X , &
                                                         self%center_of_mass%V,  &
                                                         self%center_of_mass%A , &
                                                         self%center_of_mass%Fv, &
                                                         self%center_of_mass%Fp, &
                                                         self%center_of_mass%Fh, &
                                                         self%center_of_mass%Fe
        else  
            call print_error_message("ERROR: output file for solid "//self%name//" must be opened.")
        endif

    end subroutine write_csv
    !===============================================================================================

end module
