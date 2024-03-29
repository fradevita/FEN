!> This module contains some I/O functions using iso_fortran_env.
module IO_mod

    use, intrinsic :: iso_fortran_env, only: stdout => OUTPUT_UNIT
    use, intrinsic :: iso_fortran_env, only: stderr => ERROR_UNIT

    implicit none

contains

    !==============================================================================================
    subroutine print_error_message(message)
        !< Subroutine to print an errore message on stderr
        
        ! In/Out variables
        character(len=*), intent(in)           :: message !< Input message

        write(stderr, '(A)') message

    end subroutine print_error_message
    !==============================================================================================

end module
