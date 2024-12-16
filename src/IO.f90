!> This module contains some I/O functions using iso_fortran_env.
module IO_mod

    use, intrinsic :: iso_fortran_env, only: stdout => OUTPUT_UNIT
    use, intrinsic :: iso_fortran_env, only: stderr => ERROR_UNIT

    implicit none

contains

    !===============================================================================================
    subroutine print_error_message(message)
        !< Subroutine to print an errore message on stderr

#ifdef MPI
        use global_mod, only : myrank
#endif

        ! In/Out variables
        character(len=*), intent(in) :: message !< Input message

#ifdef MPI
        if (myrank == 0) then
#endif
        write(stderr, '(A)') message
#ifdef MPI
        endif
#endif

    end subroutine print_error_message
    !===============================================================================================

end module
