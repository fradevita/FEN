module class_Vector

  ! This module contains the vector type definition and procedures.
  ! The vector field type is a type with three scalar fields.

  use precision, only : dp
  use class_Scalar
  
  implicit none

  private
  public :: vector

  type vector
     ! Three scalar fields
     type(scalar) :: x, y, z

   contains
     procedure :: allocate
     procedure :: apply_bc
     procedure :: destroy
  end type vector

contains

  !========================================================================================
  subroutine allocate(v, l)

    ! Procedure to allocate memory for the vector v,
    ! l is the number of ghost node per side.

    ! In/Out variables
    class(vector), intent(inout)           :: v
    integer      , intent(in   ), optional :: l

    ! By default the field has zero ghost nodes, if l is present set the number of
    ! ghost nodes per side equal to l
    if (present(l)) then
       v%x%gl = l
       v%y%gl = l
       v%z%gl = l
    endif

    ! Allocate memory
    call v%x%allocate(v%x%gl)
    call v%y%allocate(v%y%gl)
#if DIM == 3
    call v%z%allocate(v%z%gl)
#endif

  end subroutine allocate
  !========================================================================================

  !========================================================================================
  subroutine destroy(v)

    ! Free all memory allocated by vector field v.

    ! In/Out variables
    class(vector), intent(inout) :: v

    call v%x%destroy()
    call v%y%destroy()
#if DIM == 3
    call v%z%destroy()
#endif

  end subroutine destroy
  !========================================================================================

  !========================================================================================
  subroutine apply_bc(v)

    ! Apply boundary condition on all three scrala field of v.
    
    use mpi
    use io        , only : print_error_message
    use class_Grid, only : base_grid
    use halo      , only : field_halo_update

    ! In/Out variables
    class(vector), intent(inout) :: v

    ! Local variables
    integer  :: ierror, lo(3), hi(3), l

    lo = base_grid%lo
    hi = base_grid%hi

    ! Apply bc on each component starting from v%x

    ! **** V%x ****
    ! First update halos
    call field_halo_update(v%x%f, v%x%gl)

    ! Left Boundary
    if (v%x%bc%type_l == 0) then ! Periodic
       do l = 1,v%x%gl
          v%x%f(lo(1)-l,:,:) = v%x%f(hi(1)-l+1,:,:)
       end do
    elseif (v%x%bc%type_l == 1) then ! Dirichlet
       v%x%f(base_grid%lo(1)-1,:,:) = v%x%bc%l(:,:)
    elseif (v%x%bc%type_l == 2) then ! Neumann
       v%x%f(base_grid%lo(1)-1,:,:) = v%x%f(base_grid%lo(1),:,:)
    else
       call print_error_message('ERROR: wrong left boundary condition type for scalar s')
    endif

    ! Right Boundary
    if (v%x%bc%type_r == 0) then ! Periodic
       do l = 1,v%x%gl
          v%x%f(hi(1)+l,:,:) = v%x%f(lo(1)+l-1,:,:)
       end do
    elseif (v%x%bc%type_r == 1) then ! Dirichlet
       v%x%f(base_grid%hi(1)  ,:,:) = v%x%bc%r(:,:)
       v%x%f(base_grid%hi(1)+1,:,:) = v%x%bc%r(:,:)
    elseif (v%x%bc%type_r == 2) then ! Neumann
       v%x%f(base_grid%hi(1)+1,:,:) = v%x%f(base_grid%hi(1),:,:)
    else
       call print_error_message('ERROR: wrong rigth boundary condition type for scalar s')
    endif

    ! Bottom Boundary
    !if (base_grid%rank == 0) then
    if (base_grid%lo(2) == 1) then
       if (v%x%bc%type_b == 0) then ! Periodic
          ! Need to update this bc only if using 1 rank in y
          if (base_grid%prow == 1) then
             do l = 1,v%x%gl
                v%x%f(:,lo(2)-1,:) = v%x%f(:,hi(2),:)
             end do
          endif
       elseif (v%x%bc%type_b == 1) then ! Dirichlet
          v%x%f(:,base_grid%lo(2)-1,:) = 2.0_dp*v%x%bc%b(:,:) - v%x%f(:,base_grid%lo(2),:)
       elseif (v%x%bc%type_b == 2) then ! Neumann
          v%x%f(:,base_grid%lo(2)-1,:) = v%x%f(:,base_grid%lo(2),:)
       else
          call print_error_message('ERROR: wrong bottom boundary condition type for scalar s')
       endif
    endif

    ! Top Boundary
    !if (base_grid%rank == base_grid%nranks - 1) then
    if (base_grid%hi(2) == base_grid%Ny) then
       if (v%x%bc%type_t == 0) then ! Periodic
          ! Need to update this bc only if using 1 rank in y
          if (base_grid%prow == 1) then
             do l = 1,v%x%gl
                v%x%f(:,hi(2)+1,:) = v%x%f(:,lo(2),:)
             end do
          endif
       elseif (v%x%bc%type_t == 1) then ! Dirichlet
          v%x%f(:,base_grid%hi(2)+1,:) = 2.0_dp*v%x%bc%t(:,:) - v%x%f(:,base_grid%hi(2),:)
       elseif (v%x%bc%type_t == 2) then ! Neumann
          v%x%f(:,base_grid%hi(2)+1,:) = v%x%f(:,base_grid%hi(2),:)
       else
          call print_error_message('ERROR: wrong top boundary condition type for scalar s')
       endif
    endif

#if DIM==3
    if (base_grid%pcol == 1) then
       v%x%f(:,:,lo(3)-1) = v%x%f(:,:,hi(3))
       v%x%f(:,:,hi(3)+1) = v%x%f(:,:,lo(3))
    end if
#endif
    
    ! This is necessary to avoid deadlock
    call mpi_barrier(mpi_comm_world, ierror)

    ! **** V%y ****
    ! First update halos
    call field_halo_update(v%y%f, v%y%gl)

    ! Left Boundary
    if (v%y%bc%type_l == 0) then ! Periodic
       do l = 1,v%y%gl
          v%y%f(lo(1)-l,:,:) = v%y%f(hi(1)-l+1,:,:)
       end do
    elseif (v%y%bc%type_l == 1) then ! Dirichlet
       v%y%f(base_grid%lo(1)-1,:,:) = 2.0_dp*v%y%bc%l(:,:) - v%y%f(base_grid%lo(1),:,:)
    elseif (v%y%bc%type_l == 2) then ! Neumann
       v%y%f(base_grid%lo(1)-1,:,:) = v%y%f(base_grid%lo(1),:,:)
    else
       call print_error_message('ERROR: wrong left boundary condition type for scalar s')
    endif

    ! Right Boundary
    if (v%y%bc%type_r == 0) then ! Periodic
       do l = 1,v%y%gl
          v%y%f(hi(1)+l,:,:) = v%y%f(lo(1)+l-1,:,:)
       end do
    elseif (v%y%bc%type_r == 1) then ! Dirichlet
       v%y%f(base_grid%hi(1)+1,:,:) = 2.0_dp*v%y%bc%r(:,:) - v%y%f(base_grid%hi(1),:,:) 
    elseif (v%y%bc%type_r == 2) then ! Neumann
       v%y%f(base_grid%hi(1)+1,:,:) = v%y%f(base_grid%hi(1),:,:)
    else
       call print_error_message('ERROR: wrong rigth boundary condition type for scalar s')
    endif

    ! Bottom Boundary
    !if (base_grid%rank == 0) then
    if (base_grid%lo(2) == 1) then
       if (v%y%bc%type_b == 0) then ! Periodic
          ! Need to update this bc only if using 1 rank in y
          if (base_grid%prow == 1) then
             do l = 1,v%y%gl
                v%y%f(:,lo(2)-1,:) = v%y%f(:,hi(2),:)
             end do
          endif
       elseif (v%y%bc%type_b == 1) then ! Dirichlet
          v%y%f(:,base_grid%lo(2)-1,:) = v%y%bc%b(:,:)
       elseif (v%y%bc%type_b == 2) then ! Neumann
          v%y%f(:,base_grid%lo(2)-1,:) = v%y%f(:,base_grid%lo(2),:)
       else
          call print_error_message('ERROR: wrong bottom boundary condition type for scalar s')
       endif
    endif

    ! Top Boundary
    !if (base_grid%rank == base_grid%nranks - 1) then
    if (base_grid%hi(2) == base_grid%Ny) then
       if (v%y%bc%type_t == 0) then ! Periodic
          ! Need to update this bc only if using 1 rank in y
          if (base_grid%prow == 1) then
             do l = 1,v%y%gl
                v%y%f(:,hi(2)+1,:) = v%y%f(:,lo(2),:)
             end do
          endif
       elseif (v%y%bc%type_t == 1) then ! Dirichlet
          v%y%f(:,base_grid%hi(2)  ,:) = v%y%bc%t(:,:)
          v%y%f(:,base_grid%hi(2)+1,:) = v%y%bc%t(:,:)
       elseif (v%y%bc%type_t == 2) then ! Neumann
          v%y%f(:,base_grid%hi(2)+1,:) = v%y%f(:,base_grid%hi(2),:)
       else
          call print_error_message('ERROR: wrong top boundary condition type for scalar s')
       endif
    endif

#if DIM==3
    if (base_grid%pcol == 1) then
       v%y%f(:,:,lo(3)-1) = v%y%f(:,:,hi(3))
       v%y%f(:,:,hi(3)+1) = v%y%f(:,:,lo(3))
    end if
#endif
    
    ! This is necessary to avoid deadlock
    call mpi_barrier(mpi_comm_world, ierror)

#if DIM==3
    ! **** V%z ****
    ! First update halos
    call field_halo_update(v%z%f, v%z%gl)

    ! Left Boundary
    if (v%z%bc%type_l == 0) then ! Periodic
       do l = 1,v%z%gl
          v%z%f(lo(1)-l,:,:) = v%z%f(hi(1)-l+1,:,:)
       end do
    elseif (v%z%bc%type_l == 1) then ! Dirichlet
       v%z%f(base_grid%lo(1)-1,:,:) = v%z%bc%l(:,:)
    elseif (v%z%bc%type_l == 2) then ! Neumann
       v%z%f(base_grid%lo(1)-1,:,:) = v%z%f(base_grid%lo(1),:,:)
    else
       call print_error_message('ERROR: wrong left boundary condition type for scalar s')
    endif

    ! Right Boundary
    if (v%z%bc%type_r == 0) then ! Periodic
       do l = 1,v%z%gl
          v%z%f(hi(1)+l,:,:) = v%z%f(lo(1)+l-1,:,:)
       end do
    elseif (v%z%bc%type_r == 1) then ! Dirichlet
       v%z%f(base_grid%hi(1)  ,:,:) = v%z%bc%r(:,:)
       v%z%f(base_grid%hi(1)+1,:,:) = v%z%bc%r(:,:)
    elseif (v%z%bc%type_r == 2) then ! Neumann
       v%z%f(base_grid%hi(1)+1,:,:) = v%z%f(base_grid%hi(1),:,:)
    else
       call print_error_message('ERROR: wrong rigth boundary condition type for scalar s')
    endif

    ! Bottom Boundary
    !if (base_grid%rank == 0) then
    if (base_grid%lo(2) == 1) then
       if (v%z%bc%type_b == 0) then ! Periodic
          ! Need to update this bc only if using 1 rank in y
          if (base_grid%prow == 1) then
             do l = 1,v%z%gl
                v%z%f(:,lo(2)-1,:) = v%z%f(:,hi(2),:)
             end do
          endif
       elseif (v%z%bc%type_b == 1) then ! Dirichlet
          v%z%f(:,base_grid%lo(2)-1,:) = 2.0_dp*v%z%bc%b(:,:) - v%z%f(:,base_grid%lo(2),:)
       elseif (v%z%bc%type_b == 2) then ! Neumann
          v%z%f(:,base_grid%lo(2)-1,:) = v%z%f(:,base_grid%lo(2),:)
       else
          call print_error_message('ERROR: wrong bottom boundary condition type for scalar s')
       endif
    endif

    ! Top Boundary
    !if (base_grid%rank == base_grid%nranks - 1) then
    if (base_grid%hi(2) == base_grid%Ny) then
       if (v%z%bc%type_t == 0) then ! Periodic
          ! Need to update this bc only if using 1 rank in y
          if (base_grid%prow == 1) then
             do l = 1,v%z%gl
                v%z%f(:,hi(2)+1,:) = v%z%f(:,lo(2),:)
             end do
          endif
       elseif (v%z%bc%type_t == 1) then ! Dirichlet
          v%z%f(:,base_grid%hi(2)+1,:) = 2.0_dp*v%z%bc%t(:,:) - v%z%f(:,base_grid%hi(2),:)
       elseif (v%z%bc%type_t == 2) then ! Neumann
          v%z%f(:,base_grid%hi(2)+1,:) = v%z%f(:,base_grid%hi(2),:)
       else
          call print_error_message('ERROR: wrong top boundary condition type for scalar s')
       endif
    endif

    if (base_grid%pcol == 1) then
       v%z%f(:,:,lo(3)-1) = v%z%f(:,:,hi(3))
       v%z%f(:,:,hi(3)+1) = v%z%f(:,:,lo(3))
    end if
    
    ! This is necessary to avoid deadlock
    call mpi_barrier(mpi_comm_world, ierror)
#endif
    
  end subroutine apply_bc
  !========================================================================================

end module class_Vector
