program main

    ! Program to test the proper creation of eulerian solids

    use mpi
    use constants            , only : pi
    use precision            , only : dp
    use class_Grid           , only : base_grid, bc_type
    use class_scalar
    use class_eulerian_circle
    use class_eulerian_solid
    use eulerian_ibm
    use json

    implicit none

    ! Parameters
    real(dp), parameter :: xc = 0.525_dp
    real(dp), parameter :: yc = 0.484_dp
    real(dp), parameter :: rc = 0.25_dp

    ! Variables
    integer       :: ierror, Nx, Ny, Nz
    real(dp)      :: Lx, Ly, Lz, origin(3)
    type(bc_type) :: bc(4)
    type(scalar)  :: temp

    ! Create one solid of type circle, it must be target
    type(circle), target :: C
    type(eulerian_solid_pointer) :: solid_list(1)

    ! Initialize MPI
    call mpi_init(ierror)

    ! The domain is a unit squared box
    Nx = 16
    Ny = 16
    Nz = 1
    Lx = 1.0_dp
    Ly = 1.0_dp
    Lz = Lx*float(Nz)/float(Nx)
    origin = [0.0_dp, 0.0_dp, 0.0_dp]
    bc(1)%s = 'Periodic'
    bc(2)%s = 'Periodic'
    bc(3)%s = 'Periodic'
    bc(4)%s = 'Periodic'
    ! Create the grid
    call base_grid%setup(Nx, Ny, Nz, Lx, Ly, Lz, origin, 1, 1, bc)

    C = circle(X = [xc, yc, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], R = rc)
    solid_list(1)%pS => C

    ! Setup eulerian fields for IBM
    call init_eulerian_ibm(solid_list)

    ! Print eulerian fields
    call temp%allocate(1)
    temp%f = ibm_index(:,:,:,0)
    call temp%write('tag_c.raw')
    temp%f = ibm_index(:,:,:,1)
    call temp%write('tag_x.raw')
    temp%f = ibm_index(:,:,:,2)
    call temp%write('tag_y.raw')

    ! Print json file for prostprocessing use
    call print_setup_json(0.0_dp)
    call case_setup

    ! free memory
    call base_grid%destroy()
    call temp%destroy()
    call destroy_ibm

   ! Finalize the simulation
   call MPI_FINALIZE(ierror)

contains 
  !=======================================================================================
   subroutine case_setup()

    integer :: json_case_id

    if (base_grid%rank == 0) then
      open(newunit = json_case_id, file = 'case.json')
      write(json_case_id,'(A1)') '{'
      write(json_case_id,'(4x,A9)') '"Case": {'
      write(json_case_id,'(8x,A6,1x,E16.8,A1)') '"xc": ', xc, ','
      write(json_case_id,'(8x,A6,1x,E16.8,A1)') '"yc": ', yc, ','
      write(json_case_id,'(8x,A6,1x,E16.8)'   ) '"rc": ', rc
      write(json_case_id,'(4x,A3)') '}'
      write(json_case_id,'(A1)') '}'
      flush(json_case_id)
      close(json_case_id)
    endif

  end subroutine case_setup
  !=======================================================================================
end program main
