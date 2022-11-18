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

    implicit none

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
    Nx = 8
    Ny = 8
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

    C = circle(X = [0.5_dp, 0.5_dp, 0.0_dp], R = 0.25_dp)
    solid_list(1)%pS => C

    ! Setup eulerian fields for IBM
    call init_eulerian_ibm(solid_list)

    ! Print eulerian fields
    call temp%allocate(1)
    temp%f = ibm_index(:,:,:,1)
    call temp%write('tag_c')
    temp%f = ibm_index(:,:,:,2)
    call temp%write('tag_x')
    temp%f = ibm_index(:,:,:,3)
    call temp%write('tag_y')

    ! free memory
    call base_grid%destroy()
    call temp%destroy()
    call destroy_ibm

   ! Finalize the simulation
   call MPI_FINALIZE(ierror)

end program main