program main

    use mpi
    use precision_mod
    use grid_mod        
    use volume_of_fluid_mod
    use IO_mod
    
    implicit none

    integer  :: ierror, Nx, Ny, Nz
    real(dp) :: Lx, Ly, Lz, origin(3)
    type(grid)    :: comp_grid
    type(scalar)  :: phi 
    
    ! Initialize MPI
    call mpi_init(ierror)

    ! The domain is a unit squared box
    Lx = 1.0_dp
    Ly = 1.0_dp
    
    ! Set the resolution
    Nx = 16
    Ny = 16
    Nz = 1
    Lz = Lx*float(Nz)/float(Nx)
    origin = [0.0_dp, 0.0_dp, 0.0_dp]

    ! Create the grid
    call comp_grid%setup(Nx, Ny, Nz, Lx, Ly, Lz, origin, 1, 1)

    ! Allocate all vof fields
    call allocate_vof_fields(comp_grid)
    
    ! Set the initial VoF from a distance function
    distance => circle
    call get_vof_from_distance

    ! Reconstruct the interface
    call get_h_from_vof

    ! Print the hyperbolic tangent H on a file
        call h%write('h.dat')

    ! Get distance function from H
    call phi%allocate(comp_grid, 1)
    call get_distance(phi, h)

    ! Finalize the simulation
    call MPI_FINALIZE(ierror)

contains

    !========================================================================================
    function circle(x,y) result(d)

        ! In/Out variables
        real(dp), intent(in) :: x, y

        ! Local variables
        real(dp) :: x0, y0, r, d

        x0 = 0.5_dp
        y0 = 0.5_dp
        r  = 0.433_dp
        d = sqrt((x - x0)**2 + (y - y0)**2) - r
        
    end function circle
    !========================================================================================

    !========================================================================================
    subroutine get_distance(phi, h)

        type(scalar), intent(in   ) :: h
        type(scalar), intent(inout) :: phi

        integer :: i, j, k
        
        do k = comp_grid%lo(3),comp_grid%hi(3)
            do j = comp_grid%lo(2),comp_grid%hi(2)
                do i = comp_grid%lo(1),comp_grid%hi(1)
                    phi%f(i,j,k) = -1.0_dp/beta*atanh(2.0_dp*h%f(i,j,k) - 1.0_dp)*comp_grid%delta
                end do
            end do
        end do

    end subroutine get_distance
    !========================================================================================
  
end program main
