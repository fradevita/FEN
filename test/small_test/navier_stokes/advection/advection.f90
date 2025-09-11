program test_advection

    use mpi
    use IO_mod
    use precision_mod     , only : dp
    use global_mod        , only : pi
    use grid_mod
    use vector_mod
    use solver_mod        , only : init_solver, destroy_solver
    use navier_stokes_mod , only : v, add_advection

    implicit none

    integer          :: ierror, Nx, Ny, Nz, n, i, j, k
    real(dp)         :: Lx, Ly, Lz, origin(3), x, y, error, Linf_x, Linf_y
    type(grid)       :: comp_grid
    type(vector)     :: Adv
    type(bc_type)    :: bc(4)

    ! Initialize MPI
    call mpi_init(ierror)

    ! The domain is a periodic squared box of size 2 pi
    Lx = 2.0_dp*pi
    Ly = Lx
    origin = [0.0_dp, 0.0_dp, 0.0_dp]
    bc(1)%s = 'Periodic'
    bc(2)%s = 'Periodic'
    bc(3)%s = 'Periodic'
    bc(4)%s = 'Periodic'

    ! Solve the problem for several different levels resolution
    refinement_loop: do n = 3,8

        ! Setup the grid
        Nx = 2**n
        Ny = Nx
        Nz = 1
        Lz = Lx*float(Nz)/float(Nx)
        call comp_grid%setup(Nx, Ny, Nz, Lx, Ly, Lz, origin, 1, 1, bc)

        ! Initialize the solver
        call init_solver(comp_grid)

        ! Set initial fields
        call init_fields

        ! Compute advective terms
        call Adv%allocate(comp_grid)
        call add_advection(comp_grid, Adv)

        ! Evaluate error
        Linf_x = 0.0_dp
        Linf_y = 0.0_dp
        do k = comp_grid%lo(3),comp_grid%hi(3)
            do j = comp_grid%lo(3),comp_grid%hi(3)
                do i = comp_grid%lo(3),comp_grid%hi(3)
                    x = float(i)*comp_grid%delta
                    y = (float(j) - 0.5_dp)*comp_grid%delta
                    error = abs(Adv%x%f(i,j,k) + Advu(x,y))
                    if (error > Linf_x) Linf_x = error
                    x = (float(i) - 0.5_dp)*comp_grid%delta
                    y = float(j)*comp_grid%delta
                    error = abs(Adv%y%f(i,j,k) + Advv(x,y))
                    if (error > Linf_y) Linf_y = error
                end do
            end do
        end do

        print *, Nx, Linf_x, Linf_y


        ! Free the memory allocated by the solver
        call destroy_solver()
        call comp_grid%destroy
        call Adv%destroy() 
        
    end do refinement_loop

    ! Finalize the simulation
    call MPI_FINALIZE(ierror)

contains

    !=======================================================================================
    real(dp) function u_x(x,y)

        real(dp) :: x, y
        u_x = sin(x)*cos(y)

    end function
    !=======================================================================================

    !=======================================================================================
    real(dp) function u_y(x,y)

        real(dp) :: x, y
        u_y = -cos(x)*sin(y)

    end function
    !=======================================================================================

    !=======================================================================================
    real(dp) function Advu(x,y)

        real(dp) :: x, y
        real(dp) :: dudx, dudy

        dudx =  cos(x)*cos(y)
        dudy = -sin(x)*sin(y)

        Advu = u_x(x,y)*dudx + u_y(x,y)*dudy

    end function
    !=======================================================================================

    !=======================================================================================
    real(dp) function Advv(x,y)

        real(dp) :: x, y
        real(dp) :: dvdx, dvdy

        dvdx =  sin(x)*sin(y)
        dvdy = -cos(x)*cos(y)

        Advv = u_x(x,y)*dvdx + u_y(x,y)*dvdy

    end function
    !=======================================================================================

    !=======================================================================================
    subroutine init_fields

        integer :: i, j, k
        real(dp) :: x, y

        do k = comp_grid%lo(3),comp_grid%hi(3)
            do j = comp_grid%lo(2),comp_grid%hi(2)
                do i = comp_grid%lo(1),comp_grid%hi(1)
                    x = float(i)*comp_grid%delta
                    y = (float(j) - 0.5_dp)*comp_grid%delta
                    v%x%f(i,j,k) = u_x(x,y)
                    x = (float(i) - 0.5_dp)*comp_grid%delta
                    y = float(j)*comp_grid%delta
                    v%y%f(i,j,k) = u_y(x,y)
                end do
            end do
        end do

        ! Update halos and bc
        call v%update_ghost_nodes

    end subroutine init_fields
    !=======================================================================================

end program test_advection
