module solver_mod

    ! This module contains all the procedures for the setup of the solver.
    ! The subroutine advance_solution is a general subroutine which will point
    ! to the requested solver.

    use precision_mod, only : dp
    use grid_mod     , only : grid

    implicit none

    interface
        subroutine advance_solver(comp_grid, step, dt)
            use precision_mod, only : dp
            use grid_mod     , only : grid
            type(grid), intent(in   ) :: comp_grid
            integer   , intent(in   ) :: step
            real(dp)  , intent(inout) :: dt
        end subroutine advance_solver

        subroutine solver_status(log_id, step, time, dt)
            use precision_mod, only : dp
            integer , intent(in) :: log_id, step
            real(dp), intent(in) :: time, dt
        end subroutine solver_status
    end interface

    procedure(advance_solver), pointer :: advance_solution => Null()
    procedure(solver_status) , pointer :: print_solver_status => Null()

contains

    !===============================================================================================
    subroutine init_solver(comp_grid)

        ! Perform all preliminary operations before start the simulation

        use poisson_mod      , only : init_poisson_solver
        use navier_stokes_mod, only : allocate_navier_stokes_fields, phi, navier_stokes_solver
        use navier_stokes_mod, only : print_navier_stokes_solver_status
#if defined(MF) || defined(NN)
        use navier_stokes_mod, only : constant_viscosity
#endif
#ifdef MF
        use navier_stokes_mod  , only : rho, mu
        use multiphase_mod     , only : allocate_multiphase_fields, update_material_properties
        use multiphase_mod     , only : rho_0, rho_1, rhomin, irhomin, advect_interface
        use volume_of_fluid_mod, only : allocate_vof_fields, get_vof_from_distance, vof, advect_vof
#endif
#ifdef IBM
        use ibm_mod
#endif
#ifdef FSI
        use fsi_mod
#endif
        ! In/out variables
        type(grid), intent(in   ) :: comp_grid

        ! Allocate fields for solving Navier-Stokes equation
        call allocate_navier_stokes_fields(comp_grid)

        ! Initialize the Poisson solver
        call init_Poisson_Solver(phi)

        ! Select the solver
#ifdef FSI      
        ! Point to the weak coupling FSI solver
        advance_solution => weak_coupling_solver

        ! Point the Navier-Stokes solver print status
        print_solver_status => print_navier_stokes_solver_status

#else
        ! Point to the Navier-Stokes solver
        advance_solution => navier_stokes_solver

        ! Point the Navier-Stokes solver print status
        print_solver_status => print_navier_stokes_solver_status
#endif

#if defined(MF) || defined(NN)
        constant_viscosity = .false.
#endif

#ifdef MF
        call allocate_vof_fields(comp_grid)
        call allocate_multiphase_fields(comp_grid)
        call get_vof_from_distance
        call update_material_properties(rho, mu, vof)
        ! set the minimum density for the Dodd & Ferrante method
        rhomin = min(rho_0, rho_1)
        irhomin = 1.0_dp/rhomin
        advect_interface => advect_vof
#endif
#ifdef IBM
        call init_ibm(comp_grid)
#endif

    end subroutine init_solver
    !===============================================================================================

    !===============================================================================================
    subroutine save_fields(step)

        use scalar_mod
        use fields_mod , only : face_to_center

        ! This subroutine save the fields of the simulation at timestep step
        use navier_stokes_mod, only : v, p
#ifdef MF
        use volume_of_fluid_mod, only : vof
#endif
#ifdef FSI
        use ibm_mod, only : print_solid
#endif

        ! In/Out variables
        integer, intent(in) :: step

        ! Local variables
        type(scalar)      :: temp
        character(len=7 ) :: sn
        character(len=20) :: filename

        ! Allocate temporary scalar field
        call temp%allocate(v%G, v%x%gl)

        write(sn,'(I0.7)') step
        
        filename = 'data/vx_'//sn//'.raw'
        call face_to_center(v%x, temp, 'x')
        call temp%write(filename)

        filename = 'data/vy_'//sn//'.raw'
        call face_to_center(v%y, temp, 'y')
        call temp%write(filename)

#if DIM==3
        filename = 'data/vz_'//sn//'.raw'
        call face_to_center(v%z, temp, 'z')
        call temp%write(filename)
#endif
        filename = 'data/p_'//sn//'.raw'
        call p%write(filename)
#ifdef MF
        filename = 'data/vof_'//sn//'.raw'
        call vof%write(filename)
#endif

#ifdef FSI
        call print_solid(step)
#endif

        call temp%destroy()

    end subroutine save_fields
    !===============================================================================================

    !===============================================================================================
    subroutine save_state(step)

        ! This subroutine save the fields of the simulation at timestep step for restart purpose
        use mpi
        use decomp_2d_io, only : decomp_2d_write_var
        use navier_stokes_mod, only : v, dv_o, p

        ! In/Out variables
        integer, intent(in) :: step

        ! Local variables
        integer                       :: fh, ierror, lo(3), hi(3)
        integer(kind=MPI_OFFSET_KIND) :: filesize, disp
        character(len=7 )             :: sn
        character(len=22)             :: filename

        lo = p%G%lo
        hi = p%G%hi

        write(sn,'(I0.7)') step
        filename = 'data/state_'//sn//'.raw'
        call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_CREATE+MPI_MODE_WRONLY, &
                                MPI_INFO_NULL, fh, ierror)
        filesize = 0_MPI_OFFSET_KIND
        call MPI_FILE_SET_SIZE(fh, filesize, ierror)  ! guarantee overwriting
        disp = 0_MPI_OFFSET_KIND
        
        call decomp_2d_write_var(fh, disp, 1,      p%f(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
        call decomp_2d_write_var(fh, disp, 1,    v%x%f(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
        call decomp_2d_write_var(fh, disp, 1,    v%y%f(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
        call decomp_2d_write_var(fh, disp, 1, dv_o%x%f(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
        call decomp_2d_write_var(fh, disp, 1, dv_o%y%f(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
#if DIM==3
        call decomp_2d_write_var(fh, disp, 1,    v%z%f(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
        call decomp_2d_write_var(fh, disp, 1, dv_o%z%f(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
#endif
        call MPI_FILE_CLOSE(fh,ierror)

    end subroutine
    !===============================================================================================

    !===============================================================================================
    subroutine load_state(step)

        ! This subroutine save the fields of the simulation at timestep step for restart purpose
        use mpi
        use decomp_2d_io, only : decomp_2d_read_var
        use navier_stokes_mod, only : v, dv_o, p
        use scalar_mod

        ! In/Out variables
        integer, intent(in) :: step

        ! Local variables
        integer                       :: fh, ierror, lo(3), hi(3)
        integer(kind=MPI_OFFSET_KIND) :: filesize, disp
        character(len=7 )             :: sn
        character(len=22)             :: filename
        type(scalar)                  :: tmp
        real(dp), allocatable         :: temp(:,:,:)
        
        lo = p%G%lo
        hi = p%G%hi

        call tmp%allocate(p%G, 1)
        allocate(temp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))

        write(sn,'(I0.7)') step
        filename = 'data/state_'//sn//'.raw'
        call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierror)
        disp = 0_MPI_OFFSET_KIND
        
        call decomp_2d_read_var(fh, disp, 1,      p%f(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
        call decomp_2d_read_var(fh, disp, 1,    v%x%f(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
        call decomp_2d_read_var(fh, disp, 1,    v%y%f(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
        call decomp_2d_read_var(fh, disp, 1, dv_o%x%f(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
        call decomp_2d_read_var(fh, disp, 1, dv_o%y%f(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
#if DIM==3
        call decomp_2d_read_var(fh, disp, 1,    v%z%f(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
        call decomp_2d_read_var(fh, disp, 1, dv_o%z%f(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
#endif
        call p%update_ghost_nodes()
        call v%update_ghost_nodes()
        call MPI_FILE_CLOSE(fh,ierror)

    end subroutine load_state
    !===============================================================================================

    !===============================================================================================
    subroutine destroy_solver

        use poisson_mod      , only : destroy_poisson_solver
        use navier_stokes_mod, only : phi, destroy_navier_stokes_solver
#ifdef MF
        use volume_of_fluid_mod, only : destroy_vof
        use multiphase_mod     , only : p_o, p_hat, grad_p_hat
#endif
#ifdef ibm
        use ibm_mod
#endif

        call destroy_poisson_solver(phi)
        call destroy_navier_stokes_solver
#ifdef MF
        call destroy_vof
        call p_hat%destroy()
        call grad_p_hat%destroy()
        call p_o%destroy()
#endif
#ifdef ibm
        call destroy_ibm()
#endif

    end subroutine destroy_solver
    !===============================================================================================

end module