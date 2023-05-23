module fsi_mod

    implicit none

contains

    !========================================================================================
    subroutine weak_coupling_solver(comp_grid, step, dt)

        use precision_mod     , only : dp
        use grid_mod
        use ibm_mod           , only : lagrangian_solid_list
        use lagrangian_ibm_mod, only : lagrangian_compute_hydrodynamic_loads => compute_hydrodynamic_loads
        use lagrangian_ibm_mod, only : advance_structure
        !use eulerian_ibm , only : tag_cells, compute_hydrodynamic_loads
        use navier_stokes_mod, only : navier_stokes_solver, v, p, mu, rho, density, g

        ! In/Out variables
        type(grid), intent(in   ) :: comp_grid
        integer   , intent(in   ) :: step
        real(dp)  , intent(inout) :: dt

        ! Local variables
        integer :: b

        ! First solve solid object motion
        ! do b = 1,size(Eulerian_Solid_list)
        !     ! Evalute hydrodynamic loads on solid object b
        !     call compute_hydrodynamic_loads(eulerian_solid_list(b)%pS, v, p, mu, rho, g)

        !     ! Advance solid object b
        !     call eulerian_solid_list(b)%pS%advance(dt)

        !     ! Check if the solid object is outside of the domain and in case traslate it
        !     call check_periodicity(Eulerian_Solid_list(b)%pS)

        !     ! Update eularian ibm variables
        !     call tag_cells(Eulerian_Solid_list)
        ! end do

        ! If solving solid body dynamics with Lagrangian IBM
        if (size(lagrangian_solid_list) > 0) then
            ! For every solid body:
            do b = 1,size(lagrangian_solid_list)
               ! Compute hydrodynamic loads with fluid fields at step n + 1 and structure
               ! variables at step n
               call lagrangian_compute_hydrodynamic_loads(lagrangian_solid_list(b)%pS, v, p, mu, rho, g)
   
               ! Solve the structure dynamics with the compute hydrodynamic loads
#ifdef MF
               ! For multiphase flows find the density of the surrounding fluid
               ie = base_grid%closest_grid_node(lagrangian_solid_array(b)%center_of_mass(), 0)
               if (ie(2) >= base_grid%lo(2) .and. ie(2) <= base_grid%hi(2)) then
                  rhof = rho%f(ie(1),ie(2),1)
               else
                  rhof = 0.0_dp
               endif
               call mpi_allreduce(mpi_in_place,rhof,1,mpi_real8,mpi_sum,mpi_comm_world,ierror)
               call advance_structure(lagrangian_solid_list(b)%pS, step, dt, g, rhof, comp_grid)
#else
               call advance_structure(lagrangian_solid_list(b)%pS, step, dt, g, density, comp_grid)
#endif
            end do
         endif
   
        ! then solve fluid motion
        call navier_stokes_solver(comp_grid, step, dt)

    end subroutine weak_coupling_solver
    !========================================================================================

    ! !========================================================================================
    ! subroutine check_periodicity(solid)

    !     use precision           , only : dp
    !     use class_Grid          , only : base_grid
    !     use class_eulerian_solid

    !     ! In/Out variables
    !     class(eulerian_solid), intent(inout) :: solid

    !     ! Local variables
    !     integer :: n

    !     if (base_grid%periodic_bc(1) .eqv. .true.) then
    !         if (solid%X(1) > base_grid%origin(1) + base_grid%Lx) then
    !             solid%X(1) = solid%X(1) - base_grid%Lx
    !             do n = 1,size(solid%surface_points)
    !                 solid%surface_points(n)%X = solid%surface_points(n)%X - [base_grid%Lx, 0.0_dp, 0.0_dp]
    !             end do
    !         elseif (solid%X(1) < base_grid%origin(1)) then
    !             solid%X(1) = solid%X(1) + base_grid%Lx
    !             do n = 1,size(solid%surface_points)
    !                 solid%surface_points(n)%X = solid%surface_points(n)%X + [base_grid%Lx, 0.0_dp, 0.0_dp]
    !             end do
    !         endif
    !     endif

    !     ! Update rotation center
    !     solid%rot_center = solid%X(1:3)

    ! end subroutine check_periodicity
    ! !========================================================================================

end module