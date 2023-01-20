module fsi

    implicit none

contains

    !========================================================================================
    subroutine weak_coupling_solver(step, dt)

        use precision    , only : dp
        use ibm          , only : Eulerian_Solid_list
        use eulerian_ibm , only : tag_cells, compute_hydrodynamic_loads
        use navier_stokes

        ! In/Out variables
        integer , intent(in   ) :: step
        real(dp), intent(inout) :: dt

        ! Local variables
        integer :: b

        ! First solve solid object motion
        do b = 1,size(Eulerian_Solid_list)
            ! Evalute hydrodynamic loads on solid object b
            call compute_hydrodynamic_loads(eulerian_solid_list(b)%pS, v, p, mu, rho, g)

            ! Advance solid object b
            call eulerian_solid_list(b)%pS%advance(dt)

            ! Check if the solid object is outside of the domain and in case traslate it
            call check_periodicity(Eulerian_Solid_list(b)%pS)

            ! Update eularian ibm variables
            call tag_cells(Eulerian_Solid_list)
        end do

        ! then solve fluid motion
        call navier_stokes_solver(step, dt)

    end subroutine weak_coupling_solver
    !========================================================================================

    !========================================================================================
    subroutine check_periodicity(solid)

        use precision           , only : dp
        use class_Grid          , only : base_grid
        use class_eulerian_solid

        ! In/Out variables
        class(eulerian_solid), intent(inout) :: solid

        ! Local variables
        integer :: n

        if (base_grid%periodic_bc(1) .eqv. .true.) then
            if (solid%X(1) > base_grid%origin(1) + base_grid%Lx) then
                solid%X(1) = solid%X(1) - base_grid%Lx
                do n = 1,size(solid%surface_points)
                    solid%surface_points(n)%X = solid%surface_points(n)%X - [base_grid%Lx, 0.0_dp, 0.0_dp]
                end do
            elseif (solid%X(1) < base_grid%origin(1)) then
                solid%X(1) = solid%X(1) + base_grid%Lx
                do n = 1,size(solid%surface_points)
                    solid%surface_points(n)%X = solid%surface_points(n)%X + [base_grid%Lx, 0.0_dp, 0.0_dp]
                end do
            endif
        endif

        ! Update rotation center
        solid%rot_center = solid%X(1:3)

    end subroutine check_periodicity
    !========================================================================================

end module fsi