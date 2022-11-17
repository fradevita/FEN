module multiphase

  ! This module contains the fields and the procedures to solve Navier-Stokes equations
  ! for a multiphase fluid-fluid flow.

  use precision, only : dp
  use class_Scalar
  use class_Vector

  implicit none

  ! Material properties are now scalar fields
  type(scalar) :: rho, mu, rho_o, mu_o

  ! For the one-fluid formulation the two fluids have the material properties
  ! rho0, mu0, rho1, mu1
  real(dp) :: rho_0 = 1.0_dp, rho_1 = 1.0_dp, mu_0 = 1.0_dp, mu_1 = 1.0_dp

  ! Surface tension coefficient
  real(dp) :: sigma = 0.0_dp

  ! To use the pressure approximation of Dodd and Ferrante 2014 JCP it is necessary
  ! to introduce two more fields for the pressure, one for the approximated pressure
  ! phat, and one for the old pressure, p_o
  type(scalar) :: p_hat, p_o
  type(vector) :: grad_p_hat
  
  ! And also the minimum density
  real(dp) :: rhomin, irhomin

  interface
     subroutine interface_advection(v, dt)
       use precision, only : dp
       use class_Vector
       real(dp)    , intent(in) :: dt
       type(vector), intent(in) :: v
     end subroutine interface_advection
  end interface

  procedure(interface_advection), pointer :: advect_interface => Null()

contains

  !========================================================================================
  subroutine allocate_multiphase_fields

    use io        , only : print_error_message
    use class_Grid, only : base_grid
    
    ! This subroutine allocate the memory for the fields requested by the multiphase solver
    
    call p_hat%allocate(1)
    call p_o%allocate(1)
    call grad_p_hat%allocate()

    ! Based on the boundary conditions on the physical domain select boundary conditions
    ! for all fields
    ! Left boundary
    if (base_grid%boundary_conditions(1)%s == 'Periodic') then
       p_hat%bc%type_l = 0
       p_o%bc%type_l = 0
    elseif (base_grid%boundary_conditions(1)%s == 'Wall') then
       p_hat%bc%type_l = 2
       p_o%bc%type_l = 2
    else
       call print_error_message('ERROR: wrong bc on left boundary')
    endif

    ! Right boundary
    if (base_grid%boundary_conditions(2)%s == 'Periodic') then
       p_hat%bc%type_r = 0
       p_o%bc%type_r = 0
    elseif (base_grid%boundary_conditions(2)%s == 'Wall') then
       p_hat%bc%type_r = 2
       p_o%bc%type_r = 2
    else
       call print_error_message('ERROR: wrong bc on left boundary')
    endif

        ! Bottom boundary
    if (base_grid%boundary_conditions(3)%s == 'Periodic') then
       p_hat%bc%type_b = 0
       p_o%bc%type_b = 0
    elseif (base_grid%boundary_conditions(3)%s == 'Wall') then
       p_hat%bc%type_b = 2
       p_o%bc%type_b = 2
    else
       call print_error_message('ERROR: wrong bc on left boundary')
    endif

    ! Top boundary
    if (base_grid%boundary_conditions(4)%s == 'Periodic') then
       p_hat%bc%type_t = 0
       p_o%bc%type_t = 0
    elseif (base_grid%boundary_conditions(4)%s == 'Wall') then
       p_hat%bc%type_t = 2
       p_o%bc%type_t = 2
    else
       call print_error_message('ERROR: wrong bc on left boundary')
    endif

  end subroutine allocate_multiphase_fields
  !========================================================================================
  
  !========================================================================================
  subroutine update_material_properties(rho, mu, vof)

    ! This subroutine update the material properties rho and mu of the two fluids
    ! based on the vof function

    ! In/Out variables
    type(scalar), intent(in   ) :: vof
    type(scalar), intent(inout) :: rho, mu
    
    ! Update density and viscosity with new VoF location
    rho%f = rho_1*vof%f + rho_0*(1.0_dp - vof%f)
    mu%f  =  mu_1*vof%f +  mu_0*(1.0_dp - vof%f)
    
    ! Update boundary conditions
    call rho%apply_bc()
    call mu%apply_bc()
    
  end subroutine update_material_properties
  !========================================================================================

end module multiphase
