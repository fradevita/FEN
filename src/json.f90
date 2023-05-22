module json

  implicit none

contains

  !========================================================================================
  subroutine print_setup_json(dt)

    use precision_mod , only : dp
    use grid_mod, only : base_grid

    ! In/Out variables
    real(dp), intent(in) :: dt
    
    ! Local variables
    integer :: json_id
    
    if (base_grid%rank == 0) then
      open(newunit = json_id, file = 'setup.json')
      write(json_id,'(A1)') '{'
      write(json_id,'(4x,A9)') '"Grid": {'
      write(json_id,'(8x,A6,1x,I7,A1)') '"Nx": ', base_grid%Nx, ','
      write(json_id,'(8x,A6,1x,I7,A1)') '"Ny": ', base_grid%Ny, ','
      write(json_id,'(8x,A6,1x,I7,A1)') '"Nz": ', base_grid%Nz, ','
      write(json_id,'(8x,A11,1x,E16.8,A1,E16.8,A1,E16.8,A2)') '"origin": [', &
           base_grid%origin(1), ',', base_grid%origin(2), ',', base_grid%origin(3), '],'
      write(json_id,'(8x,A6,1x,E16.8,A1)') '"Lx": ', base_grid%Lx, ','
      write(json_id,'(8x,A6,1x,E16.8,A1)') '"Ly": ', base_grid%Ly, ','
      write(json_id,'(8x,A6,1x,E16.8)'   ) '"Lz": ', base_grid%Lz
      write(json_id,'(4x,A3)') '},'

      write(json_id,'(4x,A12)') '"Solvers": {'
#ifdef MF
      write(json_id,'(8x,A8)') '"MF": 1,'
#else
      write(json_id,'(8x,A8)') '"MF": 0,'
#endif
#ifdef IBM
      write(json_id,'(8x,A8)') '"IBM": 1'
#else
      write(json_id,'(8x,A8)') '"IBM": 0'
#endif
      write(json_id,'(4x,A2)') '},'
      write(json_id,'(4x,A12,E18.6)') '"Timestep": ', dt
      write(json_id,'(A1)') '}'

      flush(json_id)
      close(json_id)
    endif

  end subroutine print_setup_json
  !========================================================================================

end module json
