program memory
   
    ! Test to check properly memory allocation and free.

    use lagrangian_solid_1D_mod
    use lagrangian_solid_2D_mod

    ! Variables
    type(lagrangian_solid_1D), target :: S1 
    type(lagrangian_solid_2D), target :: S2
    
    ! Test 1D lagrangian solid
    call S1%create('mesh.txt')
    call S1%destroy()

    ! Test 2D lagrangian solid
    call S2%loadFromFile('mesh2D.txt')
    call S2%destroy()

end program
