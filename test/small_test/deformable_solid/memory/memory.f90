program memory

    use lagrangian_solid_mod

    type(lagrangian_solid), target :: test_solid
    
    test_solid%is_open = .true.
    call test_solid%create('mesh.txt')
    
    call test_solid%destroy()

end program
