program mem_test

    use lagrangian_solid

    type(solid), target :: test_solid
    
    test_solid%is_open = .true.
    call test_solid%create('mesh.txt')
    
    call test_solid%destroy()

end program mem_test