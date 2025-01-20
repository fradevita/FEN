module lagrangian_solid_2D_mod

    ! This module contain the definition of the solid class.

    use precision_mod, only : dp
    use marker_mod
    use edge_mod
    use triangle_mod
    use lagrangian_solid_mod

    implicit none
    private
    
    public :: lagrangian_solid_2D

    type, extends(lagrangian_solid) :: lagrangian_solid_2D
        
        real(dp)                    :: A0 = 1.0_dp
        real(dp)                    :: V0 = 1.0_dp

        integer                     :: numberOfTriangles = 0
        real(dp)                    :: kb = 0.0_dp          
        
        type(triangle), allocatable :: triangles(:)
        
        integer, allocatable        :: edgeMassPointsIndex(:,:)
        logical, allocatable        :: duplicateEdges(:)
        integer, allocatable        :: edgeDuplicate(:)
        integer, allocatable        :: triangleEdgesIndex(:,:)
        integer, allocatable        :: edgeTrianglesIndex(:,:)
        integer, allocatable        :: oppositeVertexIndex(:,:)
        
        logical                     :: evaluateStretching = .false.
        logical                     :: evaluateBending = .false.
        logical                     :: evaluateAreaConstraint = .false.
        logical                     :: evaluateVolumeConstraint = .false.

        
        procedure(selfSub), pass(self), pointer :: computeAdditionalForces => Null()
        procedure(realFun), pass(self), pointer :: getAdditionalEnergy => Null()
        procedure(selfSub), pass(self), pointer :: computeBendingForces => Null()
        procedure(realFun), pass(self), pointer :: getBendingEnergy => Null()
    contains
    
        procedure, pass(self) :: loadFromFile

        procedure, pass(self) :: setMass
        procedure, pass(self) :: setStretchingConstant
        procedure, pass(self) :: setBendingConstant
        procedure, pass(self) :: setAreaConstant
        procedure, pass(self) :: setVolumeConstant
        
        procedure, pass(self) :: update
        procedure, pass(self) :: update_center_of_mass

        procedure, pass(self) :: compute_internal_forces
        
        procedure, pass(self) :: computeStretchingForces
        
        procedure, pass(self) :: computeAreaForces
        
        procedure, pass(self) :: getTotalArea
        procedure, pass(self) :: getTotalVolume
        procedure, pass(self) :: getKineticEnergy
        procedure, pass(self) :: getStretchingEnergy
        procedure, pass(self) :: getAreaEnergy
        procedure, pass(self) :: getPotentialEnergy
        procedure, pass(self) :: getBendingEnergyTheta2
        procedure, pass(self) :: getBendingEnergyCosTheta
        
        procedure, pass(Self) :: update_forcing_elements
        procedure, pass(self) :: interpolate_from_forcing_element_to_mass_point => interpolateFromCentroidToMassPoints

        procedure, pass(self) :: internalForcesIntegral
        procedure, pass(self) :: writeSTL
        procedure, pass(self) :: writeVTK
        procedure, pass(self) :: writeGNU
        procedure, pass(Self) :: printInfo
        
        procedure, pass(self) :: destroy
       
    end type lagrangian_solid_2D

    interface isPresent
        module procedure markerIsPresent
        module procedure edgeIsPresent
    end interface

    interface
        subroutine selfSub(self)
            import lagrangian_solid_2D
            class(lagrangian_solid_2D), intent(inout), target :: self  
        end subroutine selfSub
        function realFun(self) result(f)
            use precision_mod, only : dp
            import lagrangian_solid_2D
            class(lagrangian_solid_2D), intent(in ), target :: self
            real(dp)                          :: f   
        end function realFun
    end interface

contains
    
    !===============================================================================================
    subroutine loadFromFile(self, filename)

        use global_mod   , only : Ndim, myrank
        use euclidean_mod, only : dotProduct, ZEROV
        use utils_mod    , only : clamp, print_progress_bar

        ! In/Out variables
        class(lagrangian_solid_2D)    , intent(inout), target :: self
        character(len=*), intent(in)            :: filename

        ! Local variables
        integer                   :: fid, nl, error, l, t, i1, i2, i3, temp(2)
        real(dp)                  :: P(3)
        type(marker)              :: M(3)
        type(edge)                :: e(3)
        type(marker), allocatable :: Mtemp(:)
        character(len=5)          :: mn

        ! Conut number of input file lines
        open(newunit = fid, file = filename)
        nl = 0
        do 
            read(fid, *, iostat = error)
            if (error == -1) exit
            nl = nl + 1
        end do
        rewind(fid)
        
        ! Number of triangles is equal to nl/3
        if (mod(nl,3) /= 0) then
            print *, 'ERROR, number of vertex is not multiple of 3.'
            stop
        endif
        self%numberOfTriangles = nl/3
        allocate(self%triangles(self%numberOfTriangles))
        allocate(self%triangleEdgesIndex(self%numberOfTriangles,3))

        ! Number of mass points could vary, must be determined inside the 
        ! following loop
        ! Set an initial size equal to nl
        allocate(self%mass_points(nl))
        allocate(self%edgeMassPointsIndex(nl,2))
        allocate(self%edges(nl))
        allocate(self%duplicateEdges(nl))
        self%duplicateEdges = .false.
        allocate(self%edgeDuplicate(nl))
        self%edgeDuplicate = -1
        allocate(self%edgeTrianglesIndex(nl,2))

        ! This is needed to start the triangle loop below
        ! Read first three points
        do l = 1,3
            read(fid, *) P(1), P(2), P(3)
            write(mn,'(I0.5)') l
            self%mass_points(l) = marker(3, l, P)
            self%mass_points(l)%index = l
        end do

        ! Read first three edges
        do l = 1,3
            write(mn,'(I0.5)') l
            self%edges(l) = edge(self%mass_points(l), &
                        self%mass_points(merge(l+1,1,l<=2)), 'e_'//mn) 
            self%edgeMassPointsIndex(l,:) = [l, merge(l+1,1,l<=2)]
        end do
        self%number_of_mass_points = 3
        self%number_of_edges = 3
        
        ! Set first triangle
        self%triangleEdgesIndex(1,1:3) = [1, 2, 3]
        self%triangles(1)%e1 => self%edges(1)
        self%triangles(1)%e2 => self%edges(2)
        self%triangles(1)%e3 => self%edges(3)

        triangle_loop: do t = 2,self%numberOfTriangles

            ! Read the three verteces of triangle t
            do l = 1,3
                read(fid, *) P(1), P(2), P(3)
                write(mn,'(I0.5)') l
                M(l) = marker(3, -1, P)
                if (isPresent(M(l), self%mass_points, self%number_of_mass_points)) then
                    ! Go to the next point in the file
                else
                    ! Increase the counter of mass points
                    self%number_of_mass_points = self%number_of_mass_points + 1
                    write(mn,'(I0.5)') self%number_of_mass_points
                    ! Add to the list of mass points
                    self%mass_points(self%number_of_mass_points) = &
                            marker(3, self%number_of_mass_points, P)
                    self%mass_points(self%number_of_mass_points)%index = &
                            self%number_of_mass_points
                endif
            end do

            ! Set the three edges of the triangle
            do l = 1,3
                ! create the edge
                e(l) = edge(M(l), M(merge(l+1,1,l<=2)))
                
                ! NOTE: for now add all edges to the list self%edges
                ! and save if the new added one is duplicated
                ! Increase the counter of edges
                self%number_of_edges = self%number_of_edges + 1
                
                ! Add to the list of edges
                i1 = find_marker(e(l)%m1, self%mass_points)
                i2 = find_marker(e(l)%m2, self%mass_points)
                self%edges(self%number_of_edges) = &
                        edge(self%mass_points(i1), self%mass_points(i2))

                self%edgeMassPointsIndex(self%number_of_edges,:) = [i1, i2]
                self%triangleEdgesIndex(t,l) = self%number_of_edges
                self%duplicateEdges(self%number_of_edges) = &
                        edgeIsPresent(self%edges(self%number_of_edges), &
                        self%edges, self%number_of_edges-1)
                
                if (self%duplicateEdges(self%number_of_edges)) & 
                    self%edgeDuplicate(self%number_of_edges) = &
                        searchEdge(self%edges(self%number_of_edges), &
                        self%edges, self%number_of_edges)
            end do

            if (myrank == 0) call print_progress_bar('Performing triangle loop', t, self%numberOfTriangles)
        
        end do triangle_loop

        ! Set proper size of mass point array
        allocate(Mtemp(self%number_of_mass_points))
        Mtemp = self%mass_points(1:self%number_of_mass_points)
        deallocate(self%mass_points)
        allocate(self%mass_points(self%number_of_mass_points))
        self%mass_points = Mtemp
        deallocate(Mtemp)

        ! Set proper size of edges array
        deallocate(self%edges)
        allocate(self%edges(self%number_of_edges))
        do l = 1,self%number_of_edges
            i1 = self%edgeMassPointsIndex(l, 1)
            i2 = self%edgeMassPointsIndex(l, 2)
            self%edges(l) = edge(self%mass_points(i1), self%mass_points(i2))
        end do

        ! Finalize the triangles array
        do t = 1,self%numberOfTriangles
            i1 = self%triangleEdgesIndex(t,1)
            i2 = self%triangleEdgesIndex(t,2)
            i3 = self%triangleEdgesIndex(t,3)
            self%triangles(t) = triangle(self%edges(i1), self%edges(i2), self%edges(i3), 6)
        end do

        ! Create the edge - triangles index array
        do l = 1,self%number_of_edges
            self%edgeTrianglesIndex(l,:) = &
                    searchEdge_triangles(l, self%triangleEdgesIndex)
            if (self%duplicateEdges(l)) then
                if (self%edgeTrianglesIndex(l,2) > -1) then
                    print *, 'Error, second index has already been set'
                else
                    temp = searchEdge_triangles(self%edgeDuplicate(l), self%triangleEdgesIndex)
                    self%edgeTrianglesIndex(l,2) = temp(1)
                endif
            endif
        end do

        ! For bending force we also need to know the opposite mass point to edge
        ! i in triangle j
        allocate(self%oppositeVertexIndex(self%numberOfTriangles,3))
        do l = 1, self%numberOfTriangles
            self%oppositeVertexIndex(l,1) = &
                        self%edges(self%triangleEdgesIndex(l,3))%m1%index
            self%oppositeVertexIndex(l,2) = &
                        self%edges(self%triangleEdgesIndex(l,1))%m1%index
            self%oppositeVertexIndex(l,3) = &
                        self%edges(self%triangleEdgesIndex(l,2))%m1%index
        end do

        ! Setup the center of mass. Note that the degree of fredom is equal to 6 in 3D.
        self%center_of_mass = marker(6, 0)
        call self%update_center_of_mass()

        ! Setup the rest theta
        do l = 1,self%number_of_edges
            if (self%duplicateEdges(l))then
                i1 = self%edgeTrianglesIndex(l,1)
                i2 = self%edgeTrianglesIndex(l,2)
                self%edges(l)%theta0 = acos(clamp(dotProduct(&
                    self%triangles(i1)%n,self%triangles(i2)%n), -1.0_dp, 1.0_dp))
            endif
        end do

        ! Set initial Area
        self%A0 = self%getTotalArea()

        ! Set initial volume
        self%V0 = self%getTotalVolume()

        ! By default set rotation center equal to center of mass
        self%rotation_center => self%center_of_mass%X(1:Ndim)

#ifdef IBM
        ! To save memory, allocate forcing_elements arrays only when performing IBM simulations.
        allocate(self%forcing_elements(self%numberOfTriangles))
        do t = 1,self%numberOfTriangles
            self%forcing_elements(t)%C => self%triangles(t)%C
            self%forcing_elements(t)%A => self%triangles(t)%A
            self%forcing_elements(t)%n => self%triangles(t)%n
        end do
        allocate(self%p_s(self%numberOfTriangles))
        allocate(self%tau_s(6,self%numberOfTriangles))
#endif

    end subroutine
    !===============================================================================================

    !===============================================================================================
    function markerIsPresent(M, mass_points, lmax)

        ! In/Out variables
        logical                            :: markerIsPresent
        type(marker), intent(in)           :: M, mass_points(:)
        integer     , intent(in), optional :: lmax

        ! Local variables
        integer :: l, nl

        if (present(lmax)) then
            nl = lmax
        else
            nl = size(mass_points)
        endif

        markerIsPresent = .false.
        do l = 1,nl
            if (markersAreEqual(M, mass_points(l))) then
                markerIsPresent = .true.
                exit
            endif
        end do

    end function
    !===============================================================================================

    !===============================================================================================
    function find_marker(M, mass_points, lmax)

         ! In/Out variables
        integer                            :: find_marker
        type(marker), intent(in)           :: M, mass_points(:)
        integer     , intent(in), optional :: lmax

        ! Local variables
        integer :: l, nl

        if (present(lmax)) then
            nl = lmax
        else
            nl = size(mass_points)
        endif

        do l = 1,nl
            if (markersAreEqual(M, mass_points(l))) then
                find_marker = l
                exit
            endif
        end do

    end function
    !===============================================================================================

    !===============================================================================================
    function edgeIsPresent(e, edges, le)

        ! In/Out variables
        logical                :: edgeIsPresent
        type(edge), intent(in) :: e, edges(:)
        integer   , intent(in) :: le

        ! Local variables
        integer :: l

        edgeIsPresent = .false.
        do l = 1,le
            if ( (markersAreEqual(e%m1, edges(l)%m1) .and. &
                  markersAreEqual(e%m2, edges(l)%m2)) .or. &
                 (markersAreEqual(e%m1, edges(l)%m2) .and. &
                  markersAreEqual(e%m2, edges(l)%m1)) ) then
                edgeIsPresent = .true.
                exit
            endif
        end do

    end function
    !===============================================================================================

    !===============================================================================================
    function searchEdge(e, edges, le)

        ! In/Out variables
        integer                :: searchEdge
        type(edge), intent(in) :: e, edges(:)
        integer   , intent(in) :: le

        ! Local variables
        integer :: l

        searchEdge = -1
        do l = 1,le
            if ( (markersAreEqual(e%m1, edges(l)%m1) .and. &
                  markersAreEqual(e%m2, edges(l)%m2)) .or. &
                 (markersAreEqual(e%m1, edges(l)%m2) .and. &
                  markersAreEqual(e%m2, edges(l)%m1)) ) then
                searchEdge = l
                exit
            endif
        end do

        if (searchEdge == -1) then
            print *, 'Error: unable to find the duplicated edge'
            stop
        endif

    end function
    !===============================================================================================

    !===============================================================================================
    function searchEdge_triangles(le, triangle_edges_index) result(index)

        ! In/Out variables
        integer, intent(in) :: le
        integer, intent(in) :: triangle_edges_index(:,:)
        integer             :: index(2)

        ! Local variables
        integer :: t, nt, ii
 
        index = [-1, -1]
        nt = size(triangle_edges_index(:,1))
        ii = 1
        do t = 1,nt
            if (triangle_edges_index(t,1) == le .or. triangle_edges_index(t,2) == le .or. &
                triangle_edges_index(t,3) == le) then
                index(ii) = t
                ii = ii + 1
            endif
        end do

        if (ii == 1) then
            print *, 'ERROR: unable to find edge - triangles connection'
            stop
        endif

    end function
    !===============================================================================================

    !===============================================================================================
    subroutine setMass(self, M)

        ! In/Out variables
        class(lagrangian_solid_2D), intent(inout) :: self
        real(dp)                  , intent(in   ) :: M

        ! Local variables
        integer  :: n
        real(dp) :: mt
        
        do n = 1,self%number_of_mass_points
            self%mass_points(n)%m = 0.0_dp
        end do

        if (self%is_open) then
            mt = M/self%numberOfTriangles
            do n = 1,self%numberOfTriangles
                self%triangles(n)%v1%m = self%triangles(n)%v1%m + mt/3.0_dp
                self%triangles(n)%v2%m = self%triangles(n)%v2%m + mt/3.0_dp
                self%triangles(n)%v3%m = self%triangles(n)%v3%m + mt/3.0_dp
            end do
        else
            do n = 1,self%number_of_mass_points
                self%mass_points(n)%m = M/self%number_of_mass_points
            end do
        end if

    end subroutine setMass
    !===============================================================================================

    !===============================================================================================
    subroutine update(self)

        class(lagrangian_solid_2D), intent(inout), target :: self

        ! Local variables
        integer :: n

        do n = 1,self%number_of_edges
            call self%edges(n)%update_length
        end do

        do n = 1,self%numberOfTriangles
            call self%triangles(n)%updateArea()
            call self%triangles(n)%updateNorm()
            call self%triangles(n)%updateCentroid()
        end do

        !call self%update_center_of_mass

    end subroutine
    !===============================================================================================

    !===============================================================================================
    subroutine update_center_of_mass(self)

        class(lagrangian_solid_2D), intent(inout) :: self

        ! Local variables
        integer :: n

        self%center_of_mass%X(1:3) = [0.0_dp, 0.0_dp, 0.0_dp]
        
        do n = 1,self%number_of_mass_points
            self%center_of_mass%X(1:3) = self%center_of_mass%X(1:3) + self%mass_points(n)%X
        end do
        self%center_of_mass%X(1:3) = self%center_of_mass%X(1:3)/self%number_of_mass_points

    end subroutine
    !===============================================================================================

    !===============================================================================================
    subroutine setStretchingConstant(self, Ks, model)

        class(lagrangian_solid_2D)    , intent(inout)           :: self
        real(dp)        , intent(in   )           :: ks
        character(len=*), intent(in   ), optional :: model

        ! Local variables
        integer :: e, indexes(2)
        real(dp) :: l, A1, A2

        if (present(model)) then
            do e = 1,self%number_of_edges
                l = self%edges(e)%l
                indexes = self%edgeTrianglesIndex(e,:)
                A1 = self%triangles(indexes(1))%A
                if (indexes(2) == -1) then
                    !A2 = 0.0_dp
                    A2 = A1
                else
                    A2 = self%triangles(indexes(2))%A
                endif
                self%edges(e)%ks = Ks*(A1 + A2)/l**2
            end do
        else
            do e = 1,self%number_of_edges
                self%edges(e)%ks = Ks
            end do
        endif
        self%evaluateStretching = .true.

    end subroutine
    !===============================================================================================

    !===============================================================================================
    subroutine setBendingConstant(self, kb, model)

        ! Select the bending model and set the bending constant.

        ! In/Out variables
        class(lagrangian_solid_2D), intent(inout)           :: self
        real(dp)                  , intent(in   )           :: kb
        character(len=*)          , intent(in   ), optional :: model 

        ! Set the constant
        if (present(model)) then
            if (model == 'cos') then
                self%computeBendingForces => computeBendingForcesCosTheta
                self%getBendingEnergy => getBendingEnergyCosTheta
            else
                print *, 'WRONG BENDING MODEL. EXIT'
                stop
            endif
        else
            self%computeBendingForces => computeBendingForcesTheta2
            self%getBendingEnergy => getBendingEnergyTheta2
        endif
        self%kb = kb

        ! Activate the bending 
        self%evaluateBending = .true.

    end subroutine
    !===============================================================================================

    !===============================================================================================
    subroutine setAreaConstant(self, ka)

        class(lagrangian_solid_2D), intent(inout) :: self
        real(dp)    , intent(in   ) :: ka

        ! local variables
        integer :: n

        do n = 1,self%numberOfTriangles
            self%triangles(n)%ka = ka
        end do
        self%evaluateAreaConstraint = .true.

    end subroutine
    !===============================================================================================

    !===============================================================================================
    subroutine setVolumeConstant(self, kv)

        class(lagrangian_solid_2D), intent(inout) :: self
        real(dp)    , intent(in   ) :: kv

        ! local variables
        integer :: n

        do n = 1,self%numberOfTriangles
            self%triangles(n)%ka = kv
        end do
        self%evaluateVolumeConstraint = .true.

    end subroutine
    !===============================================================================================

    !===============================================================================================
    function getTotalArea(self) result(A)

        class(lagrangian_solid_2D), intent(in) :: self

        integer :: t
        real(dp) :: A

        A = 0.0_dp
        do t = 1,size(self%triangles)
            A = A + self%triangles(t)%A
        end do

    end function getTotalArea
    !===============================================================================================

    !===============================================================================================
    function getTotalVolume(self) result(V)

        class(lagrangian_solid_2D), intent(in) :: self

        integer :: t
        real(dp) :: V

        V = 0.0_dp
        do t = 1,size(self%triangles)
            V = V + self%triangles(t)%SignedVolume()
        end do

    end function getTotalVolume
    !===============================================================================================
 
    !===============================================================================================
    subroutine compute_internal_forces(self)

        use euclidean_mod, only : ZEROV

        ! In/Out variables
        class(lagrangian_solid_2D), intent(inout), target :: self

        ! Local variables
        integer :: n

        ! Set internal forces to zero
        do n = 1,self%number_of_mass_points
            self%mass_points(n)%Fi = ZEROV
        end do

        ! Elastic in-plane forces
        if (self%evaluateStretching) call self%computeStretchingForces()
        
        ! Bending forces
        if (self%evaluateBending) call self%computeBendingForces()
 
        ! Area-conservation forces
        if (self%evaluateAreaConstraint) call self%computeAreaForces()

        ! Any other potentials
        if (associated(self%computeAdditionalForces)) call self%computeAdditionalForces()

    end subroutine
    !===============================================================================================

    !===============================================================================================
    subroutine computeStretchingForces(self)

        use euclidean_mod, only : distanceVector

        ! In/Out variables
        class(lagrangian_solid_2D), intent(inout), target :: self
    
        ! Local variables
        integer             :: i
        type(edge), pointer :: e 
        real(dp)            :: dWsdL
        real(dp)            :: r12(3), r21(3)

        do i = 1,self%number_of_edges
            ! Add in-plane forces only for non-duplicated edges
            if (self%duplicateEdges(i) .eqv. .false.) then
                ! Select edge i
                e => self%edges(i)

                ! Distance vector from r1 to r2
                r12 = distanceVector(e%m2%X, e%m1%X)

                ! Distance vector from r2 to r1
                r21 = distanceVector(e%m1%X, e%m2%X)
                
                ! Stretching potential derivative w.r.t L
                dWsdL = e%ks*(e%l - e%l0)
                
                ! Add forces to mass points
                e%m1%Fi = e%m1%Fi - dWsdL*r12/e%l
                e%m2%Fi = e%m2%Fi - dWsdL*r21/e%l
            endif
        end do
    
    end subroutine
    !===============================================================================================

    !===============================================================================================
    subroutine computeBendingForcesTheta2(self)

        use euclidean_mod, only : distanceVector, crossProduct, distance, ZEROV
        use euclidean_mod, only : dotProduct, tilda_matrix

        ! The potential is defined as Wb = 0.5*kb*Lb*theta**2

        use utils_mod, only : clamp

        ! In/Out variables
        class(lagrangian_solid_2D), intent(inout), target :: self

        ! Local variables
        integer                  :: i, t1I, t2I
        type(edge), pointer      :: e
        type(triangle), pointer  :: t1, t2
        type(marker)  , pointer  :: mp1, mp2, mp3, mp4
        real(dp), dimension(3)   :: r31, r21, r24, r34, r23, r32
        real(dp), dimension(3)   :: m1, n1, m2, n2, dLbdr2, dLbdr3
        real(dp)                 :: im1norm, im2norm, n1Dotn2, theta, Lb
        real(dp)                 :: dthetadn1Dotn2, dtheta, dWbdLb, dWbdtheta
        real(dp), dimension(3,3) :: dn1dr21, dn1dr31, dn2dr24, dn2dr34
        real(dp)                 :: prefactor, an1(3,1), an2(3,1), tmp(3,1)

        do i = 1,self%number_of_edges
            
            ! Add bending forces only for duplicated edges
            if (self%duplicateEdges(i) .eqv. .true.) then

                ! Select local edge
                e => self%edges(i)

                ! Select the two triangles connected by edge e
                t1I = self%edgeTrianglesIndex(i,1)
                t1 => self%triangles(t1I)
                t2I = self%edgeTrianglesIndex(i,2)
                t2 => self%triangles(t2I)

                ! Select the four points
                mp2 => e%m2
                mp3 => e%m1

                ! Select mp1 from triangle 1
                if (self%triangleEdgesIndex(t1I,1) == i .or. &
                    self%triangleEdgesIndex(t1I,1) == self%edgeDuplicate(i)) then    
                    mp1 => self%mass_points(self%oppositeVertexIndex(t1I,1))
                else
                    if (self%triangleEdgesIndex(t1I,2) == i .or. &
                        self%triangleEdgesIndex(t1I,2) == self%edgeDuplicate(i)) then    
                        mp1 => self%mass_points(self%oppositeVertexIndex(t1I,2))
                    else
                        mp1 => self%mass_points(self%oppositeVertexIndex(t1I,3))
                    endif
                endif

                ! Select mp4 from triangle 2
                if (self%triangleEdgesIndex(t2I,1) == i .or. &
                    self%triangleEdgesIndex(t2I,1) == self%edgeDuplicate(i)) then
                    mp4 => self%mass_points(self%oppositeVertexIndex(t2I,1))
                else
                    if (self%triangleEdgesIndex(t2I,2) == i .or. &
                        self%triangleEdgesIndex(t2I,2) == self%edgeDuplicate(i)) then
                        mp4 => self%mass_points(self%oppositeVertexIndex(t2I,2))
                    else
                        mp4 => self%mass_points(self%oppositeVertexIndex(t2I,3))
                    endif
                endif

                ! Evaluate distasnce  vectors and normal vectors
                r31 = distanceVector(mp1%X, mp3%X)
                r21 = distanceVector(mp1%X, mp2%X)
                m1 = crossProduct(r31,r21)
                im1norm = 1.0_dp/distance(m1, ZEROV)
                n1 = m1*im1norm
                an1(:,1) = n1

                r24 = distanceVector(mp4%X, mp2%X)
                r34 = distanceVector(mp4%X, mp3%X)
                m2 = crossProduct(r24,r34)
                im2norm = 1.0_dp/distance(m2, ZEROV)
                n2 = m2*im2norm
                an2(:,1) = n2 

                ! Angle between two triangles
                n1Dotn2 = clamp(dotProduct(n1,n2), -1.0_dp, 1.0_dp)
                theta = acos(n1Dotn2)

                ! Derivative of theta w.r.t n1_dot_n2
                dthetadn1Dotn2 = -1.0_dp*sqrt(1.0_dp - n1Dotn2**2)/ &
                                    (sqrt(1.0_dp - n1Dotn2**2)**2 + 1.0e-12)

                ! lenght of the shared edge
                Lb = e%l

                ! Bengind potential derivatives w.r.t Lb and theta
                dtheta = theta - e%theta0
                dWbdLb = 0.5_dp*self%kb*(dtheta)**2
                dWbdtheta = self%kb*Lb*(dtheta)

                tmp(:,1) = crossProduct(n1,r31)
                dn1dr21 = -(tilda_matrix(r31) + & 
                            matmul(an1, transpose(tmp)))*im1norm
                
                tmp(:,1) = crossProduct(n1, r21)
                dn1dr31 = (tilda_matrix(r21) + &
                            matmul(an1, transpose(tmp)))*im1norm
                
                tmp(:,1) = crossProduct(n2, r34)
                dn2dr24 = (tilda_matrix(r34) + &
                            matmul(an2, transpose(tmp)))*im2norm
                
                tmp(:,1) = crossProduct(n2, r24)
                dn2dr34 = -(tilda_matrix(r24) + &
                            matmul(an2, transpose(tmp)))*im2norm

                ! Evaluate forces
                prefactor = dWbdtheta*dthetadn1Dotn2
                mp1%Fi = mp1%Fi - &
                        prefactor*(matmul(transpose(-dn1dr21 - dn1dr31), n2))
                r23 = distanceVector(mp3%X, mp2%X)
                dLbdr2 = r23/distance(r23, ZEROV)
                mp2%Fi = mp2%Fi - (dWbdLb*dLbdr2 + &
                        prefactor*(matmul(transpose(dn1dr21), n2) +&
                                   matmul(transpose(dn2dr24), n1)))
                r32 = distanceVector(mp2%X, mp3%X)
                dLbdr3 = r32/distance(r32, ZEROV)
                mp3%Fi = mp3%Fi - (dWbdLb*dLbdr3 + &
                         prefactor*(matmul(transpose(dn1dr31), n2) +&
                                    matmul(transpose(dn2dr34), n1)))
                mp4%Fi = mp4%Fi - &
                        prefactor*(matmul(transpose(-dn2dr24 - dn2dr34), n1))

            endif
        end do
    
    end subroutine
    !===============================================================================================

    !===============================================================================================
    subroutine computeBendingForcesCosTheta(self)

        use euclidean_mod, only : distanceVector, crossProduct, distance, ZEROV
        use euclidean_mod, only : dotProduct, tilda_matrix

        ! The potential is defined as Wb = kb*[1. - cos(theta - theta_0)]

        use utils_mod, only : clamp

        ! In/Out variables
        class(lagrangian_solid_2D), intent(inout), target :: self

        ! Local variables
        integer                  :: i, t1I, t2I
        type(edge), pointer      :: e
        type(triangle), pointer  :: t1, t2
        type(marker)  , pointer  :: mp1, mp2, mp3, mp4
        real(dp), dimension(3)   :: r31, r21, r24, r34, r23, r32
        real(dp), dimension(3)   :: m1, n1, m2, n2, dLbdr2, dLbdr3
        real(dp)                 :: im1norm, im2norm, n1Dotn2, theta, Lb
        real(dp)                 :: dthetadn1Dotn2, dtheta, dWbdLb, dWbdtheta
        real(dp), dimension(3,3) :: dn1dr21, dn1dr31, dn2dr24, dn2dr34
        real(dp)                 :: prefactor, an1(3,1), an2(3,1), tmp(3,1)

        do i = 1,self%number_of_edges
            
            ! Add bending forces only for duplicated edges
            if (self%duplicateEdges(i) .eqv. .true.) then

                ! Select local edge
                e => self%edges(i)

                ! Select the two triangles connected by edge e
                t1I = self%edgeTrianglesIndex(i,1)
                t1 => self%triangles(t1I)
                t2I = self%edgeTrianglesIndex(i,2)
                t2 => self%triangles(t2I)

                ! Select the four points
                mp2 => e%m2
                mp3 => e%m1

                ! Select mp1 from triangle 1
                if (self%triangleEdgesIndex(t1I,1) == i .or. &
                    self%triangleEdgesIndex(t1I,1) == self%edgeDuplicate(i)) then    
                    mp1 => self%mass_points(self%oppositeVertexIndex(t1I,1))
                else
                    if (self%triangleEdgesIndex(t1I,2) == i .or. &
                        self%triangleEdgesIndex(t1I,2) == self%edgeDuplicate(i)) then    
                        mp1 => self%mass_points(self%oppositeVertexIndex(t1I,2))
                    else
                        mp1 => self%mass_points(self%oppositeVertexIndex(t1I,3))
                    endif
                endif

                ! Select mp4 from triangle 2
                if (self%triangleEdgesIndex(t2I,1) == i .or. &
                    self%triangleEdgesIndex(t2I,1) == self%edgeDuplicate(i)) then
                    mp4 => self%mass_points(self%oppositeVertexIndex(t2I,1))
                else
                    if (self%triangleEdgesIndex(t2I,2) == i .or. &
                        self%triangleEdgesIndex(t2I,2) == self%edgeDuplicate(i)) then
                        mp4 => self%mass_points(self%oppositeVertexIndex(t2I,2))
                    else
                        mp4 => self%mass_points(self%oppositeVertexIndex(t2I,3))
                    endif
                endif

                ! Evaluate distasnce vectors and normal vectors
                r31 = distanceVector(mp1%X, mp3%X)
                r21 = distanceVector(mp1%X, mp2%X)
                m1 = crossProduct(r31,r21)
                im1norm = 1.0_dp/distance(m1, ZEROV)
                n1 = m1*im1norm
                an1(:,1) = n1

                r24 = distanceVector(mp4%X, mp2%X)
                r34 = distanceVector(mp4%X, mp3%X)
                m2 = crossProduct(r24,r34)
                im2norm = 1.0_dp/distance(m2, ZEROV)
                n2 = m2*im2norm
                an2(:,1) = n2 

                ! Angle between two triangles
                n1Dotn2 = clamp(dotProduct(n1,n2), -1.0_dp, 1.0_dp)
                theta = acos(n1Dotn2)

                ! Derivative of theta w.r.t n1_dot_n2
                dthetadn1Dotn2 = -1.0_dp*sqrt(1.0_dp - n1Dotn2**2)/ &
                                    (sqrt(1.0_dp - n1Dotn2**2)**2 + 1.0e-12)

                ! lenght of the shared edge
                Lb = e%l

                ! Bending potential derivatives w.r.t Lb and theta
                dtheta = theta - e%theta0
                dWbdtheta = self%kb*sin(dtheta)
                dWbdLb = 0.0_dp

                block
                    real(dp)               :: betab, b11, b12, b22, modxsi, modzed, costheta, sintheta
                    real(dp), dimension(3) :: xsi, zed, a21, a31, a34, a24, a32, a13, a42, a23
                    type(marker), pointer  :: temp

                    mp2 => e%m1
                    mp3 => e%m2
    
                    a21 = distanceVector(mp1%X, mp2%x)
                    a31 = distanceVector(mp1%X, mp3%x)
                    a34 = distanceVector(mp4%X, mp3%x)
                    a24 = distanceVector(mp4%X, mp2%x)

                    a32 = distanceVector(mp2%X, mp3%X)
                    a13 = distanceVector(mp3%X, mp1%X)
                    a42 = distanceVector(mp2%X, mp4%X)
                    a23 = distanceVector(mp3%X, mp2%X)

                    xsi = crossProduct(a21, a31)
                    modxsi = distance(xsi, ZEROV)
                    zed = crossProduct(a34, a24)
                    modzed = distance(zed, ZEROV)

                    costheta = clamp(dotProduct(xsi, zed)/modxsi/modzed, -1.0_dp, 1.0_dp)
                    theta = acos(costheta)
                    sintheta = sqrt(1.0_dp - costheta*costheta)
                 
                    if (dotProduct(distanceVector(zed, xsi), &
                            distanceVector(t2%C%X,t1%C%X)) < 0.0_dp ) then
                        sintheta = -sintheta
                    end if
                    sintheta = sin(theta)
                    
                    betab = self%kb*(sintheta*cos(e%theta0) - costheta*sin(e%theta0))* &
                                sqrt(1. - costheta**2)/((sqrt(1. - costheta**2))**2 + 1.0e-12_dp)
                    
                    b11 = -betab*costheta/modxsi**2
                    b12 =  betab/(modxsi*modzed)
                    b22 = -betab*costheta/modzed**2

                    mp1%Fi = mp1%Fi + b11*crossProduct(xsi, a32) + b12*crossProduct(zed, a32)
                    mp2%Fi = mp2%Fi + b11*crossProduct(xsi, a13) + b12*(crossProduct(xsi, a34) + &
                                    crossProduct(zed, a13)) + b22*crossProduct(zed, a34)
                    mp3%Fi = mp3%Fi + b11*crossProduct(xsi, a21) + b12*(crossProduct(xsi, a42) + &
                                    crossProduct(zed, a21)) + b22*crossProduct(zed, a42)
                    mp4%Fi = mp4%Fi + b12*crossProduct(xsi, a23) + b22*crossProduct(zed, a23)
                
                end block

            endif
        end do
    
    end subroutine
    !===============================================================================================

    !===============================================================================================
    subroutine computeAreaForces(self)

        use euclidean_mod, only : distanceVector, crossProduct, ex, ey, ez
  
        ! In/Out variables
        class(lagrangian_solid_2D), intent(inout), target :: self

        ! Local variables
        integer                 :: i
        real(dp)                :: dWadA, B(3,3)
        type(marker), pointer   :: mp1, mp2, mp3
        real(dp), dimension(3)  :: r21, r31
        type(triangle), pointer :: t

        do i = 1,self%numberOfTriangles

            ! Select local triangle
            t => self%triangles(i)

            ! Potential derivative with respect to A
            dWadA = t%ka*(t%A - t%A0)/t%A0

            ! Select the three vertex
            mp1 => t%e1%m1
            mp2 => t%e1%m2
            if (markersAreEqual(t%e2%m1, t%e1%m1) .or. &
                markersAreEqual(t%e2%m1, t%e1%m2)) then
                ! Select the other point
                mp3 => t%e2%m2
                ! Check that is a valid one
                if (markersAreEqual(t%e2%m2, t%e1%m1) .or. &
                    markersAreEqual(t%e2%m2, t%e1%m2)) then
                    print *, 'ERROR: also m2 is equal to other vertex'
                    stop
                endif
            else
                mp3 => t%e2%m1
            endif

            ! Distance vector
            r21 = distanceVector(mp1%X, mp2%X)
            r31 = distanceVector(mp1%X, mp3%X)

            ! Evaluate force on mp1
            B(1,:) = crossProduct(-ex, r31)
            B(2,:) = crossProduct(-ey, r31)
            B(3,:) = crossProduct(-ez, r31)
            mp1%Fi = mp1%Fi - 0.5_dp*dWadA*(matmul(B, t%n))

            B(1,:) = crossProduct(r21, -ex)
            B(2,:) = crossProduct(r21, -ey)
            B(3,:) = crossProduct(r21, -ez)
            mp1%Fi = mp1%Fi - 0.5_dp*dWadA*(matmul(B, t%n))
            
            ! Evaluate force on r2
            B(1,:) = crossProduct(ex, r31)
            B(2,:) = crossProduct(ey, r31)
            B(3,:) = crossProduct(ez, r31)
            mp2%Fi = mp2%Fi - 0.5_dp*dWadA*(matmul(B, t%n))
           
            ! Evaluate force on r3
            B(1,:) = crossProduct(r21, ex)
            B(2,:) = crossProduct(r21, ey)
            B(3,:) = crossProduct(r21, ez)
            mp3%Fi = mp3%Fi - 0.5_dp*dWadA*(matmul(B, t%n))
        
        end do

    end subroutine
    !===============================================================================================
    
    !===============================================================================================
    function getKineticEnergy(self) result(Ke)

        use euclidean_mod, only : distance, ZEROV

        ! In/Out variables
        class(lagrangian_solid_2D), intent(in), target :: self
        real(dp)                         :: Ke

        ! Local variables
        integer                :: i
        type(marker), pointer  :: mp

        Ke = 0.0_dp
        do i = 1,self%number_of_mass_points
            mp => self%mass_points(i)
            Ke = Ke + 0.5_dp*mp%m*(distance(mp%V, ZEROV))**2
        end do

    end function
    !===============================================================================================

    !===============================================================================================
    function getStretchingEnergy(self) result(Ws)

        ! In/Out variables
        class(lagrangian_solid_2D), intent(in), target :: self
        real(dp)                         :: Ws

        ! Local variables
        integer             :: i
        type(edge), pointer :: e

        Ws = 0.0_dp
        do i = 1,self%number_of_edges
            if (self%duplicateEdges(i) .eqv. .false.) then
                e => self%edges(i)
                Ws = Ws + 0.5_dp*e%ks*(e%L - e%L0)**2
            endif
        enddo 

    end function
    !===============================================================================================

    !===============================================================================================
    function getAreaEnergy(self) result(Wa)

        ! In/Out variables
        class(lagrangian_solid_2D), intent(in), target :: self
        real(dp)                         :: Wa

        ! Local variables
        integer                 :: i
        type(triangle), pointer :: t

        Wa = 0.0_dp
        do i = 1,self%numberOfTriangles
            t => self%triangles(i)
            Wa = Wa + 0.5_dp*t%ka*((t%A - t%A0)/t%A0)**2*t%A0
        end do

    end function
    !===============================================================================================

    !===============================================================================================
    function getBendingEnergyTheta2(self) result(Wb)

        use utils_mod    , only : clamp
        use euclidean_mod, only : distanceVector, crossProduct, distance, ZEROV
        use euclidean_mod, only : dotProduct

        ! In/Out variables
        class(lagrangian_solid_2D), intent(in), target :: self
        real(dp)                         :: Wb

        ! Local variables
        integer                 :: i, t1I, t2I
        type(edge)    , pointer :: e
        type(triangle), pointer :: t1, t2
        type(marker)  , pointer :: mp1, mp2, mp3, mp4
        real(dp), dimension(3)  :: r31, r21, m1, n1, r24, r34, m2, n2
        real(dp)                :: theta

        Wb = 0.0_dp
        do i = 1,self%number_of_edges
            
            if (self%duplicateEdges(i) .eqv. .true.) then
                  ! Select local edge
                e => self%edges(i)

                ! Select the two triangles connected by edge e
                t1I = self%edgeTrianglesIndex(i,1)
                t1 => self%triangles(t1I)
                t2I = self%edgeTrianglesIndex(i,2)
                t2 => self%triangles(t2I)

                ! Select the four points
                mp2 => e%m2
                mp3 => e%m1

                ! Select mp1 from triangle 1
                if (self%triangleEdgesIndex(t1I,1) == i .or. &
                    self%triangleEdgesIndex(t1I,1) == self%edgeDuplicate(i)) then    
                    mp1 => self%mass_points(self%oppositeVertexIndex(t1I,1))
                else
                    if (self%triangleEdgesIndex(t1I,2) == i .or. &
                        self%triangleEdgesIndex(t1I,2) == self%edgeDuplicate(i)) then    
                        mp1 => self%mass_points(self%oppositeVertexIndex(t1I,2))
                    else
                        mp1 => self%mass_points(self%oppositeVertexIndex(t1I,3))
                    endif
                endif

                ! Select mp4 from triangle 2
                if (self%triangleEdgesIndex(t2I,1) == i .or. &
                    self%triangleEdgesIndex(t2I,1) == self%edgeDuplicate(i)) then
                    mp4 => self%mass_points(self%oppositeVertexIndex(t2I,1))
                else
                    if (self%triangleEdgesIndex(t2I,2) == i .or. &
                        self%triangleEdgesIndex(t2I,2) == self%edgeDuplicate(i)) then
                        mp4 => self%mass_points(self%oppositeVertexIndex(t2I,2))
                    else
                        mp4 => self%mass_points(self%oppositeVertexIndex(t2I,3))
                    endif
                endif

                ! Evaluate distasnce  vectors and normal vectors
                r31 = distanceVector(mp1%X, mp3%X)
                r21 = distanceVector(mp1%X, mp2%X)
                m1 = crossProduct(r31,r21)
                n1 = m1/distance(m1, ZEROV)

                r24 = distanceVector(mp4%X, mp2%X)
                r34 = distanceVector(mp4%X, mp3%X)
                m2 = crossProduct(r24,r34)
                n2 = m2/distance(m2, ZEROV) 

                ! Angle between two triangles
                theta = acos(clamp(dotProduct(n1, n2), -1.0_dp, 1.0_dp))
                
                ! Add energy contribution of this edge
                Wb = Wb + 0.5_dp*self%kb*e%L*(theta - e%theta0)**2
            endif
        end do

    end function
    !===============================================================================================

    !===============================================================================================
    function getBendingEnergyCosTheta(self) result(Wb)

        use utils_mod    , only : clamp
        use euclidean_mod, only : distanceVector, crossProduct, distance, ZEROV
        use euclidean_mod, only : dotProduct

        ! In/Out variables
        class(lagrangian_solid_2D), intent(in), target :: self
        real(dp)                         :: Wb

        ! Local variables
        integer                 :: i, t1I, t2I
        type(edge)    , pointer :: e
        type(triangle), pointer :: t1, t2
        type(marker)  , pointer :: mp1, mp2, mp3, mp4
        real(dp), dimension(3)  :: r31, r21, m1, n1, r24, r34, m2, n2
        real(dp)                :: theta

        Wb = 0.0_dp
        do i = 1,self%number_of_edges
            
            if (self%duplicateEdges(i) .eqv. .true.) then
                  ! Select local edge
                e => self%edges(i)

                ! Select the two triangles connected by edge e
                t1I = self%edgeTrianglesIndex(i,1)
                t1 => self%triangles(t1I)
                t2I = self%edgeTrianglesIndex(i,2)
                t2 => self%triangles(t2I)

                ! Select the four points
                mp2 => e%m2
                mp3 => e%m1

                ! Select mp1 from triangle 1
                if (self%triangleEdgesIndex(t1I,1) == i .or. &
                    self%triangleEdgesIndex(t1I,1) == self%edgeDuplicate(i)) then    
                    mp1 => self%mass_points(self%oppositeVertexIndex(t1I,1))
                else
                    if (self%triangleEdgesIndex(t1I,2) == i .or. &
                        self%triangleEdgesIndex(t1I,2) == self%edgeDuplicate(i)) then    
                        mp1 => self%mass_points(self%oppositeVertexIndex(t1I,2))
                    else
                        mp1 => self%mass_points(self%oppositeVertexIndex(t1I,3))
                    endif
                endif

                ! Select mp4 from triangle 2
                if (self%triangleEdgesIndex(t2I,1) == i .or. &
                    self%triangleEdgesIndex(t2I,1) == self%edgeDuplicate(i)) then
                    mp4 => self%mass_points(self%oppositeVertexIndex(t2I,1))
                else
                    if (self%triangleEdgesIndex(t2I,2) == i .or. &
                        self%triangleEdgesIndex(t2I,2) == self%edgeDuplicate(i)) then
                        mp4 => self%mass_points(self%oppositeVertexIndex(t2I,2))
                    else
                        mp4 => self%mass_points(self%oppositeVertexIndex(t2I,3))
                    endif
                endif

                ! Evaluate distasnce  vectors and normal vectors
                r31 = distanceVector(mp1%X, mp3%X)
                r21 = distanceVector(mp1%X, mp2%X)
                m1 = crossProduct(r31,r21)
                n1 = m1/distance(m1, ZEROV)

                r24 = distanceVector(mp4%X, mp2%X)
                r34 = distanceVector(mp4%X, mp3%X)
                m2 = crossProduct(r24,r34)
                n2 = m2/distance(m2, ZEROV) 

                ! Angle between two triangles
                theta = acos(clamp(dotProduct(n1, n2), -1.0_dp, 1.0_dp))
                
                ! Add energy contribution of this edge
                Wb = Wb + self%kb*(1.0_dp - cos(theta - e%theta0))
            endif
        end do

    end function
    !===============================================================================================

    !===============================================================================================
    function getPotentialEnergy(self) result(W)

        ! In/Out variables
        class(lagrangian_solid_2D), intent(in), target :: self
        real(dp)                         :: W

        ! Local variables
        real(dp) :: Ws, Wb, Wa, Wadd

        W = 0.0_dp
        Ws = 0.0_dp
        Wa = 0.0_dp
        Wb = 0.0_dp

        ! Add each potential contribution
        if (self%evaluateStretching) Ws = self%getStretchingEnergy()
        if (self%evaluateBending) Wb = self%getBendingEnergy()
        if (self%evaluateAreaConstraint) Wa = self%getAreaEnergy()

        W = Ws + Wa + Wb

        if (associated(self%getAdditionalEnergy)) then
            Wadd = self%getAdditionalEnergy()
            W = W + Wadd
        endif

    end function
    !===============================================================================================

    !===============================================================================================
    subroutine interpolateFromCentroidToMassPoints(self)

        use euclidean_mod, only : distance, ZEROV

        ! In/Out variables
        class(lagrangian_solid_2D), intent(inout), target :: self

        ! Local variables
        integer                 :: i
        real(dp)                :: w1, w2, w3, wtot
        type(triangle), pointer :: t

        do i = 1,self%number_of_mass_points
            self%mass_points(i)%Fh = ZEROV
            self%mass_points(i)%Fe = ZEROV
        end do

        do i = 1,self%numberOfTriangles
            t => self%triangles(i)
            w1 = distance(t%v1%X, t%c%X)
            w2 = distance(t%v2%X, t%c%X)
            w3 = distance(t%v3%X, t%c%X)
            wtot = w1 + w2 + w3
            t%v1%Fe(1:3) = t%v1%Fe(1:3) + t%c%Fe(1:3)*w1/wtot
            t%v2%Fe(1:3) = t%v2%Fe(1:3) + t%c%Fe(1:3)*w2/wtot
            t%v3%Fe(1:3) = t%v3%Fe(1:3) + t%c%Fe(1:3)*w3/wtot

            t%v1%Fh(1:3) = t%v1%Fh(1:3) + t%c%Fh(1:3)*w1/wtot
            t%v2%Fh(1:3) = t%v2%Fh(1:3) + t%c%Fh(1:3)*w2/wtot
            t%v3%Fh(1:3) = t%v3%Fh(1:3) + t%c%Fh(1:3)*w3/wtot
        end do

    end subroutine
    !===============================================================================================
    
    !===============================================================================================
    function internalForcesIntegral(self) result(int_fi)

        class(lagrangian_solid_2D), intent(in) :: self
        real(dp)                 :: int_fi(3)

        ! Local varibales
        integer :: i

        int_fi = 0.0_dp
        do i = 1,self%number_of_mass_points
            int_fi = int_fi + self%mass_points(i)%Fi
        end do
        int_fi = int_fi/self%number_of_mass_points

    end function
    !===============================================================================================
   
    ! !===============================================================================================
    ! subroutine velocityVerlet(self, dt)

    !     ! In/Out variables
    !     class(lagrangian_solid_2D), intent(inout), target :: self
    !     real(dp)    , intent(in   )         :: dt

    !     ! Local variables
    !     integer               :: i
    !     type(marker), pointer :: mp

    !     ! **** First step of the velocity verlet

    !     ! Evaluate internal forces from potential
    !     call self%computeInternalForces()
       
    !     do i = 1,self%number_of_mass_points
    !         ! Select current mass point
    !         mp => self%massPoints(i) 
    !         ! Acceleration at n
    !         mp%a = (mp%Fe + mp%Fi)/mp%m
    !         ! Velocity at n + 1/2
    !         mp%v = mp%v + 0.5_dp*dt*mp%a
    !         ! Position at n + 1
    !         mp%X = mp%X + dt*mp%v 
    !     end do

    !     ! Enforce BC
    !     if (associated(self%applyConstraints)) then
    !        call self%applyConstraints()
    !     endif

    !     ! Update the solid with the new mass points position
    !     call self%update()

    !     ! **** Second step of the velocity verlet

    !     ! Get new internal forces
    !     call self%computeInternalForces()

    !     do i = 1,self%number_of_mass_points
    !         ! Select current mass point
    !         mp => self%massPoints(i)
    !         ! Accleration at n + 1
    !         mp%a = (mp%Fe + mp%Fi)/mp%m
    !         ! Veloicyt at n + 1
    !         mp%v = mp%v*(1.0_dp - self%gamma) + 0.5_dp*dt*mp%a
    !     end do

    ! end subroutine
    ! !===============================================================================================

    ! !===============================================================================================
    ! subroutine expImp(self, dt)

    !     ! In/Out variables
    !     class(lagrangian_solid_2D), intent(inout), target :: self
    !     real(dp)    , intent(in   )         :: dt

    !     ! Local variables
    !     integer               :: i
    !     real(dp)              :: a, b, c
    !     type(marker), pointer :: mp

    !     ! Some precomputed coefficient
    !     a = 0.5_dp*dt/self%mass_points(1)%m !assuming all points with same mass
    !     b = a*self%gamma
    !     c = 1.0_dp/(1.0_dp + b)

    !     ! **** First step is explicit
    !     ! Evaluate internal forces from potential at step n
    !     call self%computeInternalForces()
       
    !     do i = 1,self%number_of_mass_points
    !         ! Select current mass point
    !         mp => self%mass_points(i) 
    !         ! Velocity at n + 1/2
    !         mp%v = mp%v*(1.0_dp - b) + a*(mp%Fe + mp%Fi)
    !         ! Position at n + 1
    !         mp%X = mp%X + dt*mp%v 
    !     end do

    !     ! Enforce BC on position at n + 1
    !     if (associated(self%applyConstraints)) then
    !        call self%applyConstraints()
    !     endif

    !     ! Update the solid with the new mass points position
    !     call self%update()

    !     ! **** Second step is implicit

    !     ! Evaluate internal forces from potential at step n+1
    !     call self%computeInternalForces()

    !     do i = 1,self%number_of_mass_points
    !         ! Select current mass point
    !         mp => self%mass_points(i)
    !         ! Veloicyt at n + 1
    !         mp%v = (mp%v + a*(mp%Fe + mp%Fi))*c
    !     end do

    ! end subroutine
    ! !===============================================================================================

    !===============================================================================================
    subroutine writeGNU(self, filename)

        ! In/Out variables
        class(lagrangian_solid_2D)    , intent(in), target :: self
        character(len=*), intent(in)         :: filename

        ! Local variables
        integer                 :: i, file_id
        type(triangle), pointer :: t

        open(newunit = file_id, file = filename)
        do i = 1,self%numberOfTriangles
            
            ! Select local triangle
            t => self%triangles(i)

            write(file_id,*) t%e1%m1%X
            write(file_id,*) t%e1%m2%X
            write(file_id,*) ''
            
            write(file_id,*) t%e2%m1%X
            write(file_id,*) t%e2%m2%X
            write(file_id,*) ''
            
            write(file_id,*) t%e3%m1%X
            write(file_id,*) t%e3%m2%X
            write(file_id,*) ''

            write(file_id,*) ''
        end do
        close(file_id)

    end subroutine
    !===============================================================================================

    !===============================================================================================
    subroutine writeSTL(self, filename)

        ! In/out varibales
        class(lagrangian_solid_2D)    , intent(in), target :: self
        character(len=*), intent(in)         :: filename

        ! Local variables
        integer                 :: i, file_id
        type(triangle), pointer :: t

        open(newunit = file_id, file = filename)
        write(file_id,'(A11)') 'solid Plate'
        do i = 1,self%numberOfTriangles
            
            ! Select local triangle
            t => self%triangles(i)

            write(file_id,'(A12,3(1x,F16.8))') 'facet normal', t%n
            write(file_id,'(A14)') '    outer loop'
            write(file_id,'(A14,3(1x,F16.8))') '        vertex', t%e1%m1%X
            write(file_id,'(A14,3(1x,F16.8))') '        vertex', t%e2%m1%X
            write(file_id,'(A14,3(1x,F16.8))') '        vertex', t%e3%m1%X
            write(file_id,'(A11)') '    endloop'
            write(file_id,'(A8)') 'endfacet'
        end do
        write(file_id,'(A14)') 'endsolid Plate'
        close(file_id)

    end subroutine writeSTl
    !===============================================================================================

    !===============================================================================================
    subroutine writeVTK(self, filename)

        ! In/out varibales
        class(lagrangian_solid_2D)    , intent(in), target :: self
        character(len=*), intent(in)         :: filename

        ! Local variables
        integer                 :: i, file_id
        type(triangle), pointer :: t

        open(newunit = file_id, file = filename)
        write(file_id,'(A26)'     ) '# vtk DataFile Version 2.0'
        write(file_id,'( A5)'     ) 'Solid'
        write(file_id,'( A5)'     ) 'ASCII'
        write(file_id,'(A16)'     ) 'DATASET POLYDATA'
        write(file_id,'(A7,I7,A6)') 'POINTS ', self%number_of_mass_points, ' float'
        do i = 1,self%number_of_mass_points
            write(file_id, '(E16.8,1x,E16.8,1x,E16.8)') self%mass_points(i)%X(1:3)
        end do
        write(file_id,'(A9,I8,1x,I8)') 'POLYGONS ', self%numberOfTriangles, &
                                                    self%numberOfTriangles*4
        do i = 1,self%numberOfTriangles
            write(file_id,'(I1,1x,I7,1x,I7,1x,I7)') 3, self%triangles(i)%e1%m1%index-1, &
                                                       self%triangles(i)%e1%m2%index-1, &
                                                       self%triangles(i)%e2%m2%index-1 
        end do                          
        write(file_id,'(A10,I7)') 'CELL_DATA ', self%numberOfTriangles
        write(file_id,'(A26)') 'NORMALS cell_normals float'
        do i = 1,self%numberOfTriangles
            write(file_id,'(E16.8,1x,E16.8,1x,E16.8)') self%triangles(i)%n(1:3) 
        end do
        write(file_id,'(A15)') 'VECTORS V float'
        do i = 1,self%numberOfTriangles
            write(file_id,'(E16.8,1x,E16.8,1x,E16.8)') self%triangles(i)%C%V(1:3)
        end do
        write(file_id,'(A15)') 'SCALARS P float'
        write(file_id,'(A20)') 'LOOKUP_TABLE default'
        do i = 1,self%numberOfTriangles
            write(file_id,'(E16.8,1x,E16.8,1x,E16.8)') self%p_s(i)
        end do
        write(file_id,'(A17)') 'TENSORS TAU float'
        do i = 1,self%numberOfTriangles
            write(file_id,'(E16.8,1x,E16.8,1x,E16.8)') self%tau_s(1,i), self%tau_s(2,i), self%tau_s(4,i)
            write(file_id,'(E16.8,1x,E16.8,1x,E16.8)') self%tau_s(2,i), self%tau_s(3,i), self%tau_s(5,i)
            write(file_id,'(E16.8,1x,E16.8,1x,E16.8)') self%tau_s(4,i), self%tau_s(5,i), self%tau_s(6,i)
        end do
        
        close(file_id)

    end subroutine writeVTK
    !===============================================================================================

    !===============================================================================================
    subroutine printInfo(self)

        class(lagrangian_solid_2D), intent(in) :: self

        integer :: n

        print *, 'Number of mass points: ', self%number_of_mass_points
        print *, 'Number of edges      : ', self%number_of_edges
        print *, 'Number of triangles  : ', self%numberOfTriangles

        do n = 1,self%numberOfTriangles
            call self%triangles(n)%print()
        end do

    end subroutine
    !===============================================================================================

    !===============================================================================================
    subroutine destroy(self)

        class(lagrangian_solid_2D), intent(inout) :: self

        deallocate(self%mass_points, self%edges, self%triangles,                     &
                   self%edgeMassPointsIndex, self%triangleEdgesIndex,                &
                   self%duplicateEdges, self%edgeTrianglesIndex, self%edgeDuplicate, &
                   self%oppositeVertexIndex)
        call self%center_of_mass%destroy()

    end subroutine
    !===============================================================================================

    !===============================================================================================
    subroutine update_forcing_elements(self)

        ! Update velocity and accelecartion of the forcin element from mass poitns

        use euclidean_mod, only : distance

        ! In/Out variables
        class(lagrangian_solid_2D), intent(inout), target :: self

        ! ! Local variables
        integer                 :: n
        real(dp)                :: w1, w2, w3, wtot
        type(triangle), pointer :: t 

        ! Cycle over forcing elements
        do n = 1,self%numberOfTriangles
            t => self%triangles(n)
            w1 = distance(t%v1%X, t%c%X)
            w2 = distance(t%v2%X, t%c%X)
            w3 = distance(t%v3%X, t%c%X)
            wtot = w1 + w2 + w3
            t%C%V = (t%v1%V*w1 + t%v2%V*w2 + t%v3%V*w3)/wtot
            t%C%A = (t%v1%A*w1 + t%v2%A*w2 + t%v3%A*w3)/wtot
        end do

    end subroutine
    !===============================================================================================

end module
