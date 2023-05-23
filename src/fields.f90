module fields_mod

    ! This module contains the fields (scalar, vector and tensor) procedures.

    use precision_mod, only : dp
    use scalar_mod
    use vector_mod
    use tensor_mod

    implicit none

    ! Define wrappers for procedures
    interface gradient
        module procedure gradient_of_scalar
        !module procedure gradient_of_vector
    end interface gradient
    
    interface divergence
        module procedure divergence_of_vector
        !module procedure divergence_of_tensor
    end interface divergence

    interface laplacian
        module procedure laplacian_of_scalar
        module procedure laplacian_of_vector
    end interface laplacian

contains

    !==============================================================================================
    subroutine gradient_of_scalar(s, grad_s)

        ! Compute the gradient grad_s of a scalar field s.

        ! In/Out variables
        type(scalar), intent(in   ) :: s      !< input scalar field
        type(vector), intent(inout) :: grad_s !< output vector field

        ! Local variables
        integer  :: i, j, jp, k
#if DIM==3
        integer :: kp
#endif
        real(dp) :: idelta

        idelta = 1.0_dp/s%G%delta

        do k = s%G%lo(3),s%G%hi(3)
#if DIM==3
            kp = k + 1
#endif
            do j = s%G%lo(2),s%G%hi(2)
                jp = j + 1
                    do i = s%G%lo(1),s%G%hi(1)
                        grad_s%x%f(i,j,k) = (s%f(i+1,j,k) - s%f(i,j,k))*idelta
                        grad_s%y%f(i,j,k) = (s%f(i,jp ,k) - s%f(i,j,k))*idelta
#if DIM==3
                        grad_s%z%f(i,j,k) = (s%f(i,j,kp) - s%f(i,j,k))*idelta
#endif
                end do
            end do
        end do

    end subroutine gradient_of_scalar
    !==============================================================================================

    !==============================================================================================
    subroutine gradient_of_vector(v, grad_v)

        ! Compute the gradient grad_v of a vector field v.
        ! The diagonal elements of the gradient are defined at cell center, the off diagonal
        ! at cell vertex.
        
        ! In/Out variables
        type(vector), intent(in   ) :: v
        type(tensor), intent(inout) :: grad_v
        
        ! Local variables
        integer  :: i, ip, im, j, jp, jm, k
#if DIM==3
        integer :: kp, km
#endif
        real(dp) :: idelta

        idelta = 1.0_dp/v%G%delta

        do k = v%G%lo(3),v%G%hi(3)
#if DIM==3
            kp = k + 1
            km = k - 1
#endif
            do j = v%G%lo(2),v%G%hi(2)
                jp = j + 1
                jm = j - 1
                do i = v%G%lo(1),v%G%hi(1)
                    ip = i + 1
                    im = i - 1
                    grad_v%x%x%f(i,j,k) = (v%x%f(i,j,k) - v%x%f(im,j,k))*idelta
                    grad_v%x%y%f(i,j,k) = 0.125_dp*(v%x%f(i ,jp,k) - v%x%f(i ,jm,k) + &
                                                    v%x%f(im,jp,k) - v%x%f(im,jm,k) + &
                                                    v%y%f(ip,j ,k) - v%y%f(im,j ,k) + &
                                                    v%y%f(ip,jm,k) - v%y%f(im,jm,k))*idelta
                    grad_v%y%x%f(i,j,k) = grad_v%x%y%f(i,j,k)
                    grad_v%y%y%f(i,j,k) = (v%y%f(i,j,k) - v%y%f(i,jm,k))*idelta             
#if DIM == 3

#endif
                end do
            end do
        end do

        ! Update boundary conditions
        if(grad_v%gl > 0) call grad_v%update_ghost_nodes()
        
    end subroutine gradient_of_vector
    !==============================================================================================

    !==============================================================================================
    subroutine divergence_of_vector(v, div_v)

        ! Compute the divergence div_v of a vector field v.

        ! In/Out variables
        type(vector), intent(in   ) :: v     !< input vector field
        type(scalar), intent(inout) :: div_v !< output scalar field

        ! Local variables
        integer :: i, j, jm, k
#if DIM==3
        integer :: km
#endif
        real(dp) :: idelta

        idelta = 1.0_dp/v%G%delta

        do k = v%G%lo(3),v%G%hi(3)
#if DIM==3
            km = k - 1
#endif
            do j = v%G%lo(2),v%G%hi(2)
                jm = j - 1
                do i = v%G%lo(1),v%G%hi(1)
                    div_v%f(i,j,k) = (v%x%f(i,j,k) - v%x%f(i-1,j,k))*idelta + &
                                     (v%y%f(i,j,k) - v%y%f(i,jm ,k))*idelta
#if DIM == 3
                    div_v%f(i,j,k) = div_v%f(i,j,k) + (v%z%f(i,j,k) - v%z%f(i,j,km))*idelta
#endif
                end do
            end do
        end do

    end subroutine divergence_of_vector
    !================================================================================================

    !========================================================================================
    subroutine divergence_of_tensor(T, div_T)

        ! Compute the divergence div_T of a tensor field T.

        ! In/Out variables
        type(tensor), intent(in   ) :: T
        type(vector), intent(inout) :: div_T

        call divergence(T%x, div_T%x) 
        call divergence(T%y, div_T%y)
#if DIM==3
        call divergence(T%z, div_T%z)
#endif

    end subroutine divergence_of_tensor
    !==============================================================================================

    !==============================================================================================
    subroutine center_to_face(s, v)

        ! Compute the face value v of the cell centered scalar s

        ! In/Out variables
        type(scalar), intent(in   ) :: s  !< input scalar field
        type(vector), intent(inout) :: v  !< input scalar face field (which is considered as vector)

        ! Local variables
        integer :: i, ip, j, jp, k
#if DIM==3
        integer :: kp
#endif
    
        do k = s%G%lo(3),s%G%hi(3)
#if DIM==3
            kp = k + 1
#endif
            do j = s%G%lo(2),s%G%hi(2)
                jp = j + 1
                do i = s%G%lo(1),s%G%hi(1)
                    ip = i + 1
                    v%x%f(i,j,k) = 0.5_dp*(s%f(ip,j,k) + s%f(i,j,k))
                    v%y%f(i,j,k) = 0.5_dp*(s%f(i,jp,k) + s%f(i,j,k))
#if DIM==3
                    v%z%f(i,j,k) = 0.5_dp*(s%f(i,j,kp) + s%f(i,j,k))
#endif
                end do
            end do
        end do

    end subroutine center_to_face
    !========================================================================================

    !==============================================================================================
    subroutine laplacian_of_scalar(s, lap_s)

        ! Compute the laplacian lap_v of the vector v.

        ! In/Out variables
        type(scalar), intent(in   ) :: s
        type(scalar), intent(inout) :: lap_s

        ! Local variables
        integer :: i, ip, im, j, jp, jm, k
#if DIM==3
        integer :: kp, km
#endif
        real(dp) :: idelta2

        idelta2 = 1.0_dp/s%G%delta**2 
    
        do k = s%G%lo(3),s%G%hi(3)
#if DIM==3
            kp = k + 1
            km = k - 1
#endif
            do j = s%G%lo(2),s%G%hi(2)
                jp = j + 1
                jm = j - 1
                do i = s%G%lo(1),s%G%hi(1)
                    ip = i + 1
                    im = i - 1
                    lap_s%f(i,j,k) = ((s%f(ip,j,k) - 2.0_dp*s%f(i,j,k) + s%f(im,j,k)) + &
                                      (s%f(i,jp,k) - 2.0_dp*s%f(i,j,k) + s%f(i,jm,k)))*idelta2
#if DIM==3
                    lap_s%f(i,j,k) = lap_s%f(i,j,k) + (s%f(i,j,kp) - 2.0_dp*s%f(i,j,k) + s%f(i,j,km))*idelta2
#endif
                end do
            end do
        end do

    end subroutine laplacian_of_scalar
    !================================================================================================

    !========================================================================================
    subroutine laplacian_of_vector(v, lap_v)

        ! Compute the laplacian lap_v of the vector v.

        ! In/Out variables
        type(vector), intent(in   ) :: v
        type(vector), intent(inout) :: lap_v

        ! Local variables
        integer :: i, ip, im, j, jp, jm, k
#if DIM==3
        integer :: kp, km
#endif
        real(dp) :: idelta2

        idelta2 = 1.0_dp/v%G%delta**2 
    
        do k = v%G%lo(3),v%G%hi(3)
#if DIM==3
            kp = k + 1
            km = k - 1
#endif
            do j = v%G%lo(2),v%G%hi(2)
                jp = j + 1
                jm = j - 1
                do i = v%G%lo(1),v%G%hi(1)
                    ip = i + 1
                    im = i - 1
                    lap_v%x%f(i,j,k) = ((v%x%f(ip,j,k) - 2.0_dp*v%x%f(i,j,k) + v%x%f(im,j,k)) + &
                                        (v%x%f(i,jp,k) - 2.0_dp*v%x%f(i,j,k) + v%x%f(i,jm,k)))*idelta2
                    lap_v%y%f(i,j,k) = ((v%y%f(ip,j,k) - 2.0_dp*v%y%f(i,j,k) + v%y%f(im,j,k)) + &
                                        (v%y%f(i,jp,k) - 2.0_dp*v%y%f(i,j,k) + v%y%f(i,jm,k)))*idelta2
#if DIM==3
                    lap_v%x%f(i,j,k) = lap_v%x%f(i,j,k) + &
                                        (v%x%f(i,j,kp) - 2.0_dp*v%x%f(i,j,k) + v%x%f(i,j,km))*idelta2
                    lap_v%y%f(i,j,k) = lap_v%y%f(i,j,k) + &
                                        (v%y%f(i,j,kp) - 2.0_dp*v%y%f(i,j,k) + v%y%f(i,j,km))*idelta2
                    lap_v%z%f(i,j,k) = ((v%z%f(ip,j,k) - 2.0_dp*v%z%f(i,j,k) + v%z%f(im,j,k)) + &
                                        (v%z%f(i,jp,k) - 2.0_dp*v%z%f(i,j,k) + v%z%f(i,jm,k)) + &
                                         (v%z%f(i,j,kp) - 2.0_dp*v%z%f(i,j,k) + v%z%f(i,j,km)))*idelta2
#endif
                end do
            end do
        end do

    end subroutine laplacian_of_vector
    !========================================================================================
  
    !========================================================================================
    subroutine curl(v, curl_v)

        ! Compute the curl curl_v of a vector field v.
        ! curl_v is defined on the top right cell vertex.
    
        ! In/Out variables
        type(vector), intent(in   ) :: v
        type(vector), intent(inout) :: curl_v

        ! Local variables
        integer  :: i, j, k, ip, jp
#if DIM==3
        integer :: kp
#endif
        real(dp) :: idelta

        idelta = 1.0_dp/v%G%delta

        do k = v%G%lo(3),v%G%hi(3)
#if DIM==3
            kp = k + 1
#endif
           do j = v%G%lo(2),v%G%hi(2)
                jp = j + 1
                do i = v%G%lo(1),v%G%hi(1)
                    ip = i + 1
#if DIM == 3
                    curl_v%x%f(i,j,k) = (v%z%f(i,jp,k) - v%z%f(i,j,k))*idelta - &
                                        (v%y%f(i,j,kp) - v%y%f(i,j,k))*idelta
                    curl_v%y%f(i,j,k) = (v%x%f(i,j,kp) - v%x%f(i,j,k))*idelta - &
                                        (v%z%f(ip,j,k) - v%z%f(i,j,k))*idelta
                    curl_v%z%f(i,j,k) = (v%y%f(ip,j,k) - v%y%f(i,j,k))*idelta - &
                                        (v%x%f(i,jp,k) - v%x%f(i,j,k))*idelta
#else
                    ! In 2D simulations the z component of vectors is not allocated so
                    ! we use the x component of curl_v.
                    curl_v%x%f(i,j,k) = (v%y%f(ip,j,k) - v%y%f(i,j,k))*idelta - &
                                        (v%x%f(i,jp,k) - v%x%f(i,j,k))*idelta
                                 
#endif
          end do
       end do
    end do

  end subroutine curl
  !========================================================================================

end module