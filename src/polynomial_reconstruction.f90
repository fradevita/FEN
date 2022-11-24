!> This module conains procedures for reconstructing any field locally with a 
!> polynomial expression of a given order.
module polynomial_reconstruction

    use precision, only : dp

    implicit none

    ! Matrix coefficient for the reconstruction of the surface polynomial.
    ! The System is:
    ! [XYZ] [a_p] = [b]
    ! where [XYZ] is the Vandermonde matrix, a the coefficients to be determined
    ! and b the interpolated field. Setting A = [XYZ] and x = [a_p] the system is
    ! transformed as
    ! AT*A*x = AT*b
    ! x = inv_(AT*A) * AT * b
    ! a_p is the solution array x
    real(dp), dimension(:)  , allocatable :: ap
    real(dp), dimension(:,:), allocatable :: inv_ATAxAT
 
    ! The matrix A has size mxn which depends on the polynomial order
    integer :: m, n

    private
    public :: init
    public :: get_coefficients
    public :: P2nd
    public :: destroy

contains

    !=====================================================================================
    subroutine init(order)

        use class_Grid, only : base_grid

        ! In/Out variables
        integer :: order

        ! Local variables
        integer :: info, row, j, i, lwork
        integer , dimension(:)  , allocatable :: ipiv
        real(dp) :: xc, yc, work_dimension(1)
        real(dp), dimension(:)  , allocatable :: work
        real(dp), dimension(:,:), allocatable :: A, AT, inv_ATA


        ! Select matrix dimension based on polynomial order
        if (order == 2) then
            m = 9
            n = 6
        else
            print *, 'ERROR: for now only 2nd order polynomial are allowd'
            stop
        endif

        allocate(A(m,n))
        allocate(AT(n,m))
        allocate(inv_ATA(n,n))
        allocate(inv_ATAxAT(n,m))
        allocate(ipiv(n))
        allocate(ap(n))
  
        ! Build the matrix A
        row = 1
        do j = -1,1
           yc = base_grid%delta*real(j, dp)
           do i = -1,1
              xc = base_grid%delta*real(i, dp)
              !          [ a00, a10, a20, a01, a11, a02]
              A(row,:) = [1.0_dp, xc, xc**2, yc, xc*yc, yc**2]
              row = row + 1
           end do
        end do
  
        ! Compute tranpose of A
        AT = transpose(A)
  
        ! Compute the matrix AT*A for now stored in inv_ATA
        inv_ATA = matmul(AT, A)
  
        ! Perform LU factorization
        call dgetrf(n, n, inv_ATA, n, ipiv, info)
        if (info /= 0) print *, 'error in dgertf'
  
        ! Query for optimal dimension
        call dgetri(n, inv_ATA, n, ipiv, work_dimension, -1, info)
        if (info /= 0) print *, 'error in dgerti'
  
        ! Allocate optimal work space
        lwork = int(work_dimension(1))
        allocate(work(lwork))
  
        ! Compute invers of AT*A
        call dgetri(n, inv_ATA, n, ipiv, work, lwork, info)
        if (info /= 0) print *, 'error in dgerti'
  
        ! Compute inv_ATA*AT
        inv_ATAxAT = matmul(inv_ATA, AT)
  
        ! Free memory
        deallocate(A, AT, inv_ATA, ipiv, work)

    end subroutine
    !=====================================================================================

    !=====================================================================================
    function get_coefficients(f) result(c)
 
        ! In/Out variables
        real(dp), intent(in ) :: f(-1:1,-1:1)
        real(dp)              :: c(n)
  
        c = matmul(inv_ATAxAT, reshape(f, (/m/)))
  
    end function get_coefficients
    !=====================================================================================
  
    !=====================================================================================
     function P2nd(x, a) result(f)
  
        real(dp), intent(in) :: x(:), a(:)
        real(dp)             :: f
  
        f = a(1) + a(2)*x(1) + a(3)*x(1)**2 + a(4)*x(2) + a(5)*x(1)*x(2) + a(6)*x(2)**2
  
    end function P2nd
    !=====================================================================================

    !=====================================================================================
    subroutine destroy()

        deallocate(ap, inv_ATAxAT)

    end subroutine destroy
    !=====================================================================================

end module polynomial_reconstruction