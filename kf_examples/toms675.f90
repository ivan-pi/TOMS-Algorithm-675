module toms675

  implicit none
  private

  public :: srcf, srif
  public :: srcfob, srifco
  public :: cohess, obhess
  public :: dp
  public :: prmt
  
  public :: upper_triang_chol
  public :: lower_triang_chol
  public :: UDU_factorization
  public :: one_step_prediction

  integer, parameter :: dp = kind(1.0d0)

  interface

!>  Square Root Covariance Filter
    subroutine srcf(s, lds, a, lda, b, ldb, q, ldq, c, ldc, r, ldr,&
                    n, m, p, k, ldk, wrk, ldw, multbq, withk, tol)
      import dp
      integer, intent(in) :: lds, lda, ldb, ldq, ldc, ldr, n, m, p, ldk, ldw
      real(dp), intent(in) :: tol
      logical, intent(in) :: multbq
      logical, intent(inout) :: withk
      real(dp), intent(inout) :: s(lds,*), k(ldk,*), wrk(ldw,*)
      real(dp), intent(in) :: a(lda,*), b(ldb,*), q(ldq,*), c(ldc,*),r(ldr,*)
    end subroutine

!>  Square Root Information Filter
    subroutine srif(t, ldt, ainv, lda, b, ldb, rinv, ldr, c, ldc,&
                    qinv, ldq, x, rinvy, w, n, m, p, wrk, ldw,&
                    multab, multrc, withx, tol)
      import dp
      integer, intent(in) :: ldt,lda,ldb,ldr,ldc,ldq,n,m,p,ldw
      real(dp), intent(in) :: tol
      logical, intent(in) :: multab,multrc
      logical, intent(inout) :: withx
      real(dp), intent(inout) :: t(ldt,*), x(*), wrk(ldw,*)
      real(dp), intent(in) :: ainv(lda,*), b(ldb,*), rinv(ldr,*), c(ldc,*),&
                              qinv(ldq,*), rinvy(*), w(*)          
    end subroutine

!>  Square Root Information Filter using the condensed observer Hessenberg form
    subroutine srcfob(s, lds, a, lda, b, ldb, q, ldq, c, ldc, r, ldr,&
                      n, m, p, k, ldk, wrk, ldw, multbq, withk, tol)
      import dp
      integer, intent(in) :: lds, lda, ldb, ldq, ldc, ldr, n, m, p, ldk, ldw
      real(dp), intent(in) :: tol
      logical, intent(in) :: multbq
      logical, intent(inout) :: withk
      real(dp), intent(inout) :: s(lds,*), k(ldk,*), wrk(ldw,*)
      real(dp), intent(in) :: a(lda,*), b(ldb,*), q(ldq,*), c(ldc,*),r(ldr,*)
    end subroutine

!>  Square Root Information Filter using the condensed controller Hessenberg form
    subroutine srifco(t, ldt, ainv, lda, ainvb, ldb, rinv, ldr, c,&
                      ldc, qinv, ldq, x, rinvy, w, n, m, p, wrk,&
                      ldw, multrc, withx, tol)
      import dp
      integer, intent(in) :: ldt,lda,ldb,ldr,ldc,ldq,n,m,p,ldw
      real(dp), intent(in) :: tol
      logical, intent(in) :: multrc
      logical, intent(inout) :: withx
      real(dp), intent(inout) :: t(ldt,*), x(*), wrk(ldw,*)
      real(dp), intent(in) :: ainv(lda,*), ainvb(ldb,*), rinv(ldr,*),&
                              c(ldc,*),qinv(ldq,*), rinvy(*), w(*)
    end subroutine

    subroutine obhess(a,lda,n,c,ldc,p,u,ldu,withu,upper)
      import dp
      integer, intent(in) :: lda, n, ldc, p, ldu
      logical, intent(in) :: withu, upper
      real(dp), intent(inout) :: a(lda,*), c(ldc,*), u(ldu,*)
    end subroutine

    subroutine cohess(a,lda,n,b,ldb,m,u,ldu,withu,upper)
      import dp
      integer, intent(in) :: lda, n, ldb, m, ldu
      logical, intent(in) :: withu, upper
      real(dp), intent(inout) :: a(lda,*), b(ldb,*), u(ldu,*)
    end subroutine

    subroutine prmt(A, LDA, M, N, TEXT, KW, L)
      import dp
      INTEGER, intent(in) :: KW, L, M, LDA, N
      CHARACTER(len=8), intent(in) :: TEXT
      real(dp), intent(in) :: A(LDA,N)
    end subroutine
  end interface

contains

  subroutine one_step_prediction(S,lds,A,lda,B,ldb,Q,ldq,n,m,multbq,tol)
    integer, intent(in) :: lds, lda, ldb, ldq, n, m
    real(dp), intent(inout) :: S(lds,n)
      !! A two-dimensional real array containing the lower triangular square root
      !! of the n x n state covariance matrix P_{k|k}. Upon return it contains
      !! the corresponding square root of P_{k+1|k}.
    real(dp), intent(in) :: A(lda,n)
      !! A two-dimensional real array containing the n x n state transition matrix
      !! A_k of the discrete-time system.
    real(dp), intent(in) :: B(ldb,m) 
      !! A two-dimensional real array containing the n x m input weight matrix B_k
      !! (or its product with Q if `multbq = .true.`)
    real(dp), intent(in) :: Q(ldq,m)
      !! A two-dimensional real array containing the m * m lower triangular Cholesky
      !! factor of the process noise covariance matrix 
    logical, intent(in) :: multbq
      !! Logical variable indicating how input matrices B and Q are transferred.
      !! If `multbq = .true.`, then the product B*Q is transferred via B, and Q is not used.
      !! If `multbq = .false.`, then B and Q are transferred via the parameters B and Q,
      !! respectively.
    real(dp), intent(in) :: tol

    real(dp), allocatable :: wrk(:,:)
    real(dp) :: dz1
    integer :: nm, ldw, i, i1, j

    nm = n+m
    ldw = n
    allocate(wrk(n,nm))

!
!   First part -  Storing B x Q in the (1,2) block of WRK.
!
    if (multbq) then
      do i = 1, m
        call dcopy(n, b(1,i), 1, wrk(1,n+i), 1)
      end do
    else
      do i = 1, m
        call dgemv('n',n,m-i+1,1._dp,B(1,i),ldb,Q(i,i),1,0.0_dp,wrk(1,n+i),1)
      end do
    end if
    ! print *, "Step 1"
    ! do i = 1, n
    !   print *, wrk(i,:)
    ! end do
!
!   Second part  -  Storing A x S in the (1,1) block of WRK.
!
    do i = 1, n
      call dgemv('n',n,n-i+1,1._dp,A(1,i),lda,S(i,i),1,0._dp,wrk(1,i),1)
    end do
!     print *, "Step 2"
!     do i = 1, n
!       print *, wrk(i,:)
!     end do
! !
!   Step 3: triangularize the remaining (1,1) and (1,2) blocks of WRK.
!
    do i = 1, n
      i1 = i + 1
      call f06fsf(n+m-i,wrk(i,i),wrk(i,i1),ldw,tol,dz1)
      if (i1 <= n) then
        do j = i1, n
          call f06fuf(n+m-i,wrk(i,i1),ldw,dz1,wrk(j,i), &
                      wrk(j,i1),ldw)
          end do
      end if
    end do
    ! print *, "Step 3"
    ! do i = 1, n
    !   print *, wrk(i,:)
    ! end do

!
!   Output S.
!
    do i = 1, n
      call dcopy(n-i+1,wrk(i,i),1,S(i,i),1)
    end do

  end subroutine

  subroutine upper_triang_chol(n,P)
    integer, intent(in) :: n
    real(dp), intent(inout) :: P(n,n)
    integer :: j,k,i
    real(dp) :: alpha, beta

    do j = n, 2, -1
      p(j,j) = sqrt(p(j,j))
      alpha = 1._dp/p(j,j)
      do k = 1, j-1
        p(k,j) = alpha*p(k,j)
        beta = p(k,j)
        do i = 1, k
          p(i,k) = p(i,k) - beta*p(i,j)
        end do
      end do
    end do

    p(1,1) = sqrt(p(1,1))
  end subroutine

!> Lower triangular Cholesky factorization
  subroutine lower_triang_chol(n,P)
    integer, intent(in) :: n
      !! Size of matrix P
    real(dp), intent(inout) :: P(n,n)
      !! On input a symmetric positive definite matrix
      !! On output the lower part of P contains the lower
      !! Cholesky factorization L, where P = L L^T
    integer :: j,k,i
    real(dp) :: alpha, beta

    do j = 1, n-1
      p(j,j) = sqrt(p(j,j))
      alpha = 1._dp/p(j,j)
      do k = n, j+1,-1
        beta = p(k,j)
        p(k,j) = alpha*beta
        do i = k, n
          p(i,k) = p(i,k) - p(i,j)*beta
        end do
      end do
    end do
    p(n,n) = sqrt(p(n,n))
  end subroutine

!> Upper triangular UDU^T factorization
  subroutine UDU_factorization(n,P)
    integer, intent(in) :: n
      !! Dimension of the matrix
    real(dp), intent(inout) :: P(n,n)
      !! On input a symmetric matrix (positive definite ?)
      !! On output the P = UDU^T factorization, where the elements
      !! of D are stored on the diagonal, and the upper portion of P
      !! contains the values of U (the diagonal values of U are implicitly 1)

    real(dp) :: alpha, beta
    integer :: i,j,k

    do j = n, 2, -1
      ! d(j) = p(j,j)
      alpha = 1._dp/p(j,j)
      do k = 1, j-1
        beta = p(k,j)
        p(k,j) = alpha*beta
        do i = 1, k
          p(i,k) = p(i,k) - beta*p(i,j)
        end do
      end do
    end do
    ! d(1) = p(1,1)

  end subroutine
end module

! program test_prediction

!   use toms675, only: dp, one_step_prediction
!   implicit none

!   real(dp) :: a(2,2), p(2,2), b(2,2), q(2,2), sol(2,2), c(2,2), r(2,2)
!   logical :: multbq
!   real(dp) :: tol
!   integer :: info, i

!   tol = epsilon(tol)

!   p = 0
!   p(1,1) = 1
!   p(2,2) = 1
!   do i = 1, 2
!     print *, p(i,:)
!   end do

!   a = reshape([1.0_dp,0.0_dp,0.01_dp,1.0_dp],[2,2])
!   b = reshape([1,0,0,1],[2,2])
!   q = reshape([0.001_dp,0.0_dp,0.0_dp,0.0001_dp],[2,2])

!   sol = matmul(a,matmul(p,transpose(a))) + matmul(b,matmul(q,transpose(b)))

!   call dpotrf('L',2,sol,2,info)
!   if (info > 0) then
!     write(*,*) "sol is not positive definite"
!     error stop
!   end if
!   print *, "sol = "
!   do i = 1, 2
!     print *, sol(i,:)
!   end do
!   sol(1,2) = 0.0_dp
!   sol = matmul(sol,transpose(sol))
!   do i = 1, 2
!     print *, sol(i,:)
!   end do

!   call dpotrf('L',2,p,2,info)
!   if (info > 0) then
!     write(*,*) "P is not positive definite"
!     error stop
!   end if
!   call dpotrf('L',2,q,2,info)
!   if (info > 0) then
!     write(*,*) "q is not positive definite"
!     error stop
!   end if
!   multbq = .false.
!   call one_step_prediction(p,2,a,2,b,2,q,2,2,2,multbq,tol)
!   do i = 1, 2
!     print *, p(i,:)
!   end do

!   p = matmul(p,transpose(p))
!   do i = 1, 2
!     print *, p(i,:)
!   end do
! end program





! program test_chol

!   use toms675, only: dp, lower_triang_chol, upper_triang_chol, &
!     UDU_factorization
!   implicit none

!   real(dp) :: p(3,3), u(3,3), d(3,3)
!   real(dp) :: vals(9)

!   integer :: i, info

!   vals = [1._dp,0.5_dp,0.1_dp,0.5_dp,1._dp,0.5_dp,0.1_dp,0.5_dp,1._dp]

!   p = reshape(vals,[3,3])

!   write(*,*)
!   do i = 1, 3
!     print *, p(i,:)
!   end do

!   call upper_triang_chol(3,p)

!   write(*,*)
!   do i = 1, 3
!     print *, p(i,:)
!   end do

!   p = reshape(vals,[3,3])
!   ! call lower_triang_chol(3,p)
!   call dpotrf('L',3,p,3,info)
!   if (info > 0) then
!     write(*,*) "P is not positive definite"
!     error stop
!   end if


!   write(*,*)
!   do i = 1, 3
!     print *, p(i,:)
!   end do

!   p = reshape(vals,[3,3])
!   call UDU_factorization(3,p)

!   u = 0
!   u(1,1) = 1
!   u(2,2) = 1
!   u(3,3) = 1
!   u(1,2) = p(1,2)
!   u(1,3) = p(1,3)
!   u(2,3) = p(2,3)

!   d = 0
!   d(1,1) = p(1,1)
!   d(2,2) = p(2,2)
!   d(3,3) = p(3,3)

!   write(*,*)
!   do i = 1, 3
!     print *, p(i,:)
!   end do

!   p = matmul(u,matmul(d,transpose(u)))
!   write(*,*)
!   do i = 1, 3
!     print *, p(i,:)
!   end do
! end program