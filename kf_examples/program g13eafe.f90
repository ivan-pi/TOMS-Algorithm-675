program g13eafe

    use toms675, only: srcf, wp => dp, prmt
    implicit none

    Real(Kind=wp), Parameter    :: one = 1.0_wp
    Real(Kind=wp), Parameter    :: zero = 0.0_wp
    Integer, Parameter          :: inc1 = 1
!   .. Local Scalars ..
    Real(Kind=wp)               :: dev, tol
    Integer                     :: i, ifail, info, istep, l, ldm, ldq,  &
                                   lds, lwk, m, n, ncall, tdq
    Logical                     :: full, is_const, read_matrix, stq, withk
!   .. Local Arrays ..
    Real(Kind=wp), Allocatable  :: a(:,:), ax(:), b(:,:), c(:,:),       &
                                   h(:,:), k(:,:), p(:,:), q(:,:),      &
                                   r(:,:), s(:,:), wk(:,:), x(:), y(:),   &
                                   ymean(:), ky(:)

    Integer, Allocatable        :: iwk(:)
!   .. Intrinsic Procedures ..
    Intrinsic                   :: log
!   .. External procedures ..
    real(wp), external :: ddot

    write(*,*) "TOMS 675 Example Program Results"
    write(*,*)

!   Skip heading
    read(*,*) 

!   Read in the problem size
    read(*,*) n, m, l, stq, is_const
    write(*,*) n, m, l, stq, is_const

    lds = n
    if (.not. stq) then
        ldq = l
        tdq = l
    else
        ldq = 1
        tdq = 1
    end if
    ldm = m
    lwk = (n+m)

    Allocate(a(lds,n),b(lds,l),q(ldq,tdq),c(ldm,n),r(ldm,m),s(lds,n), &
             k(lds,m),h(ldm,m),iwk(m),wk(lwk,(n+m+l)),x(n),ymean(m),y(m),ax(n),ky(n),p(lds,n))

!   Read in the state covariance matrix, S
    read(*,*) (s(i,1:n),i=1,n)
    write(*,"(4(F6.4,:,X))") (s(i,1:n),i=1,n)

!   Read in flag indicating whether S is the full matrix, or its
!   Cholesky decomposition
    read(*,*) full
    write(*,*) full
    if (full) then
        call dpotrf('L',n,s,lds,info)
        if (info > 0) then
            write(*,*) ' S is not positive definite'
            error stop 1
        end if
    end if

!   Read in initial state vector
    read(*,*) x(1:n)
    write(*,*) x(1:n)

!   Read in mean of the series
    read(*,*) ymean(1:m)
    write(*,*) ymean(1:m)

!   Read in control parameter
    read(*,*) ncall, tol
    write(*,*) ncall, tol

!   Display titles
    write(*,*) '        Residuals'
    write(*,*)

!   Initialize variables
    dev = zero
    read_matrix = .true.

!   Loop through data
    Do istep = 1, ncall
!       Read in the various matrices. If the series is constant
!       then this only happens at the first call
        If (read_matrix) Then

!           Read in transition matrix, A
            Read (*,*)(a(i,1:n),i=1,n)
!           Read in noise coefficient matrix, B
            Read (*,*)(b(i,1:l),i=1,n)
!           Read in measurement coefficient matrix, C
            Read (*,*)(c(i,1:n),i=1,m)

!           Read in measurement noise covariance matrix, R
            Read (*,*)(r(i,1:m),i=1,m)
!           Read in flag indicating whether R is the full matrix, or its
!           Cholesky decomposition
            Read (*,*) full
!           If required (full), perform the Cholesky decomposition on R
            If (full) Then
!               The NAG name equivalent of dpotrf is f07fdf
                Call dpotrf('L',m,r,ldm,info)
                If (info>0) Then
                    Write (*,*) ' R not positive definite'
                    error stop 1
                End If
            End If

!           Read in state noise matrix Q, if not assume to be identity matrix
            If (.Not. stq) Then
                Read (*,*)(q(i,1:l),i=1,l)
!               Read in flag indicating whether Q is the full matrix, or its
!               Cholesky decomposition
                Read (*,*) full
!               Perform Cholesky factorization on Q, if full matrix is supplied
                If (full) Then
!                   The NAG name equivalent of dpotrf is f07fdf
                    Call dpotrf('L',l,q,ldq,info)
                    If (info>0) Then
                        Write (*,*) ' Q not positive definite'
                        error stop 1
                    End If
                End If
            End If

!           If series is constant set flag to false
            read_matrix = .Not. is_const
        End If

!       Read in observed values
        Read (*,*) y(1:m)
        write (*,*) "y = ", y(1:m)

!       Call SRCF
        withk = .true.
        call srcf(s,lds,a,lds,b,lds,q,ldq,c,ldm,r,ldm,n,l,m,k,lds,wk,lwk,.false.,withk,tol)

!       Subtract the mean y:= y-ymean
!       The NAG name equivalent of daxpy is f06ecf
        Call daxpy(m,-one,ymean,inc1,y,inc1)

!       Perform time and measurement update x <= Ax + K(y-Cx)
!       The NAG name equivalent of dgemv is f06paf
        Call dgemv('N',m,n,-one,c,ldm,x,inc1,one,y,1) ! y := -Cx + y
        Call dgemv('N',n,n,one,a,lds,x,inc1,zero,ax,1) ! ax := Ax
        ! call dgemv('N',n,m,one,k,lds,y,inc1,zero,ky,1) ! ky := Ky
        ! Call dgemv('N',n,n,one,a,lds,ky,inc1,one,ax,1) ! ax := Aky + ax
        Call dgemv('N',n,m,one,k,lds,y,inc1,one,ax,1) ! ax := Ky + ax
        x(1:n) = ax(1:n)

!       Display the residuals
        Write (*,*) y(1:m)

!       Extract H from the work array
        do i = 1, m
            h(1:ldm,i) = -wk(1:m,i) ! Why is minus needed here!?
        end do

        call prmt(h,ldm,m,m,'H chol',6,m)
        call prmt(wk,lwk,lwk,m+n+l,'work arr',6,m+n+l)


!       Update log-likelihood.
!       The NAG name equivalent of dtrsv is f06pjf
        Call dtrsv('L','N','N',m,h,ldm,y,1)
!       The NAG name equivalent of ddot is f06eaf
        dev = dev + ddot(m,y,1,y,1)
        Do i = 1, m
          dev = dev + 2.0_wp*log(h(i,i))
        End Do
    End Do

!   Compute P from S
!   The NAG name equivalent of dtrmv is f06pff
    Do i = 1, n
      p(1:i,i) = s(i,1:i)
      Call dtrmv('L','N','N',i,s,lds,p(1,i),inc1)
      p(i,1:i-1) = p(1:i-1,i)
    End Do

!   Display final results
    Write (*,*)
    Write (*,*) ' Final X(I+1:I) '
    Write (*,*)
    Write (*,"(6F12.4)") x(1:n)
    Write (*,*)

    print *, "p = "
    call prmt(p,lds,n,n,'fval of P',6,n)

    p = matmul(s,transpose(s))
    print *, "p = "
    call prmt(p,lds,n,n,'fval of P',6,n)

    Write (*,*)
    Write (*,'(A,E13.4)') ' Deviance = ', dev
end program