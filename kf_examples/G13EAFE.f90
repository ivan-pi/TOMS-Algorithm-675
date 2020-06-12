    Program g13eafe

!     G13EAF Example Program Text

!     Mark 27.0 Release. NAG Copyright 2019.

!     .. Use Statements ..
      Use nag_library, Only: daxpy, ddot, dgemv, dpotrf, dtrmv, dtrsv, g13eaf, &
                             nag_wp, x04caf
!     .. Implicit None Statement ..
      Implicit None
!     .. Parameters ..
      Real (Kind=nag_wp), Parameter    :: one = 1.0_nag_wp
      Real (Kind=nag_wp), Parameter    :: zero = 0.0_nag_wp
      Integer, Parameter               :: inc1 = 1, nin = 5, nout = 6
!     .. Local Scalars ..
      Real (Kind=nag_wp)               :: dev, tol
      Integer                          :: i, ifail, info, istep, l, ldm, ldq,  &
                                          lds, lwk, m, n, ncall, tdq
      Logical                          :: full, is_const, read_matrix, stq
!     .. Local Arrays ..
      Real (Kind=nag_wp), Allocatable  :: a(:,:), ax(:), b(:,:), c(:,:),       &
                                          h(:,:), k(:,:), p(:,:), q(:,:),      &
                                          r(:,:), s(:,:), wk(:), x(:), y(:),   &
                                          ymean(:)
      Integer, Allocatable             :: iwk(:)
!     .. Intrinsic Procedures ..
      Intrinsic                        :: log
!     .. Executable Statements ..
      Write (nout,*) 'G13EAF Example Program Results'
      Write (nout,*)

!     Skip heading in data file
      Read (nin,*)

!     Read in the problem size
      Read (nin,*) n, m, l, stq, is_const

      lds = n
      If (.Not. stq) Then
        ldq = l
        tdq = l
      Else
        ldq = 1
        tdq = 1
      End If
      ldm = m
      lwk = (n+m)*(n+m+l)
      Allocate (a(lds,n),b(lds,l),q(ldq,tdq),c(ldm,n),r(ldm,m),s(lds,n),       &
        k(lds,m),h(ldm,m),iwk(m),wk(lwk),x(n),ymean(m),y(m),ax(n),p(lds,n))

!     Read in the state covariance matrix, S
      Read (nin,*)(s(i,1:n),i=1,n)

!     Read in flag indicating whether S is the full matrix, or its
!     Cholesky decomposition
      Read (nin,*) full

!     If required (full), perform the Cholesky decomposition on S
      If (full) Then
!       The NAG name equivalent of dpotrf is f07fdf
        Call dpotrf('L',n,s,lds,info)
        If (info>0) Then
          Write (nout,*) ' S not positive definite'
          Go To 100
        End If
      End If

!     Read in initial state vector
      Read (nin,*) x(1:n)

!     Read in mean of the series
      Read (nin,*) ymean(1:m)

!     Read in control parameter
      Read (nin,*) ncall, tol

!     Display titles
      Write (nout,*) '         Residuals'
      Write (nout,*)

!     Initialize variables
      dev = zero
      read_matrix = .True.

!     Loop through data
      Do istep = 1, ncall
!       Read in the various matrices. If the series is constant
!       then this only happens at the first call
        If (read_matrix) Then

!         Read in transition matrix, A
          Read (nin,*)(a(i,1:n),i=1,n)
!         Read in noise coefficient matrix, B
          Read (nin,*)(b(i,1:l),i=1,n)
!         Read in measurement coefficient matrix, C
          Read (nin,*)(c(i,1:n),i=1,m)

!         Read in measurement noise covariance matrix, R
          Read (nin,*)(r(i,1:m),i=1,m)
!         Read in flag indicating whether R is the full matrix, or its
!         Cholesky decomposition
          Read (nin,*) full
!         If required (full), perform the Cholesky decomposition on R
          If (full) Then
!           The NAG name equivalent of dpotrf is f07fdf
            Call dpotrf('L',m,r,ldm,info)
            If (info>0) Then
              Write (nout,*) ' R not positive definite'
              Go To 100
            End If
          End If

!         Read in state noise matrix Q, if not assume to be identity matrix
          If (.Not. stq) Then
            Read (nin,*)(q(i,1:l),i=1,l)
!           Read in flag indicating whether Q is the full matrix, or its
!           Cholesky decomposition
            Read (nin,*) full
!           Perform Cholesky factorization on Q, if full matrix is supplied
            If (full) Then
!             The NAG name equivalent of dpotrf is f07fdf
              Call dpotrf('L',l,q,ldq,info)
              If (info>0) Then
                Write (nout,*) ' Q not positive definite'
                Go To 100
              End If
            End If
          End If

!         If series is constant set flag to false
          read_matrix = .Not. is_const
        End If

!       Read in observed values
        Read (nin,*) y(1:m)

!       Call G13EAF
        ifail = 0
        Call g13eaf(n,m,l,a,lds,b,stq,q,ldq,c,ldm,r,s,k,h,tol,iwk,wk,ifail)

!       Subtract the mean y:= y-ymean
!       The NAG name equivalent of daxpy is f06ecf
        Call daxpy(m,-one,ymean,inc1,y,inc1)

!       Perform time and measurement update x <= Ax + K(y-Cx)
!       The NAG name equivalent of dgemv is f06paf
        Call dgemv('N',m,n,-one,c,ldm,x,inc1,one,y,1)
        Call dgemv('N',n,n,one,a,lds,x,inc1,zero,ax,1)
        Call dgemv('N',n,m,one,k,lds,y,inc1,one,ax,1)
        x(1:n) = ax(1:n)

!       Display the residuals
        Write (nout,99999) y(1:m)

!       Update log-likelihood.
!       The NAG name equivalent of dtrsv is f06pjf
        Call dtrsv('L','N','N',m,h,ldm,y,1)
!       The NAG name equivalent of ddot is f06eaf
        dev = dev + ddot(m,y,1,y,1)
        Do i = 1, m
          dev = dev + 2.0_nag_wp*log(h(i,i))
        End Do
      End Do

!     Compute P from S
!     The NAG name equivalent of dtrmv is f06pff
      Do i = 1, n
        p(1:i,i) = s(i,1:i)
        Call dtrmv('L','N','N',i,s,lds,p(1,i),inc1)
        p(i,1:i-1) = p(1:i-1,i)
      End Do

!     Display final results
      Write (nout,*)
      Write (nout,*) ' Final X(I+1:I) '
      Write (nout,*)
      Write (nout,99999) x(1:n)
      Write (nout,*)
      Flush (nout)
      ifail = 0
      Call x04caf('Lower','Non-Diag',n,n,p,lds,'Final Value of P',ifail)
      Write (nout,*)
      Write (nout,99998) ' Deviance = ', dev

100   Continue

99999 Format (6F12.4)
99998 Format (A,E13.4)
    End Program g13eafe