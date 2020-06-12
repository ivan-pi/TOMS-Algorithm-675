      SUBROUTINE F06PJF( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
*     .. Entry Points ..
      ENTRY      DTRSV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
*     .. Scalar Arguments ..
      INTEGER            INCX, LDA, N
      CHARACTER*1        DIAG, TRANS, UPLO
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * )
*     ..
*
*  Purpose
*  =======
*
*  F06PJF solves one of the systems of equations
*
*     A*x = b,   or   A'*x = b,
*
*  where b and x are n element vectors and A is an n by n unit, or
*  non-unit, upper or lower triangular matrix.
*
*  No test for singularity or near-singularity is included in this
*  routine. Such tests must be performed before calling this routine.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the equations to be solved as
*           follows:
*
*              TRANS = 'N' or 'n'   A*x = b.
*
*              TRANS = 'T' or 't'   A'*x = b.
*
*              TRANS = 'C' or 'c'   A'*x = b.
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit
*           triangular as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry with  UPLO = 'U' or 'u', the leading n by n
*           upper triangular part of the array A must contain the upper
*           triangular matrix and the strictly lower triangular part of
*           A is not referenced.
*           Before entry with UPLO = 'L' or 'l', the leading n by n
*           lower triangular part of the array A must contain the lower
*           triangular matrix and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U' or 'u', the diagonal elements of
*           A are not referenced either, but are assumed to be unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least n.
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element right-hand side vector b. On exit, X is overwritten
*           with the solution vector x.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X.
*           Unchanged on exit.
*
*
*  Note that UPLO, TRANS, DIAG, N and LDA must be such that the value of
*  the LOGICAL variable OK in the following statement is true.
*
*     OK = ( ( UPLO.EQ.'U' ).OR.( UPLO.EQ.'u' ).OR.
*    $       ( UPLO.EQ.'L' ).OR.( UPLO.EQ.'l' )     )
*    $     .AND.
*    $     ( ( TRANS.EQ.'N' ).OR.( TRANS.EQ.'n' ).OR.
*    $       ( TRANS.EQ.'T' ).OR.( TRANS.EQ.'t' ).OR.
*    $       ( TRANS.EQ.'C' ).OR.( TRANS.EQ.'c' )     )
*    $     .AND.
*    $     ( ( DIAG.EQ.'U' ).OR.( DIAG.EQ.'u' ).OR.
*    $       ( DIAG.EQ.'N' ).OR.( DIAG.EQ.'n' )     )
*    $     .AND.
*    $     ( N.GE.0 )
*    $     .AND.
*    $     ( LDA.GE.N )
*
*
*
*  Level 2 Blas routine.
*
*  -- Written on 30-September-1985.
*     Sven Hammarling, Nag Central Office.
*
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
*     .. Local Scalars ..
      INTEGER            I, IX, J, JX, KX
      LOGICAL            NOUNIT
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible.
*
      IF( N.EQ.0 )
     $   RETURN
      NOUNIT = ( DIAG.EQ.'N' ).OR.( DIAG.EQ.'n' )
*
*     Set up the start point in X if the increment is not unity. This
*     will be  ( N - 1 )*INCX  too small for descending loops.
*
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX.NE.1 )THEN
         KX = 1
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF( ( TRANS.EQ.'N' ).OR.( TRANS.EQ.'n' ) )THEN
*
*        Form  x := inv( A )*x.
*
         IF( ( UPLO.EQ.'U' ).OR.( UPLO.EQ.'u' ) )THEN
            IF( INCX.EQ.1 )THEN
               DO 20, J = N, 1, -1
                  IF( X( J ).NE.ZERO )THEN
                     IF( NOUNIT )
     $                  X( J ) = X( J )/A( J, J )
                     DO 10, I = J - 1, 1, -1
                        X( I ) = X( I ) - X( J )*A( I, J )
   10                CONTINUE
                  END IF
   20          CONTINUE
            ELSE
               JX = KX + ( N - 1 )*INCX
               DO 40, J = N, 1, -1
                  IF( X( JX ).NE.ZERO )THEN
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )/A( J, J )
                     IX = JX
                     DO 30, I = J - 1, 1, -1
                        IX      = IX      - INCX
                        X( IX ) = X( IX ) - X( JX )*A( I, J )
   30                CONTINUE
                  END IF
                  JX = JX - INCX
   40          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 60, J = 1, N
                  IF( X( J ).NE.ZERO )THEN
                     IF( NOUNIT )
     $                  X( J ) = X( J )/A( J, J )
                     DO 50, I = J + 1, N
                        X( I ) = X( I ) - X( J )*A( I, J )
   50                CONTINUE
                  END IF
   60          CONTINUE
            ELSE
               JX = KX
               DO 80, J = 1, N
                  IF( X( JX ).NE.ZERO )THEN
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )/A( J, J )
                     IX = JX
                     DO 70, I = J + 1, N
                        IX      = IX      + INCX
                        X( IX ) = X( IX ) - X( JX )*A( I, J )
   70                CONTINUE
                  END IF
                  JX = JX + INCX
   80          CONTINUE
            END IF
         END IF
      ELSE
*
*        Form  x := inv( A' )*x.
*
         IF( ( UPLO.EQ.'U' ).OR.( UPLO.EQ.'u' ) )THEN
            IF( INCX.EQ.1 )THEN
               DO 100, J = 1, N
                  DO 90, I = 1, J - 1
                     X( J ) = X( J ) - A( I, J )*X( I )
   90             CONTINUE
                  IF( NOUNIT )
     $               X( J ) = X( J )/A( J, J )
  100          CONTINUE
            ELSE
               JX = KX
               DO 120, J = 1, N
                  IX = KX
                  DO 110, I = 1, J - 1
                     X( JX ) = X( JX ) - A( I, J )*X( IX )
                     IX      = IX      + INCX
  110             CONTINUE
                  IF( NOUNIT )
     $               X( JX ) = X( JX )/A( J, J )
                  JX = JX + INCX
  120          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 140, J = N, 1, -1
                  DO 130, I = N, J + 1, -1
                     X( J ) = X( J ) - A( I, J )*X( I )
  130             CONTINUE
                  IF( NOUNIT )
     $               X( J ) = X( J )/A( J, J )
  140          CONTINUE
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO 160, J = N, 1, -1
                  IX = KX
                  DO 150, I = N, J + 1, -1
                     X( JX ) = X( JX ) - A( I, J )*X( IX )
                     IX      = IX      - INCX
  150             CONTINUE
                  IF( NOUNIT )
     $               X( JX ) = X( JX )/A( J, J )
                  JX = JX - INCX
  160          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of F06PJF. ( DTRSV )
*
      END