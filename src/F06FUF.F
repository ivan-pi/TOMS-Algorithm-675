      SUBROUTINE F06FUF( N, Z, INCZ, Z1, ALPHA, X, INCX )
*     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA, Z1
      INTEGER            INCX, INCZ, N
*     .. Array Arguments ..
      DOUBLE PRECISION   X( * ), Z( * )
*     ..
*
*  F06FUF performs a Householder reflection given by
*
*     ( alpha ) = P*( alpha ) ,
*     (   x   )     (   x   )
*
*  where the orthogonal matrix p is given in the form
*
*     P = I - ( 1/z( 1 ) )*z*z'.
*
*  z( 1 ) must be supplied in Z1 and the remaining n elements in Z.
*  If Z1 is zero then P is assumed to be the unit matrix and the
*  transformation is skipped, otherwise Z1 must be in the range
*  ( 1.0, 2.0 ). Z1 and Z will usually be supplied by routine F06FSF.
*
*
*  Nag Fortran 77 O( n ) basic linear algebra routine.
*
*  -- Written on 2-November-1982.
*     Sven Hammarling, Nag Central Office.
*
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
*     .. Local Scalars ..
      DOUBLE PRECISION   BETA
*     .. External Functions ..
      DOUBLE PRECISION   DDOT
      EXTERNAL           DDOT
*     .. External Subroutines ..
      EXTERNAL           DAXPY
*     ..
*     .. Executable Statements ..
      IF( Z1.NE.ZERO )THEN
         BETA  = ALPHA*Z1 + DDOT( N, Z, INCZ, X, INCX )
         ALPHA = ALPHA    - BETA
         CALL DAXPY( N, -BETA/Z1, Z, INCZ, X, INCX )
      END IF
*
      RETURN
*
*     End of F06FUF. ( SREF )
*
      END