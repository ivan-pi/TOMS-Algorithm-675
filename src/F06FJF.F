      SUBROUTINE F06FJF( N, X, INCX, SCALE, SUMSQ )
*     .. Scalar Arguments ..
      DOUBLE PRECISION   SCALE, SUMSQ
      INTEGER            INCX, N
*     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
*     ..
*
*  F06FJF returns the values scl and smsq such that
*
*     ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
*
*  where y( i ) = X( 1 + ( i - 1 )*INCX ). The value of sumsq is assumed
*  to be at least unity and the value of smsq will then satisfy
*
*     1.0 .le. smsq .le. ( sumsq + n ) .
*
*  scale is assumed to be non-negative and scl returns the value
*
*     scl = max( scale, abs( x( i ) ) ) .
*
*  scale and sumsq must be supplied in SCALE and SUMSQ respectively.
*  scl and smsq are overwritten on SCALE and SUMSQ respectively.
*
*  The routine makes only one pass through the vector X.
*
*
*  Nag Fortran 77 O( n ) basic linear algebra routine.
*
*  -- Written on 22-October-1982.
*     Sven Hammarling, Nag Central Office.
*
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
*     .. Local Scalars ..
      DOUBLE PRECISION   ABSXI
      INTEGER            IX
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
      IF( N.GT.0 )THEN
         DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
            IF( X( IX ).NE.ZERO )THEN
               ABSXI = ABS( X( IX ) )
               IF( SCALE.LT.ABSXI )THEN
                  SUMSQ = 1     + SUMSQ*( SCALE/ABSXI )**2
                  SCALE = ABSXI
               ELSE
                  SUMSQ = SUMSQ +       ( ABSXI/SCALE )**2
               END IF
            END IF
   10    CONTINUE
      END IF
      RETURN
*
*     End of F06FJF. ( SSSQ )
*
      END