      SUBROUTINE F06FSF( N, ALPHA, X, INCX, TOL, Z1 )
*     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA, TOL, Z1
      INTEGER            INCX, N
*     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
*     ..
*
*  F06FSF generates details of a Householder reflection, P, such that
*
*     P*( alpha ) = ( beta ),   P'*P = I.
*       (   X   )   (   0  )
*
*  P is given in the form
*
*     P = I - ( 1/z( 1 ) )*z*z',
*
*  where z is an ( n + 1 ) element vector.
*
*  z( 1 ) is returned in Z1. If the elements of x are all zero, or if
*  the elements of x are all less than tol*abs( alpha ) in absolute
*  value, then Z1 is returned as zero and P can be taken to be the
*  unit matrix. Otherwise Z1 always lies in the range ( 1.0, 2.0 ).
*
*  If TOL is not in the range ( 0.0, 1.0 ) then the value 0.0 is used in
*  place of TOL.
*
*  The remaining elements of z are overwritten on X and beta is
*  overwritten on ALPHA.
*
*
*  Nag Fortran 77 O( n ) basic linear algebra routine.
*
*  -- Written on 2-November-1982.
*     Sven Hammarling, Nag Central Office.
*
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     .. Local Scalars ..
      DOUBLE PRECISION   BETA, SCALE, SSQ, TL
*     .. Local Arrays ..
      DOUBLE PRECISION   WORK( 1 )
*     .. External Functions ..
      DOUBLE PRECISION   F06BMF
      EXTERNAL           F06BMF
*     .. External Subroutines ..
      EXTERNAL           F06FJF, DSCAL
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN
*     ..
*     .. Executable Statements ..
      IF( N.LT.1 )THEN
         Z1 = ZERO
      ELSE
         IF( ( TOL.LE.ZERO ).OR.( TOL.GT.ONE ) )THEN
            TL = ZERO
         ELSE
            TL = ABS( ALPHA )*TOL
         END IF
         SSQ   = ONE
         SCALE = ZERO
         CALL F06FJF( N, X, INCX, SCALE, SSQ )
         IF( ( SCALE.EQ.ZERO ).OR.( SCALE.LT.TL ) )THEN
            Z1 = ZERO
         ELSE
            IF( ALPHA.NE.ZERO )THEN
               WORK( 1 ) = ALPHA
               CALL F06FJF( 1, WORK( 1 ), 1, SCALE, SSQ )
               BETA = -SIGN( F06BMF( SCALE, SSQ ), ALPHA )
               Z1   =  ( BETA - ALPHA )/BETA
            ELSE
               BETA = -F06BMF( SCALE, SSQ )
               Z1   =  ONE
            END IF
            CALL DSCAL( N, -ONE/BETA, X, INCX )
            ALPHA = BETA
         END IF
      END IF
*
      RETURN
*
*     End of F06FSF. ( SREFG )
*
      END