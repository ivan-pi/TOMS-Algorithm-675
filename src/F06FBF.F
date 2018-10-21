      SUBROUTINE F06FBF( N, CONST, X, INCX )
*     .. Scalar Arguments ..
      DOUBLE PRECISION   CONST
      INTEGER            INCX, N
*     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
*     ..
*
*  F06FBF performs the operation
*
*     x = const*e,   e' = ( 1  1 ... 1 ).
*
*
*  Nag Fortran 77 O( n ) basic linear algebra routine.
*
*  -- Written on 22-September-1983.
*     Sven Hammarling, Nag Central Office.
*
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
*     .. Local Scalars ..
      INTEGER            IX
*     ..
*     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( CONST.NE.ZERO )THEN
            DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = CONST
   10       CONTINUE
         ELSE
            DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = ZERO
   20       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of F06FBF. ( SLOAD )
*
      END