      DOUBLE PRECISION FUNCTION F06BMF( SCALE, SSQ )
*     .. Scalar Arguments ..
      DOUBLE PRECISION                  SCALE, SSQ
*     ..
*
*  F06BMF returns the value norm given by
*
*     norm = scale*sqrt( ssq )
*
*  via the function name.
*
*
*  Nag Fortran 77 O( 1 ) basic linear algebra routine.
*
*  -- Written on 22-October-1982.
*     Sven Hammarling, Nag Central Office.
*     Modified by M. Vanbegin and P. Van Dooren, PRLB
*
*     .. Intrinsic Functions ..
      INTRINSIC             SQRT
*     ..
*     .. Executable Statements ..
*
      F06BMF = SCALE*SQRT( SSQ )
      RETURN
*
*     End of F06BMF. ( SNORM )
*
      END