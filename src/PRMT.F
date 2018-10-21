      SUBROUTINE PRMT(A, LDA, M, N, TEXT, KW, L)
C
C     LIBRARY INDEX:
C
C     1.4. General input/output routines.
C
C     PURPOSE:
C
C     This routine prints out the M by N matrix A row by row.
C     The elements of A are printed out with FORMAT(D15.7).
C     The routine first prints the contents of TEXT as a title
C     and next the elements of the matrix A in the following way.
C     - if N <= L, the M x L block is printed.
C     - if N = k L + p, k > 0, then k  M x L blocks of
C       consecutive columns of A are printed one after the other
C       followed by the M x p block of the p last columns of A.
C     Row numbers are printed on the left of each row and a column
C     number on top of each column.
C     If M <= 0 or N <= 0 or L <= 0 then the subroutine call is
C     an empty statement.
C     The routine uses 2 + (k + 1) x (m + 1) lines and 5 + c x 15
C     positions on each line where c is the actual number of columns,
C     (i.e., c = L or c = p).
C
C     CONTRIBUTOR:
C
C     H. Willemsen, Eindhoven University of Technology.
C
C     REVISIONS:
C
C     1987, November 23.
C
C     .. Scalar arguments ..
C
      INTEGER KW, L, M, LDA, N
      CHARACTER*8 TEXT
C
C     .. Array arguments ..
C
      DOUBLE PRECISION A(LDA,N)
C
C     Local variables:
C
      INTEGER I, J, J1, J2, JJ, N1
*
      WRITE (KW, FMT=999) TEXT, M, N
      IF (M.LE.0 .OR. N.LE.0 .OR. L.LE.0) RETURN
      IF (L.GE.9) L = 8
      N1 = (N - 1)/L
      J1 = 1
      J2 = L
      DO 20 J = 1, N1
         WRITE (KW, FMT=996) (JJ, JJ=J1,J2)
         DO 10 I = 1, M
            WRITE (KW, FMT=998) I, (A(I,JJ), JJ=J1,J2)
   10    CONTINUE
         WRITE (KW, FMT=997)
         J1 = J1 + L
         J2 = J2 + L
   20 CONTINUE
      WRITE (KW, FMT=996) (J, J=J1,N)
      DO 30 I = 1, M
         WRITE (KW, FMT=998) I, (A(I,JJ), JJ=J1,N)
   30 CONTINUE
      WRITE (KW, FMT=997)
*
  996 FORMAT (5X, 8(6X, I2, 7X))
  997 FORMAT (1H )
  998 FORMAT (X, I2, 2X, 8D15.7)
  999 FORMAT (1X, A8, 2H (, I2, 1HX, I2, 1H), /)
*
      RETURN
C *** Last line of PRMT *****************************************
      END