C
C     This is the driver for SUBROUTINE SRIF.
C     The routine is in TOMS.FOR and uses the file SRIF.IN as
C     input file from which N, M, Ainv, C, AinvB, Qinv, Rinv, Z, X
C     and Y are read. Output is sent to file SRIF.OUT .
C
      INTEGER I, ISTEP, J, N, M, P, LDS, LDA, LDB, LDQ, LDC, LDR, LDW,
     *        NMP
      DOUBLE PRECISION AINV(10,10), C(5,10), RINV(5,5), QINV(7,7)
      DOUBLE PRECISION WRK(22,22), SINV(10,10), Y(5)
      DOUBLE PRECISION SPS(10,10), Z(7), X(10)
      DOUBLE PRECISION AINVB(10,7), RINVY(5)
      DOUBLE PRECISION DDOT, TOL
      LOGICAL MULTAB, MULTRC, WITHX
C
      WITHX = .TRUE.
      MULTAB = .TRUE.
      MULTRC = .FALSE.
      LDS = 10
      LDA = 10
      LDB = 10
      LDQ = 7
      LDC = 5
      LDR = 5
      LDW = 22
      TOL = 1.D-15
C
C     READ N, M, Ainv, C, AinvB, Qinv, Rinv, Z, X and Y from the file
C     SRIF.IN. The matrices Ainv, C, AinvB, Qinv and Rinv are to be
C     read columnwise.
C     The input dimension P is specified by the loop index.
C

      OPEN(UNIT=1,FILE='SRIF.IN',ACTION='READ')
      OPEN(UNIT=2,FILE='SRIF.OUT',ACTION='WRITE')

      DO 750 P = 2, 3
         READ (1,51) N, M
         NMP = N + M + P
         IF (LDW .LT. NMP) WRITE(6,901)
  51  FORMAT(2I5)
  53  FORMAT(1H ,'*** N =', I3,' M =', I3,' P = ', I3)
         IF (P .EQ. 2) THEN
            WRITE(2,940)
         ELSE
            WRITE(2,941)
         END IF
         DO 50 J = 1, N
            DO 50 I = 1, N
               READ(1,52) AINV(J,I)
  50     CONTINUE
  52  FORMAT(D25.15)
         DO 60 J = 1, P
            DO 60 I = 1, N
               READ(1,52) C(J,I)
  60     CONTINUE
         DO 68 I = 1, M
            DO 68 J = 1, N
               READ(1,52) AINVB(J,I)
  68     CONTINUE
         DO 76 J = 1, M
            DO 76 I = 1, M
               READ(1,52) QINV(J,I)
  76     CONTINUE
         DO 81 J = 1, P
            DO 81 I = 1, P
               READ(1,52) RINV(J,I)
  81     CONTINUE
         DO 86 J = 1, M
            READ(1,52) Z(J)
  86     CONTINUE
         DO 88 J = 1, N
            READ(1,52) X(J)
  88     CONTINUE
         DO 90 J = 1, P
            READ(1,52) Y(J)
  90     CONTINUE
C
C     Initialize the SINV matrix to be the unit matrix.
C
         DO 100 J = 1, N
            DO 100 I = 1, N
               SINV(I,J) = 0.0D0
  100    CONTINUE
         DO 110 I = 1, N
            SINV(I,I) = 1.0D0
  110    CONTINUE
         WRITE(2,950)
         WRITE(2,53) N, M, P
         CALL PRMT(AINV,LDA,N,N,'Ainv    ',2,4)
         CALL PRMT(AINVB,LDA,N,M,'AinvB   ',2,4)
         CALL PRMT(C,LDC,P,N,'C matrix',2,4)
         CALL PRMT(SINV,LDS,N,N,'Sinv    ',2,4)
         CALL PRMT(QINV,LDQ,M,M,'Qinv    ',2,4)
         CALL PRMT(RINV,LDR,P,P,'Rinv    ',2,4)
         CALL PRMT(Z,LDQ,M,1,'Z vector',2,4)
         CALL PRMT(X,LDA,N,1,'X vector',2,4)
         CALL PRMT(Y,LDC,P,1,'Y vector',2,4)
C
C     Calculation of RINV x Y.
C
         DO 260 I = 1, P
            RINVY(I) = DDOT(P, RINV(I,1), LDR, Y, 1)
  260    CONTINUE
C
C     Now perform four steps of the Kalman filter recursion
C     (in square root covariance form).
C
         WRITE(2,951)
         DO 500 ISTEP = 1, 4
            WRITE(2,935)ISTEP
            CALL SRIF(SINV, LDS, AINV, LDA, AINVB, LDB, RINV, LDR, C,
     *                LDC, QINV, LDQ, X, RINVY, Z, N, M, P, WRK, LDW,
     *                MULTAB, MULTRC, WITHX, TOL)
            CALL PRMT(SINV,LDS,N,N,'Sinv    ',2,4)
            CALL PRMT(X,LDA,N,1,'X vector',2,4)
            DO 450 J = 1, N
               DO 450 I = 1, N
                  SPS(I,J) = DDOT(N, SINV(1,I), 1, SINV(1,J), 1)
  450       CONTINUE
            CALL PRMT(SPS,LDS,N,N,"Sin'Sin ",2,4)
  500    CONTINUE
         IF (P .EQ. 2) THEN
            WRITE(2,945)
         ELSE
            WRITE(2,946)
         END IF
  750 CONTINUE

      CLOSE(1)
      CLOSE(2)

  900 FORMAT(' ', 8(D12.6,1X), D12.6)
  901 FORMAT(' *********** Dimensions of WRK are less than N+M+P ')
  935 FORMAT(' *** ISTEP =', I3)
  940 FORMAT(
     * '1*** First example : Square root information filter with'/
     * ' *** dense A, B, C and upper triangular Q, R. ***')
  941 FORMAT(
     * '1*** Second example : Square root information filter with'/
     * ' *** dense A, B, C and upper triangular Q, R. ***')
  945 FORMAT(
     * '1*** In both these tests we start with Sinv=I and perform'/
     * ' *** four iterations of the filter.'/
     * ' *** The Sinv''Sinv matrices and X vectors are meant for '/
     * ' *** comparison with SRIFCO.')
  946 FORMAT(
     * '1*** In both these tests we start with Sinv=I and perform'/
     * ' *** four iterations of the filter.'/
     * ' *** The Sinv''Sinv matrices and X vectors are meant for '/
     * ' *** comparison with SRIFCO.')
  950 FORMAT(// ' *** Input SRIF' //)
  951 FORMAT(// ' *** Output SRIF' //)
C
C *** Last line of the program SRIF ***********************************
C
      END