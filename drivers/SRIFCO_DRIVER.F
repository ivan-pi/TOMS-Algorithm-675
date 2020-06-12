C     This is the driver for SUBROUTINE SRIFCO.
C     The routine is in TOMS.FOR and uses the file SRIFCO.IN as
C     input file from which N, M, Ainv, C, AinvB, Qinv, Rinv, Z, X,
C     and Y are read. Output is sent to file SRIFCO.OUT .
C
      INTEGER I, ISTEP, J, N, M, P, LDS, LDA, LDB, LDC, LDR, LDQ, LDW,
     *        LDU, NMP
      DOUBLE PRECISION AINV(10,10), AINVB(10,7), C(5,10), RINV(5,5),
     *                 AINVU(10,10), AINVBU(10,7), CU(5,10),
     *                 WRK(22,22), SINV(10,10), QINV(7,7), Y(5),
     *                 SPS(10,10), Z(7), X(10), XU(10),
     *                 RINVY(5), U(10,10),
     *                 SINVU(10,10)
      DOUBLE PRECISION DDOT, TOL
      LOGICAL MULTRC, WITHU, UPPER, WITHX
C
      WITHX = .TRUE.
      UPPER = .TRUE.
      MULTRC = .FALSE.
      WITHU = .TRUE.
      LDS = 10
      LDA = 10
      LDB = 10
      LDQ = 7
      LDC = 5
      LDR = 5
      LDW = 22
      LDU = 10
      TOL = 1.D-15
C
C     READ N, M, Ainv, C, AinvB, Qinv, Rinv, Z, X, and Y from the file
C     SRIFCO.IN.
C     The matrices Ainv, C, AinvB, Qinv and Rinv are to be read column-
C     wise.
C     The input dimension P is specified by the loop index.
C
C
      OPEN(UNIT=1,FILE='SRIFCO.IN',ACTION='READ')
      OPEN(UNIT=2,FILE='SRIFCO.OUT',ACTION='WRITE')
      DO 750 P = 2, 3
         READ (1,51) N, M
         NMP = N + M + P
         IF (LDW .LT. NMP) WRITE(6,901)
  51  FORMAT(2I5)
  53  FORMAT(1H , '*** N =', I3, ' M =', I3, ' P = ', I3)
         IF (P .EQ. 2) THEN
            WRITE(2,940)
         ELSE
            WRITE(2,941)
         END IF
         DO 50 J = 1, N
            DO 50 I =1, N
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
C     Initialize the U matrix to be the unity matrix.
C
         DO 100 J = 1, N
            DO 100 I = 1, N
               U(I,J) = 0.0D0
  100    CONTINUE
         DO 110 I = 1, N
            U(I,I) = 1.0D0
  110    CONTINUE
C
         WRITE(2,950)
         WRITE(2,947)
         WRITE(2,53) N, M, P
         WRITE(2,944)
         CALL PRMT(AINV,LDA,N,N,'Ainv    ',2,4)
         CALL PRMT(AINVB,LDA,N,M,'AinvB   ',2,4)
         CALL PRMT(C,LDC,P,N,'C matrix',2,4)
C
C     Transform (A,B,C) to upper Hessenberg form.
C
         CALL COHESS(AINV, LDA, N, AINVB, LDB, M, U, LDU, WITHU, UPPER)
         DO 150 I = 1, P
            CALL DGEMV('N', N, N, 1.0D0, U, LDU, C(I,1), LDC, 0.0D0,
     *                 CU(I,1), LDC)
  150    CONTINUE
         WRITE(2,948)
         CALL PRMT(AINV,LDA,N,N,'Ainvu   ',2,4)
         CALL PRMT(AINVB,LDA,N,M,'AinvBu   ',2,4)
         CALL PRMT(CU,LDC,P,N,'Cu matr.',2,4)
C
C     Initialize the SINVu matrix to be the unity matrix.
C
         DO 200 J = 1, N
            DO 190 I = 1, N
               SINVU(I,J) = 0.0D0
  190       CONTINUE
            SINVU(J,J) = 1.0D0
  200    CONTINUE
C
         CALL PRMT(SINVU,LDS,N,N,'Sinvu   ',2,4)
         WRITE(2,949)
         CALL PRMT(QINV,LDQ,M,M,'QINV    ',2,4)
         CALL PRMT(RINV,LDR,P,P,'RINV    ',2,4)
         CALL PRMT(Z,LDQ,M,1,'Z vector',2,4)
         CALL DGEMV('N', N, N, 1.0D0, U, LDU, X, 1, 0.0D0, XU, 1)
         CALL PRMT(XU,LDA,N,1,'XU vect.',2,4)
         CALL PRMT(Y,LDC,P,1,'Y vector',2,4)
C
C     Calculation of RINV x Y.
C
         DO 270 I = 1, P
            RINVY(I) = DDOT(P, RINV(I,1), LDR, Y, 1)
  270    CONTINUE
C
C     Now perform four steps of the Kalman filter recursion
C     (in square root covariance form).
C
         WRITE(2,951)
         DO 500 ISTEP = 1, 4
            WRITE(2,935)ISTEP
            CALL SRIFCO(SINVU, LDS, AINV, LDA, AINVB, LDB, RINV, LDR,
     *                  CU, LDC, QINV, LDQ, XU, RINVY, Z, N, M, P,
     *                  WRK, LDW, MULTRC, WITHX, TOL)
            CALL PRMT(SINVU,LDS,N,N,'Sinvu   ',2,4)
            DO 440 I = 1, N
               CALL DGEMV('N', N, N, 1.0D0, SINVU, LDS, U(1,I), 1,
     *                    0.0D0, SINV(1,I), 1)
  440       CONTINUE
            DO 460 J = 1, N
               DO 450 I = 1, N
                  SPS(I,J) = DDOT(N, SINV(1,I), 1, SINV(1,J), 1)
  450          CONTINUE
  460       CONTINUE
            CALL DGEMV('T', N, N, 1.0D0, U, LDU, XU, 1, 0.0D0, X, 1)
            CALL PRMT(X,LDA,N,1,'X vector',2,4)
            CALL PRMT(SPS,LDS,N,N,'Sin''Sin',2,4)
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
     * ' *** A, B, C in (upper) controller Hessenberg form'/
     * ' *** and upper triangular Q, R.')
  941 FORMAT(
     * '1*** Second example : Square root information filter with'/
     * ' *** A, B, C in (upper) controller Hessenberg form'/
     * ' *** and upper triangular Q, R.')
  944 FORMAT(// ' *** UPPER = .TRUE.'  //)
  945 FORMAT(
     * '1*** In both these tests we start with Sinvu=I and perform'/
     * ' *** four iterations of the filter.'/
     * ' *** The Sinv''Sinv matrices and X vectors are meant for '/
     * ' *** comparison with SRIF.')
  946 FORMAT(
     * '1*** In both these tests we start with Sinvu=I and perform'/
     * ' *** four iterations of the filter.'/
     * ' *** The Sinv''Sinv matrices and X vectors are meant for '/
     * ' *** comparison with SRIF.')
  947 FORMAT(// ' *** Input COHESS'  //)
  948 FORMAT(// ' *** Output COHESS'  //)
  949 FORMAT(// ' *** Additional input to SRIFCO'  //)
  950 FORMAT(// ' *** Input SRIFCO'  //)
  951 FORMAT(// ' *** Output SRIFCO MULTRC=.TRUE.'  //)
  952 FORMAT(// ' *** Output SRIFCO MULTRC=.FALSE.'  //)
C
C *** Last line of the program SRIFCO *********************************
C
      END