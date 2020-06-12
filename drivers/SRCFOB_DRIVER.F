C     This is the driver for SUBROUTINE SRCFOB.
C     The routine is in TOMS.FOR and uses the file SRCFOB.IN as
C     input file from which N, P, A, B, C, R and Q are read.
C     Output is sent to file SRCFOB.OUT .
C
      INTEGER I, ISTEP, J, N, M, P, LDS, LDA, LDB, LDQ, LDC, LDR, LDK,
     *        LDW, LDU
      DOUBLE PRECISION A(10,10), B(10,5), C(7,10), Q(5,5), R(7,7),
     *                 WRK(22,22), S(10,10), UB(10,5),
     *                 K(10,7), U(10,10), SSP(10,10),
     *                 UPSSP(10,10), UPSSPU(10,10)
      DOUBLE PRECISION DDOT, TOL
      LOGICAL WITHU, UPPER, MULTBQ, WITHK
C
      WITHK = .TRUE.
      UPPER = .FALSE.
      WITHU = .TRUE.
      MULTBQ = .TRUE.
      LDS = 10
      LDA = 10
      LDB = 10
      LDC = 7
      LDQ = 5
      LDR = 7
      LDK = 10
      LDW = 22
      LDU = 10
      TOL = 1.D-15
C
C     READ N, P, A, B, C, Q and R from the file SRCFOB.IN .
C     The matrices A, B, C, Q and R are to be read columnwise.
C     The input dimension M is specified by the loop index.
C
C
      OPEN(UNIT=1,FILE='SRCFOB.IN',ACTION='READ')
      OPEN(UNIT=2,FILE='SRCFOB.OUT',ACTION='WRITE')
      DO 750 M = 2, 3
         READ (1,51) N, P
  51  FORMAT(2I5)
  53  FORMAT(1H ,'*** N =', I3,' M =', I3,' P = ', I3)
         IF (M .EQ. 2) THEN
            WRITE(2,940)
         ELSE
            WRITE(2,941)
         END IF
         DO 50 J = 1, N
            DO 50 I = 1, N
               READ(1,52) A(I,J)
  50     CONTINUE
  52  FORMAT(D25.15)
         DO 60 J = 1, M
            DO 60 I = 1, N
               READ(1,52) B(I,J)
  60     CONTINUE
         DO 68 J = 1, N
            DO 68 I = 1, P
               READ(1,52) C(I,J)
  68     CONTINUE
         DO 76 J = 1, P
            DO 76 I = 1, P
               READ(1,52) R(I,J)
  76     CONTINUE
         DO 84 J = 1, M
            DO 84 I = 1, M
               READ(1,52) Q(I,J)
  84     CONTINUE
         DO 90 J = 1, N
            DO 90 I = 1, N
               U(I,J) = 0.0D0
   90    CONTINUE
         DO 100 I = 1, N
            U(I,I) = 1.0D0
  100    CONTINUE
         WRITE(2,950)
         WRITE(2,947)
         WRITE(2,53) N, M, P
         WRITE(2,944)
         CALL PRMT(A,LDA,N,N,'A matrix',2,4)
         CALL PRMT(B,LDB,N,M,'B matrix',2,4)
         CALL PRMT(C,LDC,P,N,'C matrix',2,4)
C
C     Transform (A,B,C) to lower observer Hessenberg form.
C
         CALL OBHESS(A, LDA, N, C, LDC, P, U, LDU, WITHU, UPPER)
         DO 140 J = 1, M
            DO 140 I = 1, N
               UB(I,J) = DDOT(N, U(I,1), LDU, B(1,J), 1)
  140    CONTINUE
         WRITE(2,948)
         CALL PRMT(A,LDA,N,N,'Au matr.',2,4)
         CALL PRMT(B,LDB,N,M,'Bu matr.',2,4)
         CALL PRMT(C,LDC,P,N,'Cu matr.',2,4)
         CALL PRMT(U,LDU,N,N,'U matr. ',2,4)
         WRITE(2,949)
         CALL PRMT(R,LDR,P,P,'R matr. ',2,4)
         CALL PRMT(Q,LDQ,M,M,'Q matr. ',2,4)
C
C     Now perform three steps of the Kalman filter recursion
C     (in square root covariance form) with MULTBQ = .TRUE.
C
         WRITE(2,951)
C
C     Initialize the S matrix to be the null matrix.
C
         DO 350 J = 1, N
            DO 350 I = 1, N
               S(I,J) = 0.0D0
  350    CONTINUE
         DO 500 ISTEP = 1, 3
            WRITE(2,935)ISTEP
            CALL SRCFOB(S, LDS, A, LDA, UB, LDU, Q, LDQ, C, LDC, R,
     *                  LDR, N, M, P, K, LDK, WRK, LDW, MULTBQ,
     *                  WITHK, TOL)
            CALL PRMT(S,LDS,N,N,'Su matr.',2,4)
            CALL PRMT(K,LDK,N,P,'Ku matr.',2,4)
            DO 450 J = 1, N
               DO 450 I = 1, N
                  SSP(I,J) = DDOT(N, S(I,1), LDS, S(J,1), LDS)
  450       CONTINUE
            DO 460 J = 1, N
               DO 460 I = 1, N
                  UPSSP(I,J) = DDOT(N, U(1,I), 1, SSP(1,J), 1)
  460       CONTINUE
            DO 470 J = 1, N
               DO 470 I = 1, N
                  UPSSPU(I,J) = DDOT(N, UPSSP(I,1), LDU, U(1,J), 1)
  470       CONTINUE
            CALL PRMT(UPSSPU,LDU,N,N,'U''SuSu''U',2,4)
  500    CONTINUE
         WRITE(2,952)
C
C     Initialize the S matrix to be the null matrix.
C
         DO 550 J = 1, N
            DO 550 I = 1, N
               S(I,J) = 0.0D0
  550    CONTINUE
C
C     Now perform three steps of the Kalman filter recursion
C     (in square root covariance form) with MULTBQ = .FALSE.
C
         MULTBQ = .FALSE.
         DO 700 ISTEP = 1, 3
            WRITE(2,935)ISTEP
            CALL SRCFOB(S, LDS, A, LDA, UB, LDU, Q, LDQ, C, LDC, R,
     *                  LDR, N, M, P, K, LDK, WRK, LDW, MULTBQ,
     *                  WITHK, TOL)
            CALL PRMT(S,LDS,N,N,'Su matr.',2,4)
            CALL PRMT(K,LDK,N,P,'Ku matr.',2,4)
            DO 650 J = 1, N
               DO 650 I = 1, N
                  SSP(I,J) = DDOT(N, S(I,1), LDS, S(J,1), LDS)
  650       CONTINUE
            DO 660 J = 1, N
               DO 660 I = 1, N
                  UPSSP(I,J) = DDOT(N, U(1,I), 1, SSP(1,J), 1)
  660       CONTINUE
            DO 670 J = 1, N
               DO 670 I = 1, N
                  UPSSPU(I,J) = DDOT(N, UPSSP(I,1), LDU, U(1,J), 1)
  670       CONTINUE
            CALL PRMT(UPSSPU,LDU,N,N,'U''SuSu''U',2,4)
  700    CONTINUE
         IF (M .EQ. 2) THEN
            WRITE(2,945)
         ELSE
            WRITE(2,946)
         END IF
  750 CONTINUE

      CLOSE(1)
      CLOSE(2)

  935 FORMAT(' *** ISTEP =', I3)
  940 FORMAT(
     * '1*** First example : Square root covariance filter with '/
     * ' *** A, B, C in (lower) observer Hessenberg form and'/
     * ' *** lower triangular Q, R.')
  941 FORMAT(
     * '1*** Second example : Square root covariance filter with'/
     * ' *** A, B, C in (lower) observer Hessenberg form and'/
     * ' *** lower triangular Q, R.')
  944 FORMAT(// ' *** UPPER = .FALSE.' //)
  945 FORMAT(
     * '1*** In both these tests we start with Su=0 and perform three'/
     * ' *** iterations of the filter.'/
     * ' *** The ranks of Su in these three steps must be equal to '/
     * ' *** M(=2), 2M(=4) and 3M(=6).'/
     * ' *** K must be 0 in the first step and nonzero afterwards.'/
     * ' *** The results for MULTBQ = .TRUE. and .FALSE. must be equal'/
     * ' *** since we chose Q=I.'/
     * ' *** The U''SuSu''U matrices are meant for comparison with '/
     * ' *** SRCF.')
  946 FORMAT(
     * '1*** In both these tests we start with Su=0 and perform three'/
     * ' *** iterations of the filter.'/
     * ' *** The ranks of Su in these three steps must be equal to '/
     * ' *** M(=3), 2M(=6) and min(N,3M)(=6).'/
     * ' *** K must be 0 in the first step and nonzero afterwards.'/
     * ' *** The results for MULTBQ = .TRUE. and .FALSE. must be equal'/
     * ' *** since we chose Q=I.'/
     * ' *** The U''SuSu''U matrices are meant for comparison with '/
     * ' *** SRCF.')
  947 FORMAT(// ' *** Input OBHESS' //)
  948 FORMAT(// ' *** Output OBHESS' //)
  949 FORMAT(// ' *** Additional input to SRCFOB' //)
  950 FORMAT(// ' *** Input SRCFOB' //)
  951 FORMAT(// ' *** Output SRCFOB for MULTBQ=.TRUE.' //)
  952 FORMAT(// ' *** Output SRCFOB for MULTBQ=.FALSE.' //)
C
C *** Last line of the program SRCFOB *********************************
C
      END