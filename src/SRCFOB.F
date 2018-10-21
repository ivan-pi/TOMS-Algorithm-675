      SUBROUTINE SRCFOB(S, LDS, A, LDA, B, LDB, Q, LDQ, C, LDC, R, LDR,
     *                 N, M, P, K, LDK, WRK, LDW, MULTBQ, WITHK, TOL)
C
C     PURPOSE:
C
C     The algorithm calculates a combined measurement and time update
C     of one iteration of the time-invariant Kalman filter. This update
C     is given for the square root covariance filter, using the
C     condensed observer-Hessenberg form.
C
C     CONTRIBUTORS:
C
C     M. Vanbegin, P. Van Dooren (PRLB)
C     M. Verhaegen (NASA Ames)
C
C     REVISIONS:
C
C     1988, Sept. 9.
C
C     Specification of parameters.
C
C     .. Scalar arguments ..
C
      INTEGER LDS, LDA, LDB, LDQ, LDC, LDR, N, M, P, LDK, LDW
      DOUBLE PRECISION TOL
      LOGICAL MULTBQ, WITHK
C
C     .. Array Arguments ..
C
      DOUBLE PRECISION S(LDS,*), A(LDA,*), B(LDB,*), Q(LDQ,*), C(LDC,*),
     *                 R(LDR,*), K(LDK,*), WRK(LDW,*)
C
C     EXTERNAL SUBROUTINES:
C
C     DTRCO from LINPACK
C     DGEMV, DTRSV from Extended-BLAS,
C     F06FBF, F06FSF, F06FUF from NAG-BLAS,
C     DAXPY, DCOPY from BLAS.
C
C     Local variables.
C
      INTEGER I, I1, IN, IND, IPM, J, MINPN, MINMP, P1, PM1, PNM, PN,
     *        PI, PI1, PJ, NAXPY
      DOUBLE PRECISION DZ1, RCOND
C
C     Construction of the pre-array WRK.
C
      MINPN = MIN(P,N)
      P1 = P + 1
      PM1 = P1 + M
      PN = P + N
      PNM = PN + M
      DO 20 J = 1, PNM
         CALL F06FBF(PN, 0.0D+0, WRK(1,J), 1)
   20 CONTINUE
C
C     First part -  Storing lower triangular factor R in the (1,1) block
C                   of WRK.
C
      DO 40 I = 1, P
         CALL DCOPY(P-I+1, R(I,I), 1, WRK(I,I), 1)
   40 CONTINUE
C
C     Second part - Storing B x Q in the (2,2) block of WRK.
C
      IF (MULTBQ) THEN
         DO 60 I = 1, M
            CALL DCOPY(N, B(1,I), 1, WRK(P1,P+I), 1)
   60    CONTINUE
      ELSE
         DO 80 I = 1, M
            CALL DGEMV('N', N, M-I+1, 1.0D0, B(1,I), LDB, Q(I,I), 1,
     *                 0.0D0, WRK(P1,P+I), 1)
   80    CONTINUE
      END IF
C
C     Third part -  Storing C x S in the (1,3) block of WRK.
C
      DO 120 I = 1, MINPN
         IPM = I + P + M
         NAXPY = P - I + 1
         DO 100 J = I, MINPN
            CALL DAXPY(NAXPY, S(J,I), C(J,J), 1, WRK(J,IPM), 1)
            NAXPY = NAXPY - 1
  100    CONTINUE
  120 CONTINUE
C
C     Fourth part  -  Storing A x S in the (2,3) block of WRK.
C
      DO 160 I = 1, N
         IPM = I + P + M
         NAXPY = N
         IND = 1
         DO 140 J = I, N
            IF (J .GT. P1) THEN
               IND = IND + 1
               NAXPY = NAXPY - 1
            END IF
            CALL DAXPY(NAXPY, S(J,I), A(IND,J), 1, WRK(IND+P,IPM), 1)
  140    CONTINUE
  160 CONTINUE
C
C     Triangularization (2 steps).
C
C     Step 1: eliminate the (1,3) block WRK.
C
      DO 200 I = 1, P
         IN = MIN(I,N)
         CALL F06FSF(IN, WRK(I,I), WRK(I,PM1), LDW, TOL, DZ1)
         I1 = I + 1
         DO 180 J = I1, PN
            CALL F06FUF(IN, WRK(I,PM1), LDW, DZ1, WRK(J,I), WRK(J,PM1),
     *                  LDW)
  180    CONTINUE
         CALL DTRCO(WRK, LDW, P, RCOND, WRK(1,PNM), 0)
         IF (RCOND .LT. TOL) WITHK = .FALSE.
  200 CONTINUE
C
C     Step 2: triangularize the remaining (2,2) and (2,3) blocks of WRK.
C
      DO 240 I = 1, N
         MINMP = MIN(M+P,M+N-I)
         PI = P + I
         PI1 = PI + 1
         CALL F06FSF(MINMP, WRK(PI,PI), WRK(PI,PI1), LDW, TOL, DZ1)
         IF (PI1 .LE. PN) THEN
            DO 220 J = PI1, PN
               CALL F06FUF(MINMP, WRK(PI,PI1), LDW, DZ1, WRK(J,PI),
     *                     WRK(J,PI1), LDW)
  220       CONTINUE
         END IF
  240 CONTINUE
C
C     Output K and S.
C
      IF (WITHK) THEN
         DO 260 J = 1, N
            CALL DCOPY(P, WRK(P+J,1), LDW, K(J,1), LDK)
            CALL DTRSV('L', 'T', 'N', P, WRK, LDW, K(J,1), LDK)
  260    CONTINUE
      END IF
      DO 280 J = 1, N
         PJ = P + J
         CALL DCOPY(N-J+1, WRK(PJ,PJ), 1, S(J,J), 1)
  280 CONTINUE
      RETURN
C
C *** Last line of the SRCFOB subroutine ******************************
C
      END