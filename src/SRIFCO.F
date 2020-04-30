      SUBROUTINE SRIFCO(T, LDT, AINV, LDA, AINVB, LDB, RINV, LDR, C,
     *                  LDC, QINV, LDQ, X, RINVY, W, N, M, P, WRK,
     *                  LDW, MULTRC, WITHX, TOL)
C
C     PURPOSE:
C
C     The algorithm calculates a combined measurement and time update
C     of one iteration of the time-invariant Kalman filter. This update
C     is given for the square root information filter, using the
C     condensed controller-Hessenberg form.
C
C     CONTRIBUTORS:
C
C     M. Vanbegin, P. Van Dooren (PRLB)
C     M. Verhaegen (NASA Ames)
C     I. Pribec (TUM)
C
C     REVISIONS:
C
C     1988, Sept. 9.
C     2020, Apr. 30.: Replaced DTRCO (LINPACK) with DTRCON (LAPACK)
C
C     Specification of parameters.
C
C     .. Scalar Arguments ..
C
      INTEGER LDT, LDA, LDB, LDR, LDC, LDQ, N, M, P, LDW
      DOUBLE PRECISION TOL
      LOGICAL MULTRC,WITHX
C
C     .. Array Arguments ..
C
      DOUBLE PRECISION T(LDT,*), AINV(LDA,*), AINVB(LDB,*),
     *                 RINV(LDR,*), C(LDC,*), QINV(LDQ,*), X(*),
     *                 RINVY(*), W(*), WRK(LDW,*)
C
C     EXTERNAL SUBROUTINES:
C
C     DTRCON from LAPACK
C     DTRMV, DTRSV from Extended-BLAS,
C     F06FBF, F06FSF, F06FUF from NAG-BLAS,
C     DAXPY, DCOPY from BLAS.
C
C     Local variables.
C
      INTEGER I, I1, IN, J, MI, MI1, MINMP, MN1, MP1, MNP
      DOUBLE PRECISION DZ1, RCOND
      INTEGER IWORK(N), INFO
      DOUBLE PRECISION WORK(3*N)
C
C     Construction of the pre-array WRK.
C
      MN1 = M + N + 1
      MP1 = M + P + 1
      MNP = M + N + P
      DO 20 J = 1, MN1
         CALL F06FBF(MNP, 0.0D+0, WRK(1,J), 1)
   20 CONTINUE
C
C     First part - Storing QINV in the (1,1) block of WRK.
C
      DO 40 J = 1, M
         CALL DCOPY(J, QINV(1,J), 1, WRK(1,J), 1)
   40 CONTINUE
C
C     Second part - Storing the process noise mean value
C                   in the (1,3) block of WRK.
C
      CALL DCOPY(M, W, 1, WRK(1,M+N+1), 1)
      CALL DTRMV('U', 'N', 'N', M, QINV, LDQ, WRK(1,M+N+1), 1)
C
C     Third part - Storing RINV x C in the (2,2) block of WRK.
C
      IF (MULTRC) THEN
         DO 80 J = 1, N
            CALL DCOPY(P, C(1,J), 1, WRK(M+1,M+J), 1)
   80    CONTINUE
      ELSE
         DO 100 I = 1, N
            CALL DCOPY(P, C(1,I), 1, WRK(M+1,M+I), 1)
            CALL DTRMV('U', 'N', 'N', P, RINV, LDR, WRK(M+1,M+I), 1)
  100    CONTINUE
      END IF
C
C     Fourth part - Storing the measurement in the (2,3) block of WRK.
C
      CALL DCOPY(P, RINVY, 1, WRK(M+1,M+N+1), 1)
C
C     Fifth part - Storing T x A and T x A x B in the (3,1) and
C                  (3,2) blocks of WRK.
*
      DO 140 I = 1, M
         DO 120 J = 1, MIN(I,N)
            CALL DAXPY(J, AINVB(J,I), T(1,J), 1, WRK(MP1,I), 1)
  120    CONTINUE
  140 CONTINUE
      DO 180 I = 1, N
         DO 160 J = 1, MIN(M+I,N)
            CALL DAXPY(J, AINV(J,I), T(1,J), 1, WRK(MP1,M+I), 1)
  160    CONTINUE
  180 CONTINUE
C
C     Sixth part - Storing T x X in the (3,3) block of WRK.
C
      CALL DCOPY(N, X, 1, WRK(MP1,M+N+1), 1)
      CALL DTRMV('U', 'N', 'N', N, T, LDT, WRK(MP1,M+N+1), 1)
C
C     Triangularization (2 steps).
C
C     Step 1: eliminate the (3,1) block  of WRK.
C
      DO 240 I = 1, M
         I1 = I + 1
         IN = MIN(I,N)
         CALL F06FSF(IN, WRK(I,I), WRK(MP1,I), 1, TOL, DZ1)
         DO 220 J = I1, MN1
            CALL F06FUF(IN, WRK(MP1,I), 1, DZ1, WRK(I,J), WRK(MP1,J), 1)
  220    CONTINUE
  240 CONTINUE
C
C     Step 2: triangularize the remaining (2,2) and (3,2) blocks of WRK.
C
      DO 280 I = 1, N
         MINMP = MIN(M+P,P+N-I)
         MI = M + I
         MI1 = MI + 1
         CALL F06FSF(MINMP, WRK(MI,MI), WRK(MI1,MI), 1, TOL, DZ1)
         DO 260 J = MI1, MN1
            CALL F06FUF(MINMP, WRK(MI1,MI), 1, DZ1, WRK(MI,J),
     *                WRK(MI1,J), 1)
  260    CONTINUE
C         CALL DTRCO(WRK(M+1,M+1), LDW, N, RCOND, WRK(M+1,1), 1)
         CALL DTRCON('1','U','N',N,WRK(M+1,M+1),LDW,RCOND,WORK,
     *               IWORK,INFO)
         IF (RCOND .LT. TOL) WITHX = .FALSE.
  280 CONTINUE
C
C     Output T and X.
C
      DO 300 J = 1, N
         CALL DCOPY(J, WRK(M+1,M+J), 1, T(1,J), 1)
  300 CONTINUE
C
      IF (WITHX) THEN
         CALL DCOPY(N, WRK(M+1,M+N+1), 1, X, 1)
         CALL DTRSV('U', 'N', 'N', N, T, LDT, X, 1)
      END IF
C
      RETURN
C
C *** Last line of the SRIFCO subroutine ******************************
C
      END