      SUBROUTINE SRIF(T, LDT, AINV, LDA, B, LDB, RINV, LDR, C, LDC,
     *                QINV, LDQ, X, RINVY, W, N, M, P, WRK, LDW,
     *                MULTAB, MULTRC, WITHX, TOL)
C
C     PURPOSE:
C
C     The algorithm calculates a combined measurement and time update
C     of one iteration of the Kalman filter. This update is given for
C     the square root information filter, using dense matrices.
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
      LOGICAL MULTAB, MULTRC, WITHX
C
C     .. Array Arguments ..
C
      DOUBLE PRECISION T(LDT,*), AINV(LDA,*), B(LDB,*), RINV(LDR,*),
     *                 C(LDC,*), QINV(LDQ,*), X(*), RINVY(*), W(*),
     *                 WRK(LDW,*)
C
C     EXTERNAL SUBROUTINES:
C
C     DTRCON from LAPACK
C     DGEMV, DTRMV, DTRSV from Extended-BLAS,
C     F06FBF, F06FSF, F06FUF from NAG-BLAS,
C     DCOPY from BLAS
C
C     Local variables.
C
      INTEGER I, I1, J, J1, MN1, MNP
      DOUBLE PRECISION DZ1, RCOND
      INTEGER IWORK(N), INFO
      DOUBLE PRECISION WORK(3*N)
C
C     Construction of the pre-array WRK.
C
      MN1 = M + N + 1
      MNP = M + N + P
      DO 20 J = 1, MN1
         CALL F06FBF(MNP, 0.0D+0, WRK(1,J), 1)
   20 CONTINUE
C
C     First part - Storing QINV in the (1,1) block of WRK.
C
C
      DO 40 J = 1, M
         CALL DCOPY(J, QINV(1,J), 1, WRK(1,J), 1)
   40 CONTINUE
C
C     Second part - Storing the process noise mean value in the (1,3)
C                   block of WRK.
C
      CALL DCOPY(M, W, 1, WRK(1,M+N+1), 1)
      CALL DTRMV('U', 'N', 'N', M, QINV, LDQ, WRK(1,M+N+1), 1)
C
C     Third part - Storing T x AINV and T x AINV x B in the
C                  (2,1) and (2,2) blocks of WRK.
C
      DO 80 I = 1, N
         CALL DCOPY(N, AINV(1,I), 1, WRK(M+1,M+I), 1)
         CALL DTRMV('U', 'N', 'N', N, T, LDT, WRK(M+1,M+I), 1)
   80 CONTINUE
      IF (MULTAB) THEN
         DO 100 I = 1, M
            CALL DCOPY(N, B(1,I), 1, WRK(M+1,I), 1)
            CALL DTRMV('U', 'N', 'N', N, T, LDT, WRK(M+1,I), 1)
  100    CONTINUE
      ELSE
         DO 120 I = 1, M
            CALL DGEMV('N', N, N, 1.0D0, WRK(M+1,M+1), LDW, B(1,I), 1,
     *                 0.0D0, WRK(M+1,I), 1)
  120    CONTINUE
      END IF
C
C     Fourth part - Storing T x X in the (2,3) block of WRK.
C
      CALL DCOPY(N, X, 1, WRK(M+1,M+N+1), 1)
      CALL DTRMV('U', 'N', 'N', N, T, LDT, WRK(M+1,M+N+1), 1)
C
C     Fifth part - Storing RINV x C in the (3,2) block of WRK.
C
      IF (MULTRC) THEN
         DO 160 J = 1, N
            CALL DCOPY(P, C(1,J), 1, WRK(M+N+1,M+J), 1)
  160    CONTINUE
      ELSE
         DO 180 I = 1, N
            CALL DCOPY(P, C(1,I), 1, WRK(M+N+1,M+I), 1)
            CALL DTRMV('U', 'N', 'N', P, RINV, LDR, WRK(M+N+1,M+I), 1)
  180    CONTINUE
      END IF
C
C     Sixth part - Storing the measurement in the (3,3) block of WRK.
C
      CALL DCOPY(P, RINVY, 1, WRK(M+N+1,M+N+1), 1)
C
C     Triangularization (2 steps).
C
C     Step 1: eliminate the (2,1) block of WRK.
C
      DO 220 I = 1, M
         I1 = I + 1
         CALL F06FSF(N, WRK(I,I), WRK(M+1,I), 1, TOL, DZ1)
         DO 200 J = I1, MN1
            CALL F06FUF(N, WRK(M+1,I), 1, DZ1, WRK(I,J), WRK(M+1,J), 1)
  200    CONTINUE
  220 CONTINUE
C
C     Step 2: triangularize the remaining (2,2) and (3,2) blocks of WRK.
C
      DO 260 I = 1, N
         CALL F06FSF(N+P-I, WRK(M+I,M+I), WRK(M+I+1,M+I), 1, TOL, DZ1)
         J1 = N - I + 1
         DO 240 J = 1, J1
            CALL F06FUF(N+P-I, WRK(M+I+1,M+I), 1, DZ1, WRK(M+I,M+I+J),
     *                  WRK(M+I+1,M+I+J), 1)
  240    CONTINUE
C         CALL DTRCO(WRK(M+1,M+1), LDW, N, RCOND, WRK(M+1,1), 1)
         CALL DTRCON('1','U','N',N,WRK(M+1,M+1),LDW,RCOND,WORK,
     *               IWORK,INFO)
         IF (RCOND .LT. TOL) WITHX = .FALSE.
  260 CONTINUE
C
C     Output T and X.
C
      DO 280 J = 1, N
         CALL DCOPY(J, WRK(M+1,M+J), 1, T(1,J), 1)
  280 CONTINUE
C
      IF (WITHX) THEN
         CALL DCOPY(N, WRK(M+1,M+N+1), 1, X, 1)
         CALL DTRSV('U', 'N', 'N', N, T, LDT, X, 1)
      END IF
C
      RETURN
C
C *** Last line of the SRIF subroutine ********************************
C
      END