      SUBROUTINE SRCF_M(S, LDS, A, LDA, B, LDB, Q, LDQ, C, LDC, R, LDR,
     *                N, M, P, K, LDK, WRK, LDW, MULTBQ, WITHK, TOL)
C
C     PURPOSE:
C
C     The algorithm calculates a combined measurement and time update
C     of one iteration of the Kalman filter. This update is given for
C     the square root covariance filter, using dense matrices.
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
C     DTRCON from LAPACK
C     DGEMV, DTRSV from EXTENDED-BLAS,
C     F06FBF, F06FUF, F06FSF from NAG-BLAS.
C     DCOPY from BLAS.
C
C     Local variables.
C
      INTEGER I, I1, J, P1, PN, PNM, PI, PI1, PJ
      DOUBLE PRECISION DZ1, RCOND
      INTEGER IWORK(P), INFO
      DOUBLE PRECISION WORK(3*P)
C
C     Construction of the pre-array WRK.
C
      P1 = P + 1
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
      print *, "First part"
      do i = 1, pn
        print *, wrk(i,:)
      end do
C
C     Second part -  Storing B x Q in the (2,3) block of WRK.
C
      IF (MULTBQ) THEN
         DO 60 I = 1, M
            CALL DCOPY(N, B(1,I), 1, WRK(P1,P+N+I), 1)
   60    CONTINUE
      ELSE
         DO 80 I = 1, M
            CALL DGEMV('N', N, M-I+1, 1.D0, B(1,I), LDB, Q(I,I), 1,
     *                 0.0D0, WRK(P1,P+N+I), 1)
   80    CONTINUE
      END IF
      print *, "Second part"
      do i = 1, pn
        print *, wrk(i,:)
      end do
C
C     Third part -  Storing C x S in the (1,2) block of WRK.
C
      DO 100 I = 1, N
         CALL DGEMV('N', P, N-I+1, 1.D0, C(1,I), LDC, S(I,I), 1, 0.D0,
     *              WRK(1,P+I), 1)
  100 CONTINUE
      print *, "Third part"
      do i = 1, pn
        print *, wrk(i,:)
      end do
C
C     Fourth part  -  Storing A x S in the (2,2) block of WRK.
C
      DO 120 I = 1, N
         CALL DGEMV('N', N, N-I+1, 1.D0, A(1,I), LDA, S(I,I), 1, 0.D0,
     *              WRK(P1,P+I), 1)
  120 CONTINUE
      print *, "Fourth part"
      do i = 1, pn
        print *, wrk(i,:)
      end do
C
C     Triangularization (2 steps).
C
C     Step 1: eliminate the (1,2) block of WRK.
C
      DO 160 I = 1, P
         CALL F06FSF(N, WRK(I,I), WRK(I,P1), LDW, TOL, DZ1)
         I1 = I + 1
         DO 140 J = I1, PN
            CALL F06FUF(N, WRK(I,P1), LDW, DZ1, WRK(J,I), WRK(J,P1),
     *                  LDW)
  140    CONTINUE
         CALL DTRCON('1','L','N',P,WRK,LDW,RCOND,WORK,IWORK,INFO)
         IF (RCOND .LT. TOL) WITHK = .FALSE.
  160 CONTINUE
      print *, "Step 1"
      do i = 1, pn
        print *, wrk(i,:)
      end do
C
C     Step 2: triangularize the remaining (2,2) and (2,3) blocks of WRK.
C
      DO 200 I = 1, N
         PI = P + I
         PI1 = PI + 1
         CALL F06FSF(N+M-I, WRK(PI,PI), WRK(PI,PI1), LDW, TOL, DZ1)
         IF (PI1 .LE. PN) THEN
            DO 180 J = PI1, PN
               CALL F06FUF(N+M-I, WRK(PI,PI1), LDW, DZ1, WRK(J,PI),
     *                     WRK(J,PI1), LDW)
  180       CONTINUE
         END IF
  200 CONTINUE
      print *, "Step 2"
      do i = 1, pn
        print *, wrk(i,:)
      end do
C
C     Output K and S.
C
      IF (WITHK) THEN
         DO 220 J = 1, N
            CALL DCOPY(P, WRK(P+J,1), LDW, K(J,1), LDK)
            CALL DTRSV('L', 'T', 'N', P, WRK, LDW, K(J,1), LDK)
  220    CONTINUE
      END IF
      DO 240 J = 1, N
         PJ = P + J
         CALL DCOPY(N-J+1, WRK(PJ,PJ), 1, S(J,J), 1)
  240 CONTINUE
      RETURN
C
C *** Last line of the SRCF subroutine ********************************
C
      END