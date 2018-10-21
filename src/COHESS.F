      SUBROUTINE COHESS(A, LDA, N, B, LDB, M, U, LDU, WITHU, UPPER)
C
C     PURPOSE:
C
C     COHESS computes a unitary state space transformation U reducing
C     the pair (A,B) to upper or lower controller Hessenberg form.
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
C     .. Scalar Arguments ..
C
      INTEGER LDA, N, LDB, M, LDU
      LOGICAL WITHU, UPPER
C
C     .. Array Arguments ..
C
      DOUBLE PRECISION A(LDA,*), B(LDB,*), U(LDU,*)
C
C     EXTERNAL ROUTINES:
C
C     F06FSF, F06FUF from NAG-BLAS,
C
C     Local variables.
C
      INTEGER M1, N1, NJ, JJ, II, J, PAR1, PAR2, PAR3, PAR4,
     *        PAR5, PAR6
      DOUBLE PRECISION TL, DZ
C
      TL = 0.0D0
      M1 = M + 1
      N1 = N - 1
C
C     Perform transformations involving both B and A.
C
      DO 120 J = 1, MIN(M,N1)
         NJ = N - J
         IF (UPPER) THEN
            PAR1 = J
            PAR2 = J
            PAR3 = J + 1
            PAR4 = M
            PAR5 = N
         ELSE
            PAR1 = M - J + 1
            PAR2 = NJ + 1
            PAR3 = 1
            PAR4 = M - J
            PAR5 = NJ
         END IF
         CALL F06FSF(NJ, B(PAR2,PAR1), B(PAR3,PAR1), 1, TL, DZ)
C
C        Update A
C
         DO 20 JJ = 1, N
            CALL F06FUF(NJ, B(PAR3,PAR1), 1, DZ, A(PAR2,JJ),
     *                  A(PAR3,JJ), 1)
   20    CONTINUE
         DO 40 II = 1, N
            CALL F06FUF(NJ, B(PAR3,PAR1), 1, DZ, A(II,PAR2),
     *                  A(II,PAR3), LDA)
   40    CONTINUE
         IF (WITHU) THEN
C
C           Update U
C
            DO 60 JJ = 1, N
               CALL F06FUF(NJ, B(PAR3,PAR1), 1, DZ, U(PAR2,JJ),
     *                  U(PAR3,JJ), 1)
   60       CONTINUE
         END IF
         IF (J .NE. M) THEN
C
C           Update B
C
            DO 80 JJ = PAR3, PAR4
               CALL F06FUF(NJ, B(PAR3,PAR1), 1, DZ, B(PAR2,JJ),
     *                     B(PAR3,JJ), 1)
   80       CONTINUE
         END IF
         DO 100 II = PAR3, PAR5
            B(II,PAR1) = 0.0D0
  100    CONTINUE
  120 CONTINUE
      IF (M1 .LE. N1) THEN
         DO 240 J = M1, N1
C
C           Perform next transformations only involving A.
C
            NJ = N - J
            IF (UPPER) THEN
               PAR1 = J - M
               PAR2 = J
               PAR3 = J + 1
               PAR4 = N
               PAR5 = J - M + 1
               PAR6 = N
            ELSE
               PAR1 = N + M1 - J
               PAR2 = NJ + 1
               PAR3 = 1
               PAR4 = NJ
               PAR5 = 1
               PAR6 = N + M - J
            END IF
            CALL F06FSF(NJ, A(PAR2,PAR1), A(PAR3,PAR1), 1, TL, DZ)
C
C           Update A
C
            DO 160 JJ = PAR5, PAR6
               CALL F06FUF(NJ, A(PAR3,PAR1), 1, DZ, A(PAR2,JJ),
     *                     A(PAR3,JJ), 1)
  160       CONTINUE
            DO 180 II = 1, N
               CALL F06FUF(NJ, A(PAR3,PAR1), 1, DZ, A(II,PAR2),
     *                     A(II,PAR3), LDA)
  180       CONTINUE
            IF (WITHU) THEN
C
C              Update U
C
               DO 200 JJ = 1, N
                  CALL F06FUF(NJ, A(PAR3,PAR1), 1, DZ, U(PAR2,JJ),
     *                        U(PAR3,JJ), 1)
  200          CONTINUE
            END IF
            DO 220 II = PAR3, PAR4
               A(II,PAR1) = 0.0D0
  220       CONTINUE
  240    CONTINUE
      END IF
      RETURN
C
C *** Last line of the COHESS subroutine ******************************
C
      END