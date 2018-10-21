C This is the driver for SUBROUTINE COHESS.
C The routine is in TOMS.FOR and uses the file COHESS.IN as
C input file from which N, UPPER, TOL, A and B are read.
C Output is sent to file COHESS.OUT .
C
C     .. Local Scalars ..
      DOUBLE PRECISION SCALE,SUMSQ,TOL
      INTEGER I,J,LDA,LDB,LDU,M,N
      LOGICAL UPPER,WITHU
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION A(10,10),A0(10,10),B(10,5),B0(10,5),U(10,10),
     *                 WORK(10)
C     ..
C     .. External Subroutines ..
      EXTERNAL COHESS,DCOPY,DGEMV,F06FBF,F06FJF
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DSQRT
C     ..
      WITHU = .TRUE.
      LDA = 10
      LDB = 10
      LDU = 10
C
C READ N, UPPER, TOL, A and B from file COHESS.IN
C The matrices A and B are to be read column by column.
C The input dimension M is specified by the loop index.
C
      OPEN(UNIT=1,FILE='COHESS.IN',ACTION='READ')
      OPEN(UNIT=2,FILE="COHESS.OUT",ACTION='WRITE')

      DO 700 M = 2,3
          READ (1,*) N,UPPER,TOL
 
   17     FORMAT (1H ,'*** N =',I3,' M =',I3,' UPPER = ',L3,' TOL =',
     *           D8.1)
 
          IF (M.EQ.2) THEN
              WRITE (2,FMT=940)
 
          ELSE
              WRITE (2,FMT=941)
          END IF
C
C The matrices A and B are read, and also stored in A0 and B0 for
C later use.
C
          READ (1,FMT=52) ((A(I,J),I=1, N),J=1, N)
          DO 51 J = 1,N
              CALL DCOPY(N,A(1,J),1,A0(1,J),1)
   51     CONTINUE
 
   52     FORMAT (D25.15)
 
          READ (1,FMT=52) ((B(I,J),I=1, N),J=1, M)
          DO 61 J = 1,M
              CALL DCOPY(N,B(1,J),1,B0(1,J),1)
   61     CONTINUE
C
C The matrix U is initialized as the identity matrix.
C
          DO 100 J = 1,N
              DO 90 I = 1,N
                  U(I,J) = 0.D0
   90         CONTINUE
              U(J,J) = 1.D0
  100     CONTINUE
          WRITE (2,FMT=950)
          WRITE (2,FMT=17) N,M,UPPER,TOL
          CALL PRMT(A0,LDA,N,N,'A0 matr. ',2,4)
          CALL PRMT(B0,LDB,N,M,'B0 matr. ',2,4)
C
C Calling the controller Hessenberg form routine.
C
          CALL COHESS(A,LDA,N,B,LDB,M,U,LDU,WITHU,UPPER)
C
C The transformed matrices A and B, and the transformation matrix U
C are output.
C
          WRITE (2,FMT=951)
          CALL PRMT(A,LDA,N,N,'A matrix',2,4)
          CALL PRMT(B,LDB,N,M,'B matrix',2,4)
          CALL PRMT(U,LDU,N,N,'U matrix',2,4)
C
C Checking the Frobenius norm of (UxA0-AxU),(UxB0-B) and (UxU'-I)
C without altering A,B,U,A0,B0.
C
          SCALE = 1.D0
          SUMSQ = 0.D0
          DO 400 I = 1,N
              CALL DGEMV('N',N,N,1.D0,A,LDA,U(1,I),1,0.D0,WORK,1)
              CALL DGEMV('N',N,N,1.D0,U,LDU,A0(1,I),1,-1.D0,WORK,1)
              CALL F06FJF(N,WORK,1,SCALE,SUMSQ)
  400     CONTINUE
          SUMSQ = SCALE*DSQRT(SUMSQ)
          WRITE (2,FMT=980) SUMSQ
          SCALE = 1.D0
          SUMSQ = 0.D0
          DO 500 I = 1,M
              CALL DCOPY(N,B(1,I),1,WORK,1)
              CALL DGEMV('N',N,N,1.D0,U,LDU,B0(1,I),1,-1.D0,WORK,1)
              CALL F06FJF(N,WORK,1,SCALE,SUMSQ)
  500     CONTINUE
          SUMSQ = SCALE*DSQRT(SUMSQ)
          WRITE (2,FMT=981) SUMSQ
          SCALE = 1.D0
          SUMSQ = 0.D0
          DO 600 I = 1,N
              CALL F06FBF(N,0.D0,WORK,1)
              WORK(I) = 1.D0
              CALL DGEMV('T',N,N,1.D0,U,LDU,U(1,I),1,-1.D0,WORK,1)
              CALL F06FJF(N,WORK,1,SCALE,SUMSQ)
  600     CONTINUE
          SUMSQ = SCALE*DSQRT(SUMSQ)
          WRITE (2,FMT=982) SUMSQ
          IF (M.EQ.2) THEN
              WRITE (2,FMT=945)
 
          ELSE
              WRITE (2,FMT=946)
          END IF
 
  700 CONTINUE

      CLOSE(1)
      CLOSE(2)
 
  940 FORMAT (' *** First example : upper Hessenberg form')
  941 FORMAT ('1*** Second example : lower Hessenberg form')
  945 FORMAT (
     * '1*** This first example is a well-conditioned one.'/
     * ' *** The above norm tests should be close to the machine'/
     * ' *** precision and A, B, U should be close to those obtained'/
     * ' *** on any other machine.')
  946 FORMAT (
     * '1*** This second example is an ill-conditioned one.'/
     * ' *** Therefore A, B and U can differ a lot from those obtained'/
     * ' *** on another machine, but the above norm tests should yet'/
     * ' *** be close to the machine precision (because of backward'/
     * ' *** stability).')
  950 FORMAT (// ' *** Input Data *************' //)
  951 FORMAT (// ' *** Output Data ************' //)
  980 FORMAT (
     * ' *** Testing norm U*A0-A*U =',D12.6,/
     * ' *** which should be of the order of the machine precision'/
     * ' *** norm(A)')
  981 FORMAT (
     * ' *** Testing norm U*B0-B =',D12.6,/
     * ' *** which should be of the order of the machine precision'/
     * ' *** norm(B)')
  982 FORMAT (
     * ' *** Testing norm U*UP-I =',D12.6,/
     * ' *** which should be of the order of the machine precision')
C
C *** Last line of the program COHESS *********************************
 
      END