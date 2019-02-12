C-------------------------------------------------------------     ************
C                                                                     SASUM
C                                                                  ************
      DOUBLE PRECISION FUNCTION SASUM(N,DX,INCX)
C
C     TAKES THE SUM OF THE ABSOLUTE VALUES.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DTEMP
      INTEGER I,INCX,M,MP1,N,IX
C
      SASUM = 0.0D0
      DTEMP = 0.0D0
      IF (N .LT. 0) STOP
      IF (N .EQ. 0) RETURN
      IF (INCX .EQ. 1) GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      IX = 1
      IF (INCX .LE. 0) IX = (-N+1)*INCX + 1
      DO 10 I = 1,N
        DTEMP = DTEMP + ABS(DX(IX))
        IX = IX + INCX
   10 CONTINUE
      SASUM = DTEMP
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,6)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
        DTEMP = DTEMP + ABS(DX(I))
   30 CONTINUE
      IF (N .LT. 6) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,6
        DTEMP = DTEMP + ABS(DX(I)) + ABS(DX(I+1)) + ABS(DX(I+2))
     *  + ABS(DX(I+3)) + ABS(DX(I+4)) + ABS(DX(I+5))
   50 CONTINUE
   60 SASUM = DTEMP
      RETURN
      END

C-------------------------------------------------------------     ************
C                                                                     SAXPY
C                                                                  ************
      SUBROUTINE SAXPY(N,DA,DX,INCX,DY,INCY)
C
C     CONSTANT TIMES A VECTOR PLUS A VECTOR.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DY(1),DA
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
C
      IF (N .LT. 0) STOP
      IF (N .EQ. 0) RETURN
      IF (DA .EQ. 0.0D0) RETURN
      IF (INCX .EQ. 1 .AND. INCY .EQ. 1) GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF (INCX .LE. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LE. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DY(IY) + DA*DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,4)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
        DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
      IF (N .LT. 4) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        DY(I) = DY(I) + DA*DX(I)
        DY(I+1) = DY(I+1) + DA*DX(I+1)
        DY(I+2) = DY(I+2) + DA*DX(I+2)
        DY(I+3) = DY(I+3) + DA*DX(I+3)
   50 CONTINUE
      RETURN
      END

C-------------------------------------------------------------     ************
C                                                                     SCOPY
C                                                                  ************
      SUBROUTINE  SCOPY(N,DX,INCX,DY,INCY)
C
C     COPIES A VECTOR, X, TO A VECTOR, Y.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DY(1)
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
C
      IF (N .LT. 0) STOP
      IF (N .EQ. 0) RETURN
      IF (INCX .EQ. 1 .AND. INCY .EQ. 1) GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF (INCX .LE. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LE. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,7)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
        DY(I) = DX(I)
   30 CONTINUE
      IF (N .LT. 7) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
        DY(I) = DX(I)
        DY(I+1) = DX(I+1)
        DY(I+2) = DX(I+2)
        DY(I+3) = DX(I+3)
        DY(I+4) = DX(I+4)
        DY(I+5) = DX(I+5)
        DY(I+6) = DX(I+6)
   50 CONTINUE
      RETURN
      END

C-------------------------------------------------------------     ************
C                                                                      SDOT
C                                                                  ************
      DOUBLE PRECISION FUNCTION SDOT(N,DX,INCX,DY,INCY)
C
C     FORMS THE DOT PRODUCT OF TWO VECTORS.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DY(1),DTEMP
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
C
      SDOT = 0.0D0
      DTEMP = 0.0D0
      IF (N .LT. 0) STOP
      IF (N .EQ. 0) RETURN
      IF (INCX .EQ. 1 .AND. INCY .EQ. 1) GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF (INCX .LE. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LE. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DTEMP = DTEMP + DX(IX)*DY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      SDOT = DTEMP
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,5)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
        DTEMP = DTEMP + DX(I)*DY(I)
   30 CONTINUE
      IF (N .LT. 5) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DTEMP = DTEMP + DX(I)*DY(I) + DX(I+1)*DY(I+1) +
     *   DX(I+2)*DY(I+2) + DX(I+3)*DY(I+3) + DX(I+4)*DY(I+4)
   50 CONTINUE
   60 SDOT = DTEMP
      RETURN
      END

C-------------------------------------------------------------     ************
C                                                                     SMACH
C                                                                  ************
      DOUBLE PRECISION FUNCTION SMACH(JOB)
      INTEGER JOB
C
C     SMACH COMPUTES MACHINE PARAMETERS OF FLOATING POINT
C     ARITHMETIC FOR USE IN TESTING ONLY.  NOT REQUIRED BY
C     LINPACK PROPER.
C
C     IF TROUBLE WITH AUTOMATIC COMPUTATION OF THESE QUANTITIES,
C     THEY CAN BE SET BY DIRECT ASSIGNMENT STATEMENTS.
C     ASSUME THE COMPUTER HAS
C
C        B = BASE OF ARITHMETIC
C        T = NUMBER OF BASE  B  DIGITS
C        L = SMALLEST POSSIBLE EXPONENT
C        U = LARGEST POSSIBLE EXPONENT
C
C     THEN
C
C        EPS = B**(1-T)
C        TINY = 100.0*B**(-L+T)
C        HUGE = 0.01*B**(U-T)
C
C     DMACH SAME AS SMACH EXCEPT T, L, U APPLY TO
C     DOUBLE PRECISION.
C
C     CMACH SAME AS SMACH EXCEPT IF COMPLEX DIVISION
C     IS DONE BY
C
C        1/(X+I*Y) = (X-I*Y)/(X**2+Y**2)
C
C     THEN
C
C        TINY = SQRT(TINY)
C        HUGE = SQRT(HUGE)
C
C
C     JOB IS 1, 2 OR 3 FOR EPSILON, TINY AND HUGE, RESPECTIVELY.
C
      DOUBLE PRECISION EPS,TINY,HUGE,S
C
      EPS = 1.0D0
   10 EPS = EPS/2.0D0
      S = 1.0D0 + EPS
      IF (S .GT. 1.0D0) GO TO 10
      EPS = 2.0D0*EPS
C
      S = 1.0D0
   20 TINY = S
      S = S/16.0D0
      IF (S*1.0 .NE. 0.0D0) GO TO 20
      TINY = (TINY/EPS)*100.0
      HUGE = 1.0D0/TINY
C
      IF (JOB .EQ. 1) SMACH = EPS
      IF (JOB .EQ. 2) SMACH = TINY
      IF (JOB .EQ. 3) SMACH = HUGE
      RETURN
      END

C-------------------------------------------------------------     ************
C                                                                     SNRM2
C                                                                  ************
      DOUBLE PRECISION FUNCTION SNRM2(N,DX,INCX)

      INTEGER NEXT
      DOUBLE PRECISION DX(1),CUTLO,CUTHI,HITEST,SUM,XMAX
      DOUBLE PRECISION ZERO,ONE
      DATA ZERO,ONE /0.0D0,1.0D0/
C
C     EUCLIDEAN NORM OF THE N-VECTOR STORED IN DX() WITH STORAGE
C     INCREMENT INCX .
C     IF    N .LE. 0 RETURN WITH RESULT = 0.
C     IF N .GE. 1 THEN INCX MUST BE .GE. 1
C
C           C.L.LAWSON, 1978 JAN 08
C
C     FOUR PHASE METHOD     USING TWO BUILT-IN CONSTANTS THAT ARE
C     HOPEFULLY APPLICABLE TO ALL MACHINES.
C         CUTLO = MAXIMUM OF  DSQRT(U/EPS)  OVER ALL KNOWN MACHINES.
C         CUTHI = MINIMUM OF  DSQRT(V)      OVER ALL KNOWN MACHINES.
C     WHERE
C         EPS = SMALLEST NO. SUCH THAT EPS + 1. .GT. 1.
C         U   = SMALLEST POSITIVE NO.   (UNDERFLOW LIMIT)
C         V   = LARGEST  NO.            (OVERFLOW  LIMIT)
C
C     BRIEF OUTLINE OF ALGORITHM..
C
C     PHASE 1    SCANS ZERO COMPONENTS.
C     MOVE TO PHASE 2 WHEN A COMPONENT IS NONZERO AND .LE. CUTLO
C     MOVE TO PHASE 3 WHEN A COMPONENT IS .GT. CUTLO
C     MOVE TO PHASE 4 WHEN A COMPONENT IS .GE. CUTHI/M
C     WHERE M = N FOR X() REAL AND M = 2*N FOR COMPLEX.
C
C     VALUES FOR CUTLO AND CUTHI..
C     FROM THE ENVIRONMENTAL PARAMETERS LISTED IN THE IMSL CONVERTER
C     DOCUMENT THE LIMITING VALUES ARE AS FOLLOWS..
C     CUTLO, S.P.   U/EPS = 2**(-102) FOR  HONEYWELL.  CLOSE SECONDS ARE
C                   UNIVAC AND DEC AT 2**(-103)
C                   THUS CUTLO = 2**(-51) = 4.44089E-16
C     CUTHI, S.P.   V = 2**127 FOR UNIVAC, HONEYWELL, AND DEC.
C                   THUS CUTHI = 2**(63.5) = 1.30438E19
C     CUTLO, D.P.   U/EPS = 2**(-67) FOR HONEYWELL AND DEC.
C                   THUS CUTLO = 2**(-33.5) = 8.23181D-11
C     CUTHI, D.P.   SAME AS S.P.  CUTHI = 1.30438D19
C     DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /
C     DATA CUTLO, CUTHI / 4.441E-16,  1.304E19 /
      DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /
C
      IF(N .GT. 0) GO TO 10
	 SNRM2  = ZERO
         GO TO 300
C
   10 ASSIGN 30 TO NEXT
      SUM = ZERO
      NN = N * INCX
C                    BEGIN MAIN LOOP
      I = 1
   20    GO TO NEXT,(30, 50, 70, 110)
   30 IF( ABS(DX(I)) .GT. CUTLO) GO TO 85
      ASSIGN 50 TO NEXT
      XMAX = ZERO
C
C                        PHASE 1.  SUM IS ZERO
C
   50 IF( DX(I) .EQ. ZERO) GO TO 200
      IF( ABS(DX(I)) .GT. CUTLO) GO TO 85
C
C                                PREPARE FOR PHASE 2.
      ASSIGN 70 TO NEXT
      GO TO 105
C
C                                PREPARE FOR PHASE 4.
C
  100 I = J
      ASSIGN 110 TO NEXT
      SUM = (SUM / DX(I)) / DX(I)
  105 XMAX = ABS(DX(I))
      GO TO 115
C
C                   PHASE 2.  SUM IS SMALL.
C                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
C
   70 IF( ABS(DX(I)) .GT. CUTLO ) GO TO 75
C
C                     COMMON CODE FOR PHASES 2 AND 4.
C                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
C
  110 IF( ABS(DX(I)) .LE. XMAX ) GO TO 115
         SUM = ONE + SUM * (XMAX / DX(I))**2
         XMAX = ABS(DX(I))
         GO TO 200
C
  115 SUM = SUM + (DX(I)/XMAX)**2
      GO TO 200
C
C
C                  PREPARE FOR PHASE 3.
C
   75 SUM = (SUM * XMAX) * XMAX
C
C
C     FOR REAL OR D.P. SET HITEST = CUTHI/N
C     FOR COMPLEX      SET HITEST = CUTHI/(2*N)
C
   85 HITEST = CUTHI/FLOAT( N )
C
C                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING.
C
      DO 95 J =I,NN,INCX
      IF(ABS(DX(J)) .GE. HITEST) GO TO 100
   95    SUM = SUM + DX(J)**2
      SNRM2 = DSQRT( SUM )
      GO TO 300
C
  200 CONTINUE
      I = I + INCX
      IF ( I .LE. NN ) GO TO 20
C
C              END OF MAIN LOOP.
C
C              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
C
      SNRM2 = XMAX * DSQRT(SUM)
  300 CONTINUE
      RETURN
      END

C-------------------------------------------------------------     ************
C                                                                      SROT
C                                                                  ************
      SUBROUTINE SROT(N,DX,INCX,DY,INCY,C,S)
C
C     APPLIES A PLANE ROTATION.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DY(1),DTEMP,C,S
      INTEGER I,INCX,INCY,IX,IY,N
C
      IF (N .LT. 0) STOP
      IF (N .EQ. 0) RETURN
      IF (INCX .EQ. 1 .AND. INCY .EQ. 1) GO TO 20
C
C       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL
C         TO 1
C
      IX = 1
      IY = 1
      IF (INCX .LE. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LE. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DTEMP = C*DX(IX) + S*DY(IY)
        DY(IY) = C*DY(IY) - S*DX(IX)
        DX(IX) = DTEMP
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C       CODE FOR BOTH INCREMENTS EQUAL TO 1
C
   20 DO 30 I = 1,N
        DTEMP = C*DX(I) + S*DY(I)
        DY(I) = C*DY(I) - S*DX(I)
        DX(I) = DTEMP
   30 CONTINUE
      RETURN
      END

C-------------------------------------------------------------     ************
C                                                                     SROTG
C                                                                  ************
      SUBROUTINE SROTG(DA,DB,C,S)
C
C     CONSTRUCT GIVENS PLANE ROTATION.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DA,DB,C,S,ROE,SCALE,R,Z
C
      ROE = DB
      IF( ABS(DA) .GT. ABS(DB) ) ROE = DA
      SCALE = ABS(DA) + ABS(DB)
      IF( SCALE .NE. 0.0D0 ) GO TO 10
         C = 1.0D0
         S = 0.0D0
         R = 0.0D0
         GO TO 20
   10 R = SCALE*DSQRT((DA/SCALE)**2 + (DB/SCALE)**2)
      R = DSIGN(1.0D0,ROE)*R
      C = DA/R
      S = DB/R
   20 Z = 1.0D0
      IF(ABS(DA) .GT. ABS(DB)) Z = S
      IF(ABS(DB) .GE. ABS(DA) .AND. C .NE. 0.0D0) Z = 1.0D0/C
      DA = R
      DB = Z
      RETURN
      END

C-------------------------------------------------------------     ************
C                                                                     SSCAL
C                                                                  ************
      SUBROUTINE  SSCAL(N,DA,DX,INCX)
C
C     SCALES A VECTOR BY A CONSTANT.
C     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DA,DX(1)
      INTEGER I,INCX,M,MP1,N,IX
C
      IF (N .LT. 0) STOP
      IF (N .EQ. 0) RETURN
      IF (INCX .EQ. 1) GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      IX = 1
      IF (INCX .LE. 0) IX = (-N+1)*INCX + 1
      DO 10 I = 1,N
        DX(IX) = DA*DX(IX)
        IX = IX + INCX
   10 CONTINUE
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,5)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
        DX(I) = DA*DX(I)
   30 CONTINUE
      IF (N .LT. 5) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DX(I) = DA*DX(I)
        DX(I+1) = DA*DX(I+1)
        DX(I+2) = DA*DX(I+2)
        DX(I+3) = DA*DX(I+3)
        DX(I+4) = DA*DX(I+4)
   50 CONTINUE
      RETURN
      END

C-------------------------------------------------------------     ************
C                                                                     SSWAP
C                                                                  ************
      SUBROUTINE  SSWAP (N,DX,INCX,DY,INCY)
C
C     INTERCHANGES TWO VECTORS.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DY(1),DTEMP
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
C
      IF (N .LT. 0) STOP
      IF (N .EQ. 0) RETURN
      IF (INCX .EQ. 1 .AND. INCY .EQ. 1) GO TO 20
C
C       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL
C         TO 1
C
      IX = 1
      IY = 1
      IF (INCX .LE. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LE. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DTEMP = DX(IX)
        DX(IX) = DY(IY)
        DY(IY) = DTEMP
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C       CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C       CLEAN-UP LOOP
C
   20 M = MOD(N,3)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
        DTEMP = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP
   30 CONTINUE
      IF (N .LT. 3) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,3
        DTEMP = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP
        DTEMP = DX(I+1)
        DX(I+1) = DY(I+1)
        DY(I+1) = DTEMP
        DTEMP = DX(I+2)
        DX(I+2) = DY(I+2)
        DY(I+2) = DTEMP
   50 CONTINUE
      RETURN
      END

C-------------------------------------------------------------     ************
C                                                                     ISAMAX
C                                                                  ************
      INTEGER FUNCTION ISAMAX(N,DX,INCX)
C
C     FINDS THE INDEX OF ELEMENT HAVING MAX. ABSOLUTE VALUE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DMAX
      INTEGER I,INCX,IX,N
C
      ISAMAX = 0
      IF (N .LT. 0) STOP
      IF (N .EQ. 0) RETURN
      ISAMAX = 1
      IF (N .EQ. 1) RETURN
      IF (INCX .EQ. 1) GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      IX = 1
      IF (INCX .LE. 0) IX = (-N+1)*INCX + 1
      DMAX = ABS(DX(IX))
      IX = IX + INCX
      DO 10 I = 2,N
         IF (ABS(DX(IX)) .LE. DMAX) GO TO 5
         ISAMAX = I
         DMAX = ABS(DX(IX))
    5    IX = IX + INCX
   10 CONTINUE
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
   20 DMAX = ABS(DX(1))
      DO 30 I = 2,N
         IF (ABS(DX(I)) .LE. DMAX) GO TO 30
	 ISAMAX = I
         DMAX = ABS(DX(I))
   30 CONTINUE
      RETURN
      END

C-------------------------------------------------------------     ************
C                                                                     SPDOT
C                                                                  ************
      DOUBLE PRECISION FUNCTION SPDOT (N,SY,INDEX,SX)
c
      DOUBLE PRECISION SY(1),SX(1),SUM
      INTEGER N,I,INDEX(1)

      SUM = 0.0D0
      DO 10 I=1,N
   10   SUM = SUM + SY(INDEX(I)) * SX(I)
      SPDOT = SUM
      RETURN
      END

C-------------------------------------------------------------     ************
C                                                                     SPAXPY
C                                                                  ************
      SUBROUTINE SPAXPY(N,SA,SX,SY,INDEX)
c
      DOUBLE PRECISION SA,SX(1),SY(1)
      INTEGER N,J,INDEX(1)

      DO 60 J = 1,N
   60   SY(INDEX(J)) = SA * SX(J) + SY(INDEX(J))
      RETURN
      END

C-------------------------------------------------------------     ************
C                                                                     CAXPY
C                                                                  ************
      SUBROUTINE CAXPY(N,CA,CX,INCX,CY,INCY)
C
C     CONSTANT TIMES A VECTOR PLUS A VECTOR.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE COMPLEX CX(1),CY(1),CA
      INTEGER I,INCX,INCY,IX,IY,N
C
      IF (N .LT. 0) STOP
      IF (N .EQ. 0) RETURN
      IF (ABS(DREAL(CA)) + ABS(DIMAG(CA)) .EQ. 0.0D0) RETURN
      IF (INCX .EQ. 1 .AND. INCY .EQ. 1) GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF (INCX .LE. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LE. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        CY(IY) = CY(IY) + CA*CX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
   20 DO 30 I = 1,N
        CY(I) = CY(I) + CA*CX(I)
   30 CONTINUE
      RETURN
      END

C-------------------------------------------------------------     ************
C                                                                     CCOPY
C                                                                  ************
      SUBROUTINE  CCOPY(N,CX,INCX,CY,INCY)
C
C     COPIES A VECTOR, X, TO A VECTOR, Y.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE COMPLEX CX(1),CY(1)
      INTEGER I,INCX,INCY,IX,IY,N
C
      IF (N .LT. 0) STOP
      IF (N .EQ. 0) RETURN
      IF(INCX .EQ. 1 .AND. INCY .EQ. 1) GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF (INCX .LE. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LE. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        CY(IY) = CX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
   20 DO 30 I = 1,N
        CY(I) = CX(I)
   30 CONTINUE
      RETURN
      END

C-------------------------------------------------------------     ************
C                                                                     CDOTC
C                                                                  ************
      DOUBLE COMPLEX FUNCTION CDOTC(N,CX,INCX,CY,INCY)
C
C     FORMS THE DOT PRODUCT OF TWO VECTORS, CONJUGATING THE FIRST
C     VECTOR.
C     JACK DONGARRA, LINPACK,  3/11/78.
C
      DOUBLE COMPLEX CX(1),CY(1),CTEMP
      INTEGER I,INCX,INCY,IX,IY,N
C
      CTEMP = (0.0D0,0.0D0)
      CDOTC = (0.0D0,0.0D0)
      IF (N .LT. 0) STOP
      IF (N .EQ. 0) RETURN
      IF (INCX .EQ. 1 .AND. INCY .EQ. 1) GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF (INCX .LE. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LE. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        CTEMP = CTEMP + CONJG(CX(IX))*CY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      CDOTC = CTEMP
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
   20 DO 30 I = 1,N
        CTEMP = CTEMP + CONJG(CX(I))*CY(I)
   30 CONTINUE
      CDOTC = CTEMP
      RETURN
      END

C-------------------------------------------------------------     ************
C                                                                     CDOTU
C                                                                  ************
      DOUBLE COMPLEX FUNCTION CDOTU(N,CX,INCX,CY,INCY)
C
C     FORMS THE DOT PRODUCT OF TWO VECTORS.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE COMPLEX CX(1),CY(1),CTEMP
      INTEGER I,INCX,INCY,IX,IY,N
C
      CTEMP = (0.0D0,0.0D0)
      CDOTU = (0.0D0,0.0D0)
      IF (N .LT. 0) STOP
      IF (N .EQ. 0) RETURN
      IF (INCX .EQ. 1 .AND. INCY .EQ. 1) GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF (INCX .LE. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LE. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        CTEMP = CTEMP + CX(IX)*CY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      CDOTU = CTEMP
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
   20 DO 30 I = 1,N
        CTEMP = CTEMP + CX(I)*CY(I)
   30 CONTINUE
      CDOTU = CTEMP
      RETURN
      END

C-------------------------------------------------------------     ************
C                                                                     CMACH
C                                                                  ************
      DOUBLE PRECISION FUNCTION CMACH(JOB)
      INTEGER JOB
C
C     SMACH COMPUTES MACHINE PARAMETERS OF FLOATING POINT
C     ARITHMETIC FOR USE IN TESTING ONLY.  NOT REQUIRED BY
C     LINPACK PROPER.
C
C     IF TROUBLE WITH AUTOMATIC COMPUTATION OF THESE QUANTITIES,
C     THEY CAN BE SET BY DIRECT ASSIGNMENT STATEMENTS.
C     ASSUME THE COMPUTER HAS
C
C        B = BASE OF ARITHMETIC
C        T = NUMBER OF BASE  B  DIGITS
C        L = SMALLEST POSSIBLE EXPONENT
C        U = LARGEST POSSIBLE EXPONENT
C
C     THEN
C
C        EPS = B**(1-T)
C        TINY = 100.0*B**(-L+T)
C        HUGE = 0.01*B**(U-T)
C
C     DMACH SAME AS SMACH EXCEPT T, L, U APPLY TO
C     DOUBLE PRECISION.
C
C     CMACH SAME AS SMACH EXCEPT IF COMPLEX DIVISION
C     IS DONE BY
C
C        1/(X+I*Y) = (X-I*Y)/(X**2+Y**2)
C
C     THEN
C
C        TINY = SQRT(TINY)
C        HUGE = SQRT(HUGE)
C
C
C     JOB IS 1, 2 OR 3 FOR EPSILON, TINY AND HUGE, RESPECTIVELY.
C
C
      DOUBLE PRECISION EPS,TINY,HUGE,S
C
      EPS = 1.0D0
   10 EPS = EPS/2.0D0
      S = 1.0 + EPS
      IF (S .GT. 1.0D0) GO TO 10
      EPS = 2.0D0*EPS
      CMACH =EPS
      IF( JOB .EQ. 1) RETURN
C
      S = 1.0D0
   20 TINY = S
      S = S/16.0D0
      IF (S*1.0D0 .NE. 0.0D0) GO TO 20
      TINY = (TINY/EPS)*100.0D0
      S = DREAL((1.0D0,0.0D0)/DCMPLX(TINY,0.0D0))
      IF (S .NE. 1.0D0/TINY) TINY = DSQRT(TINY)
      HUGE = 1.0D0/TINY
      IF (JOB .EQ. 1) CMACH = EPS
      IF (JOB .EQ. 2) CMACH = TINY
      IF (JOB .EQ. 3) CMACH = HUGE
      RETURN
      END

C-------------------------------------------------------------     ************
C                                                                      CROT
C                                                                  ************
      SUBROUTINE CROT(N,CX,INCX,CY,INCY,CC,CS)
C
C     APPLIES A PLANE ROTATION.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE COMPLEX CX(1),CY(1),CTEMP,CC,CS
      INTEGER I,INCX,INCY,IX,IY,N
C
      IF (N .LT. 0) STOP
      IF (N .EQ. 0) RETURN
      IF (INCX .EQ. 1 .AND. INCY .EQ. 1) GO TO 20
C
C       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL
C         TO 1.
C
      IX = 1
      IY = 1
      IF (INCX .LE. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LE. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        CTEMP = CC*CX(IX) + CS*CY(IY)
        CY(IY) = CC*CY(IY) - CONJG(CS)*CX(IX)
        CX(IX) = CTEMP
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C       CODE FOR BOTH INCREMENTS EQUAL TO 1
C
   20 DO 30 I = 1,N
        CTEMP = CC*CX(I) + CS*CY(I)
        CY(I) = CC*CY(I) - CONJG(CS)*CX(I)
        CX(I) = CTEMP
   30 CONTINUE
      RETURN
      END

C-------------------------------------------------------------     ************
C                                                                     CROTG
C                                                                  ************
      SUBROUTINE CROTG(CA,CB,C,S)

      DOUBLE COMPLEX CA,CB,S
      DOUBLE PRECISION C
      DOUBLE PRECISION NORM,SCALE
      DOUBLE COMPLEX ALPHA

      IF (ABS(CA) .NE. 0.0D0) GO TO 10
         C = 0.0D0
         S = (1.0D0,0.0D0)
         CA = CB
         GO TO 20
   10 CONTINUE
         SCALE = ABS(CA) + ABS(CB)
         NORM=SCALE*DSQRT((ABS(CA/SCALE))**2+(ABS(CB/SCALE))**2)
         ALPHA = CA/ABS(CA)
         C = ABS(CA)/NORM
         S = ALPHA * CONJG(CB)/NORM
         CA = ALPHA * NORM
   20 CONTINUE
      RETURN
      END

C-------------------------------------------------------------     ************
C                                                                     CSCAL
C                                                                  ************
      SUBROUTINE CSCAL(N,CA,CX,INCX)
C
C     SCALES A VECTOR BY A CONSTANT.
C     JACK DONGARRA, LINPACK,  3/11/78.
C
      DOUBLE COMPLEX CA,CX(1)
      INTEGER I,INCX,N,IX
C
      IF (N .LT. 0) STOP
      IF (N .EQ. 0) RETURN
      IF (INCX .EQ. 1) GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1.
C
      IX = 1
      IF (INCX .LE. 0) IX = (-N+1)*INCX + 1
      DO 10 I = 1,N
        CX(IX) = CA*CX(IX)
        IX = IX + INCX
   10 CONTINUE
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1.
C
   20 DO 30 I = 1,N
        CX(I) = CA*CX(I)
   30 CONTINUE
      RETURN
      END

C-------------------------------------------------------------     ************
C                                                                     CSSCAL
C                                                                  ************
      SUBROUTINE  CSSCAL(N,SA,CX,INCX)
C
C     SCALES A COMPLEX VECTOR BY A REAL CONSTANT.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE COMPLEX CX(1)
      DOUBLE PRECISION SA
      INTEGER I,INCX,N,IX
C
      IF (N .LT. 0) STOP
      IF (N .EQ. 0) RETURN
      IF (INCX .EQ. 1) GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1.
C
      IX = 1
      IF (INCX .LE. 0) IX = (-N+1)*INCX + 1
      DO 10 I = 1,N
        CX(IX) = DCMPLX(SA*DREAL(CX(IX)),SA*DIMAG(CX(IX)))
        IX = IX + INCX
   10 CONTINUE
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1.
C
   20 DO 30 I = 1,N
        CX(I) = DCMPLX(SA*DREAL(CX(I)),SA*DIMAG(CX(I)))
   30 CONTINUE
      RETURN
      END

C-------------------------------------------------------------     ************
C                                                                      CSUM
C                                                                  ************
      DOUBLE COMPLEX FUNCTION CSUM(N,CX,INCX)
C
C        TAKES THE SUM OF THE VALUES OF A VECTOR.
C        USES UNROLLED LOOPS FOR INCREMENT EQUAL TO ONE.
C        JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE COMPLEX CX(1),CTEMP
      INTEGER I,INCX,M,MP1,N,IX
C
      CSUM = (0.0D0,0.0D0)
      CTEMP = (0.0D0,0.0D0)
      IF (N .LT. 0) STOP
      IF (N .EQ. 0) RETURN
      IF (INCX .EQ. 1) GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      IX = 1
      IF (INCX .LE. 0) IX = (-N+1)*INCX + 1
      DO 10 I = 1,N
        CTEMP = CTEMP + CX(IX)
        IX = IX + INCX
   10 CONTINUE
      CSUM = CTEMP
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,6)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
        CTEMP = CTEMP + CX(I)
   30 CONTINUE
      IF (N .LT. 6) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,6
        CTEMP = CTEMP + CX(I) + CX(I+1) + CX(I+2) + CX(I+3)
     *  + CX(I+4) + CX(I+5)
   50 CONTINUE
   60 CSUM = CTEMP
      RETURN
      END

C-------------------------------------------------------------     ************
C                                                                     CSWAP
C                                                                  ************
      SUBROUTINE CSWAP(N,CX,INCX,CY,INCY)
C
C     INTERCHANGES TWO VECTORS.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE COMPLEX CX(1),CY(1),CTEMP
      INTEGER I,INCX,INCY,IX,IY,N
C
      IF (N .LT. 0) STOP
      IF (N .EQ. 0) RETURN
      IF (INCX .EQ. 1 .AND. INCY .EQ. 1) GO TO 20
C
C          CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C            NOT EQUAL TO 1.
C
      IX = 1
      IY = 1
      IF (INCX .LE. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LE. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
         CTEMP = CX(IX)
         CX(IX) = CY(IY)
         CY(IY) = CTEMP
         IX = IX + INCX
         IY = IY + INCY
   10 CONTINUE
      RETURN
C
C          CODE FOR BOTH INCREMENTS EQUAL TO 1.
C
   20 DO 30 I = 1,N
         CTEMP = CX(I)
         CX(I) = CY(I)
         CY(I) = CTEMP
   30 CONTINUE
      RETURN
      END

C-------------------------------------------------------------     ************
C                                                                     ICAMAX
C                                                                  ************
      INTEGER FUNCTION ICAMAX(N,CX,INCX)
C
C     FINDS THE INDEX OF ELEMENT HAVING MAX. ABSOLUTE VALUE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE COMPLEX CX(1)
      DOUBLE PRECISION SMAX
      INTEGER I,INCX,IX,N
      DOUBLE COMPLEX ZDUM
      DOUBLE PRECISION CABS1
      CABS1(ZDUM) = ABS(DREAL(ZDUM)) + ABS(DIMAG(ZDUM))
C
      ICAMAX = 0
      IF (N .LT. 0) STOP
      IF (N .EQ. 0) RETURN
      ICAMAX = 1
      IF (N .EQ. 1) RETURN
      IF (INCX .EQ. 1) GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      IX = 1
      IF (INCX .LE. 0) IX = (-N+1)*INCX + 1
      SMAX = CABS1(CX(IX))
      IX = IX + INCX
      DO 10 I = 2,N
         IF(CABS1(CX(IX)) .LE. SMAX) GO TO 5
         ICAMAX = I
         SMAX = CABS1(CX(IX))
    5    IX = IX + INCX
   10 CONTINUE
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
   20 SMAX = CABS1(CX(1))
      DO 30 I = 2,N
         IF(CABS1(CX(I)) .LE. SMAX) GO TO 30
         ICAMAX = I
         SMAX = CABS1(CX(I))
   30 CONTINUE
      RETURN
      END

C-------------------------------------------------------------     ************
C                                                                     ISAMIN
C                                                                  ************
      INTEGER FUNCTION ISAMIN(N,DX,INCX)
C
C        FINDS THE INDEX OF ELEMENT HAVING MIN. ABSOLUTE VALUE.
C        JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DMIN
      INTEGER I,INCX,IX,N
C
      ISAMIN = 0
      IF (N .LT. 0) STOP
      IF (N .EQ. 0) RETURN
      ISAMIN = 1
      IF (N .EQ. 1) RETURN
      IF (INCX .EQ. 1) GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      IX = 1
      IF (INCX .LE. 0) IX = (-N+1)*INCX + 1
      DMIN = ABS(DX(IX))
      IX = IX + INCX
      DO 10 I = 2,N
         IF (ABS(DX(IX)) .GE. DMIN) GO TO 5
         ISAMIN = I
         DMIN = ABS(DX(IX))
    5    IX = IX + INCX
   10 CONTINUE
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
   20 DMIN = ABS(DX(1))
      DO 30 I = 2,N
         IF (ABS(DX(I)) .GE. DMIN) GO TO 30
         ISAMIN = I
         DMIN = ABS(DX(I))
   30 CONTINUE
      RETURN
      END

C-------------------------------------------------------------     ************
C                                                                     ISMAX
C                                                                  ************
      INTEGER FUNCTION ISMAX(N,DX,INCX)
C
C        FINDS THE INDEX OF ELEMENT WITH MAX. VALUE.
C        JACK DONGARRA, LINPACK, 3/11/78.
C
      INTEGER I,INCX,IX,N
      DOUBLE PRECISION DX(1),DMAX
C
      ISMAX = 0
      IF (N .LT. 0) STOP
      IF (N .EQ. 0) RETURN
      ISMAX = 1
      IF (N .EQ. 1) RETURN
      IF (INCX .EQ. 1) GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      IX = 1
      IF (INCX .LE. 0) IX = (-N+1)*INCX + 1
      DMAX = DX(IX)
      IX = IX + INCX
      DO 10 I = 2,N
        IF (DX(IX) .LE. DMAX) GO TO 5
        ISMAX = I
        DMAX = DX(IX)
    5   IX = IX + INCX
   10 CONTINUE
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
   20 DMAX = DX(1)
      DO 30 I = 2,N,1
        IF (DX(I) .LE. DMAX) GO TO 30
        ISMAX = I
        DMAX = DX(I)
   30 CONTINUE
      RETURN
      END

C-------------------------------------------------------------     ************
C                                                                     ISMIN
C                                                                  ************
      INTEGER FUNCTION ISMIN(N,DX,INCX)
C
C        FINDS THE INDEX OF ELEMENT WITH MIN. VALUE.
C        JACK DONGARRA, LINPACK, 3/11/78.
C
      INTEGER I,INCX,IX,N
      DOUBLE PRECISION DX(1),DMIN
C
      ISMIN = 0
      IF (N .LT. 0) STOP
      IF (N .EQ. 0) RETURN
      ISMIN = 1
      IF (N .EQ. 1) RETURN
      IF (INCX .EQ. 1) GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      IX = 1
      IF (INCX .LE. 0) IX = (-N+1)*INCX + 1
      DMIN = DX(IX)
      IX = IX + INCX
      DO 10 I = 2,N
        IF (DX(IX) .GE. DMIN) GO TO 5
        ISMIN = I
        DMIN = DX(IX)
    5   IX = IX + INCX
   10 CONTINUE
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
   20 DMIN = DX(1)
      DO 30 I = 2,N,1
        IF (DX(I) .GE. DMIN) GO TO 30
        ISMIN = I
        DMIN = DX(I)
   30 CONTINUE
      RETURN
      END

C-------------------------------------------------------------     ************
C                                                                     SCASUM
C                                                                  ************
      DOUBLE PRECISION FUNCTION SCASUM(N,CX,INCX)
C
C     TAKES THE SUM OF THE ABSOLUTE VALUES OF A COMPLEX VECTOR AND
C     RETURNS A SINGLE PRECISION RESULT.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE COMPLEX CX(1)
      DOUBLE PRECISION STEMP
      INTEGER I,INCX,N,IX
C
      SCASUM = 0.0D0
      STEMP = 0.0D0
      IF (N .LT. 0) STOP
      IF (N .EQ. 0) RETURN
      IF (INCX .EQ. 1) GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      IX = 1
      IF (INCX .LE. 0) IX = (-N+1)*INCX + 1
      DO 10 I = 1,N
        STEMP = STEMP+ABS(DREAL(CX(IX)))+ABS(DIMAG(CX(IX)))
        IX = IX + INCX
   10 CONTINUE
      SCASUM = STEMP
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
   20 DO 30 I = 1,N
        STEMP = STEMP+ABS(DREAL(CX(I)))+ABS(DIMAG(CX(I)))
   30 CONTINUE
      SCASUM = STEMP
      RETURN
      END

C-------------------------------------------------------------     ************
C                                                                     SCNRM2
C                                                                  ************
      DOUBLE PRECISION FUNCTION SCNRM2(N,CX,INCX)
C
      LOGICAL IMAG, SCALE
      INTEGER NEXT
      DOUBLE PRECISION CUTLO,CUTHI,HITEST,SUM,XMAX,ABSX,ZERO,ONE
      DOUBLE COMPLEX CX(1)
      DATA ZERO, ONE /0.0D0, 1.0D0/
C
C     UNITARY NORM OF THE COMPLEX N-VECTOR STORED IN CX() WITH STORAGE
C     INCREMENT INCX .
C     IF    N .LE. 0 RETURN WITH RESULT = 0.
C     IF N .GE. 1 THEN INCX MUST BE .GE. 1
C
C           C.L.LAWSON , 1978 JAN 08
C
C     FOUR PHASE METHOD     USING TWO BUILT-IN CONSTANTS THAT ARE
C     HOPEFULLY APPLICABLE TO ALL MACHINES.
C         CUTLO = MAXIMUM OF  SQRT(U/EPS)  OVER ALL KNOWN MACHINES.
C         CUTHI = MINIMUM OF  SQRT(V)      OVER ALL KNOWN MACHINES.
C     WHERE
C         EPS = SMALLEST NO. SUCH THAT EPS + 1. .GT. 1.
C         U   = SMALLEST POSITIVE NO.   (UNDERFLOW LIMIT)
C         V   = LARGEST  NO.            (OVERFLOW  LIMIT)
C
C     BRIEF OUTLINE OF ALGORITHM..
C
C     PHASE 1    SCANS ZERO COMPONENTS.
C     MOVE TO PHASE 2 WHEN A COMPONENT IS NONZERO AND .LE. CUTLO
C     MOVE TO PHASE 3 WHEN A COMPONENT IS .GT. CUTLO
C     MOVE TO PHASE 4 WHEN A COMPONENT IS .GE. CUTHI/M
C     WHERE M = N FOR X() REAL AND M = 2*N FOR COMPLEX.
C
C     VALUES FOR CUTLO AND CUTHI..
C     FROM THE ENVIRONMENTAL PARAMETERS LISTED IN THE IMSL CONVERTER
C     DOCUMENT THE LIMITING VALUES ARE AS FOLLOWS..
C     CUTLO, S.P.   U/EPS = 2**(-102) FOR  HONEYWELL.  CLOSE SECONDS ARE
C                   UNIVAC AND DEC AT 2**(-103)
C                   THUS CUTLO = 2**(-51) = 4.44089E-16
C     CUTHI, S.P.   V = 2**127 FOR UNIVAC, HONEYWELL, AND DEC.
C                   THUS CUTHI = 2**(63.5) = 1.30438E19
C     CUTLO, D.P.   U/EPS = 2**(-67) FOR HONEYWELL AND DEC.
C                   THUS CUTLO = 2**(-33.5) = 8.23181D-11
C     CUTHI, D.P.   SAME AS S.P.  CUTHI = 1.30438D19
C     DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /
C     DATA CUTLO, CUTHI / 4.441E-16,  1.304E19 /
      DATA CUTLO, CUTHI / 4.441D-16,  1.304D19 /
C
      IF(N .GT. 0) GO TO 10
         SCNRM2  = ZERO
         GO TO 300
C
   10 ASSIGN 30 TO NEXT
      SUM = ZERO
      NN = N * INCX
C                    BEGIN MAIN LOOP
      DO 210 I=1,NN,INCX
         ABSX = ABS(DREAL(CX(I)))
         IMAG = .FALSE.
         GO TO NEXT,(30, 50, 70, 90, 110)
   30 IF(ABSX .GT. CUTLO) GO TO 85
      ASSIGN 50 TO NEXT
      SCALE = .FALSE.
C
C                        PHASE 1.  SUM IS ZERO
C
   50 IF (ABSX .EQ. ZERO) GO TO 200
      IF (ABSX .GT. CUTLO) GO TO 85
C
C                                PREPARE FOR PHASE 2.
      ASSIGN 70 TO NEXT
      GO TO 105
C
C                                PREPARE FOR PHASE 4.
C
  100 ASSIGN 110 TO NEXT
      SUM = (SUM / ABSX) / ABSX
  105 SCALE = .TRUE.
      XMAX = ABSX
      GO TO 115
C
C                   PHASE 2.  SUM IS SMALL.
C                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
C
   70 IF( ABSX .GT. CUTLO ) GO TO 75
C
C                     COMMON CODE FOR PHASES 2 AND 4.
C                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
C
  110 IF (ABSX .LE. XMAX) GO TO 115
         SUM = ONE + SUM * (XMAX / ABSX)**2
         XMAX = ABSX
         GO TO 200
C
  115 SUM = SUM + (ABSX/XMAX)**2
      GO TO 200
C
C
C                  PREPARE FOR PHASE 3.
C
   75 SUM = (SUM * XMAX) * XMAX
C
   85 ASSIGN 90 TO NEXT
      SCALE = .FALSE.
C
C     FOR REAL OR D.P. SET HITEST = CUTHI/N
C     FOR COMPLEX      SET HITEST = CUTHI/(2*N)
C
      HITEST = CUTHI/FLOAT(N)
C
C                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING.
C
   90 IF(ABSX .GE. HITEST) GO TO 100
         SUM = SUM + ABSX**2
  200 CONTINUE
C                  CONTROL SELECTION OF REAL AND IMAGINARY PARTS.
C
      IF(IMAG) GO TO 210
         ABSX = ABS(DIMAG(CX(I)))
         IMAG = .TRUE.
      GO TO NEXT,(50, 70, 90, 110)
C
  210 CONTINUE
C
C              END OF MAIN LOOP.
C              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
C
      SCNRM2 = DSQRT(SUM)
      IF (SCALE) SCNRM2 = SCNRM2 * XMAX
  300 CONTINUE
      RETURN
      END

C-------------------------------------------------------------     ************
C                                                                      SSUM
C                                                                  ************
      DOUBLE PRECISION FUNCTION SSUM(N,DX,INCX)
C
C        TAKES THE SUM OF THE VALUES OF A VECTOR.
C        USES UNROLLED LOOPS FOR INCREMENT EQUAL TO ONE.
C        JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DTEMP
      INTEGER I,INCX,M,MP1,N,IX
C
      SSUM = 0.0D0
      DTEMP = 0.0D0
      IF (N .LT. 0) STOP
      IF (N .EQ. 0) RETURN
      IF (INCX .EQ. 1) GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      IX = 1
      IF (INCX .LE. 0) IX = (-N+1)*INCX + 1
      DO 10 I = 1,N
        DTEMP = DTEMP + DX(IX)
        IX = IX + INCX
   10 CONTINUE
      SSUM = DTEMP
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,6)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
        DTEMP = DTEMP + DX(I)
   30 CONTINUE
      IF (N .LT. 6) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,6
        DTEMP = DTEMP + DX(I) + DX(I+1) + DX(I+2) + DX(I+3)
     *  + DX(I+4) + DX(I+5)
   50 CONTINUE
   60 SSUM = DTEMP
      RETURN
      END

C-------------------------------------------------------------     ************
C                                                                     SROTMG
C                                                                  ************
      SUBROUTINE SROTMG(D1,D2,B1,B2,PARAM)
C
C     This subroutine computes a modified Givens plane rotation matrix.
C     Martin J. McBride, 7/10/85.
C     General Electric CRD, Information System Operation.
C
      DOUBLE PRECISION D1,D2,B1,B2,PARAM(5)
      DOUBLE PRECISION H(2,2),U,M,TEMP
      INTEGER I
C
      M = 4096.0D0
C  Case 4:  If D2 or B2 are equal to 0, then PARAM(1) = 2 is the only
C           change.
      IF (D2*B2 .EQ. 0.0D0) THEN
         PARAM(1) = -2.0D0
         RETURN
      ENDIF

C  If D1 is less than 0, then PARAM(1) = -1 and PARAM(2) to PARAM(5)
C           are set to 0.
      IF (D1 .LT. 0.0D0) THEN
         D1 = 0.0D0
         D2 = 0.0D0
         B1 = 0.0D0
         PARAM(1) = -1.0D0
         DO 5 I = 2,5
            PARAM(I) = 0.0D0
    5    CONTINUE
         RETURN
      ENDIF

      IF (ABS(D1*B1*B1) .GT. ABS(D2*B2*B2)) THEN
C  Case 1:  D1 and D2 are updated by a factor of 1/U, where U is the
C           determinant of matrix H.
         H(1,1) = 1.0D0
         H(1,2) = (D2*B2)/(D1*B1)
         H(2,1) = -B2/B1
         H(2,2) = 1.0D0
         PARAM(1) = 0.0D0
         PARAM(3) = H(2,1)
         PARAM(4) = H(1,2)
         U = 1.0D0 + (D2*B2*B2)/(D1*B1*B1)
         D1 = D1/U
         D2 = D2/U
         B1 = B1*U
         IF (D1 .EQ. 0.0D0 .OR. D2 .EQ. 0.0D0) RETURN
         IF (ABS(D1) .LT. 1.0/(M*M) .OR. ABS(D1) .GT. M*M) THEN
            CALL SROTMGW(D1,1,B1,H,PARAM,M)
         ENDIF
         IF (ABS(D2) .LT. 1.0/(M*M) .OR. ABS(D2) .GT. M*M) THEN
            CALL SROTMGW(D2,2,B1,H,PARAM,M)
         ENDIF
      ELSE
C    If D2 is less than 0, then PARAM(1) = -1 and PARAM(2) to PARAM(5)
C             are set to 0.
         IF (D2 .LT. 0.0D0) THEN
            D1 = 0.0D0
            D2 = 0.0D0
            B1 = 0.0D0
            PARAM(1) = -1.0D0
            DO 15 I = 2,5
               PARAM(I) = 0.0D0
   15       CONTINUE
            RETURN
         ENDIF
C  Case 2:  D1 and D2 are updated by a factor of 1/U, where U is the
C           determinant of matrix H, and then are interchanged.
         H(1,1) = (D1*B1)/(D2*B2)
         H(1,2) = 1.0D0
         H(2,1) = -1.0D0
         H(2,2) = B1/B2
         PARAM(1) = 1.0D0
         PARAM(2) = H(1,1)
         PARAM(5) = H(2,2)
         U = 1.0D0 + (D1*B1*B1)/(D2*B2*B2)
         TEMP = D1/U
         D1 = D2/U
         D2 = TEMP
         B1 = B2*U
         IF (D1 .EQ. 0.0D0 .OR. D2 .EQ. 0.0D0) RETURN
         IF (ABS(D1) .LT. 1.0/(M*M) .OR. ABS(D1) .GT. M*M) THEN
            CALL SROTMGW(D1,1,B1,H,PARAM,M)
         ENDIF
         IF (ABS(D2) .LT. 1.0/(M*M) .OR. ABS(D2) .GT. M*M) THEN
            CALL SROTMGW(D2,2,B1,H,PARAM,M)
         ENDIF
      ENDIF
      RETURN
      END

      SUBROUTINE SROTMGW(DI,I,B1,H,PARAM,M)

      DOUBLE PRECISION DI,B1,PARAM(5),H(2,2),M
      INTEGER I

C  Case 3:  D1 or D2 is updated until it is within a given window,
C           m**-2 <= abs(D(i)) <= m**2 where m=4096.
   10    IF (ABS(DI) .LT. 1.0/(M*M)) THEN
            DI = DI*M*M
            IF (I .EQ. 1) B1 = B1/M
            H(I,1) = H(I,1)/M
            H(I,2) = H(I,2)/M
         ELSE IF (ABS(DI) .GT. M*M) THEN
                 DI = DI*1.0/(M*M)
                 IF (I .EQ. 1) B1 = B1*M
                 H(I,1) = H(I,1)*M
                 H(I,2) = H(I,2)*M
         ENDIF
      IF (ABS(DI) .LT. 1.0/(M*M) .OR. ABS(DI) .GT. M*M) GO TO 10

      PARAM(1) = -1.0D0
      PARAM(2) = H(1,1)
      PARAM(3) = H(2,1)
      PARAM(4) = H(1,2)
      PARAM(5) = H(2,2)

      RETURN
      END

C-------------------------------------------------------------     ************
C                                                                     SROTM
C                                                                  ************
      SUBROUTINE SROTM(N,SX,INCX,SY,INCY,PARAM)
C
C     This subroutine applies the modified Givens rotation matrix.
C     Martin J. McBride.  7/11/85.
C     General Electric CRD, Information System Operation.
C
      INTEGER N,INCX,INCY,IX,IY,I
      DOUBLE PRECISION SX(1),SY(1),PARAM(1),H(2,2)

      IF (N .LT. 0) STOP
      IF (N .EQ. 0) RETURN
      IF (PARAM(1) .EQ. -2.0D0) RETURN

C  Conditions for setting up matrix H for multiplication.
      IF (PARAM(1) .EQ. 1.0D0) THEN
         H(1,1) = PARAM(2)
         H(2,1) = -1.0D0
         H(1,2) = 1.0D0
         H(2,2) = PARAM(5)
      ELSE IF (PARAM(1) .EQ. 0.0D0) THEN
         H(1,1) = 1.0D0
         H(2,1) = PARAM(3)
         H(1,2) = PARAM(4)
         H(2,2) = 1.0D0
      ELSE IF (PARAM(1) .EQ. -1.0D0) THEN
         H(1,1) = PARAM(2)
         H(2,1) = PARAM(3)
         H(1,2) = PARAM(4)
         H(2,2) = PARAM(5)
      ELSE
         PRINT*
         PRINT*,'     SROTM called with incorrect parameter key'
         PRINT*
         RETURN
      ENDIF

      IF (INCX .EQ. 1 .AND. INCY .EQ. 1) GO TO 20
C
C        Code for unequal increments of vectors X and Y or
C          equal increments not equal to 1.
C
      IX = 1
      IY = 1
      IF (INCX .LE. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LE. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
         X = SX(IX)
         Y = SY(IY)
         SX(IX) = X*H(1,1) + Y*H(1,2)
         SY(IY) = X*H(2,1) + Y*H(2,2)
         IX = IX + INCX
         IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        Code for equal increments of vectors X and Y.
C
   20 DO 30 I = 1,N
         X = SX(I)
         Y = SY(I)
         SX(I) = X*H(1,1) + Y*H(1,2)
         SY(I) = X*H(2,1) + Y*H(2,2)
   30 CONTINUE
      RETURN
      END
