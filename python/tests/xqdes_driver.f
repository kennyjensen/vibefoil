      PROGRAM XQDESDR
      INCLUDE 'XFOIL.INC'
      INTEGER NIN, NSPIN, IQ1IN, IQ2IN, NITERQ
      INTEGER I
      INTEGER LCPXXI, LIMAGEI, SHARPI
      REAL ALFAI, MINFI, QINFI, XCMREFI, YCMREFI
      REAL PSIOI, QDOF0I, QDOF1I, QDOF2I, QDOF3I
      REAL ANTEI, ASTEI, DSTEI, XTEI, YTEI, YIMAGEI
      REAL PI_VAL

      READ(*,*,END=900) NIN, NSPIN, IQ1IN, IQ2IN, NITERQ
      READ(*,*) ALFAI, MINFI, QINFI, XCMREFI, YCMREFI, PSIOI,
     &          QDOF0I, QDOF1I, QDOF2I, QDOF3I
      READ(*,*) LCPXXI, LIMAGEI, SHARPI
      READ(*,*) ANTEI, ASTEI, DSTEI, XTEI, YTEI, YIMAGEI

      N = NIN
      NSP = NSPIN
      IQ1 = IQ1IN
      IQ2 = IQ2IN

      ALFA = ALFAI
      MINF = MINFI
      QINF = QINFI
      XCMREF = XCMREFI
      YCMREF = YCMREFI
      PSIO = PSIOI
      QDOF0 = QDOF0I
      QDOF1 = QDOF1I
      QDOF2 = QDOF2I
      QDOF3 = QDOF3I

      LCPXX = (LCPXXI.NE.0)
      LIMAGE = (LIMAGEI.NE.0)
      SHARP = (SHARPI.NE.0)

      ANTE = ANTEI
      ASTE = ASTEI
      DSTE = DSTEI
      XTE = XTEI
      YTE = YTEI
      YIMAGE = YIMAGEI

      PI_VAL = 4.0*ATAN(1.0)
      PI = PI_VAL
      HOPI = 0.5/PI_VAL
      QOPI = 0.25/PI_VAL
      DTOR = PI_VAL/180.0

      Z_QINF = 0.0
      Z_ALFA = 0.0
      Z_QDOF0 = 0.0
      Z_QDOF1 = 0.0
      Z_QDOF2 = 0.0
      Z_QDOF3 = 0.0

      DO I=1, N
        READ(*,*) X(I)
      END DO
      DO I=1, N
        READ(*,*) Y(I)
      END DO
      DO I=1, N
        READ(*,*) S(I)
      END DO
      DO I=1, N
        READ(*,*) XP(I)
      END DO
      DO I=1, N
        READ(*,*) YP(I)
      END DO
      DO I=1, N
        READ(*,*) NX(I)
      END DO
      DO I=1, N
        READ(*,*) NY(I)
      END DO
      DO I=1, N
        READ(*,*) APANEL(I)
      END DO

      DO I=1, NSP
        READ(*,*) SSPEC(I)
      END DO
      DO I=1, NSP
        READ(*,*) QSPEC(I,1)
      END DO
      DO I=1, N
        READ(*,*) GAM(I)
      END DO
      DO I=1, N
        READ(*,*) SIG(I)
      END DO

C---- derive leading edge and chord length for consistent TE logic
      CALL LEFIND(SLE,X,XP,Y,YP,S,N)
      XLE = SEVAL(SLE,X,XP,S,N)
      YLE = SEVAL(SLE,Y,YP,S,N)
      CHORD = SQRT((XTE - XLE)**2 + (YTE - YLE)**2)

C---- Clear derived arrays used by MIXED/PSILIN
      DO I=1, N
        GAM_A(I) = 0.0
        GAMU(I,1) = 0.0
        GAMU(I,2) = 0.0
        QF0(I) = 0.0
        QF1(I) = 0.0
        QF2(I) = 0.0
        QF3(I) = 0.0
        DZDG(I) = 0.0
        DZDN(I) = 0.0
        DQDG(I) = 0.0
        DZDM(I) = 0.0
        DQDM(I) = 0.0
      END DO

      CALL MIXED(1, NITERQ)

      WRITE(*,*) N
      DO I=1, N
        WRITE(*,'(3(1X,ES24.16))') X(I), Y(I), GAM(I)
      END DO
      WRITE(*,'(10(1X,ES24.16))') PSIO, QDOF0, QDOF1, QDOF2, QDOF3,
     &    CL, CM, CDP, CL_ALF, CL_MSQ

  900 CONTINUE
      END
