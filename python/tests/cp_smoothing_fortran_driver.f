C***********************************************************************
C  Fortran driver to compute inviscid Cp from panel data via GGCALC/QISET.
C***********************************************************************
      PROGRAM CPSMOOTH
      INCLUDE 'XFOIL.INC'
      INTEGER I, NIN
      INTEGER SHARPF
      REAL ALFA_IN, MINF_IN, QINF_IN
      REAL ANTE_IN, ASTE_IN, DSTE_IN
      REAL QOPI_IN, HOPI_IN, PI_IN
      REAL CPLOC(10000)

      READ(*,*,END=900) NIN
      READ(*,*) ALFA_IN, MINF_IN, QINF_IN
      READ(*,*) SHARPF
      READ(*,*) ANTE_IN, ASTE_IN, DSTE_IN
      READ(*,*) QOPI_IN, HOPI_IN, PI_IN

      N = NIN
      ALFA = ALFA_IN
      MINF = MINF_IN
      QINF = QINF_IN
      SHARP = (SHARPF .NE. 0)
      ANTE = ANTE_IN
      ASTE = ASTE_IN
      DSTE = DSTE_IN
      QOPI = QOPI_IN
      HOPI = HOPI_IN
      PI = PI_IN

      DO I=1, N
        READ(*,*) X(I), Y(I), S(I), XP(I), YP(I), NX(I), NY(I)
      END DO
      DO I=1, N
        READ(*,*) SIG(I)
      END DO
      DO I=1, N
        READ(*,*) QF0(I), QF1(I), QF2(I), QF3(I)
      END DO

      XTE = 0.5*(X(1)+X(N))
      YTE = 0.5*(Y(1)+Y(N))

      CALL GGCALC
      CALL QISET
      CALL CPCALC(N,QINV,QINF,MINF,CPLOC)

      WRITE(*,'(I0)') N
      DO I=1, N
        WRITE(*,'(3ES24.15)') X(I), Y(I), CPLOC(I)
      END DO

  900 CONTINUE
      END
