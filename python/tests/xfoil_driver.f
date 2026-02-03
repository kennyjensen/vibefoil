C***********************************************************************
C  Xfoil driver for TECALC/CPCLAC/CLCALC.
C***********************************************************************
      PROGRAM XFOILDR
      INCLUDE 'XFOIL.INC'
      INTEGER NIN, I
      REAL ALFAI, MINFI, QINFI, XREFI, YREFI
      REAL QLOC(10000), CPLOC(10000)
      REAL SHV

      READ(*,*,END=900) NIN
      READ(*,*) ALFAI, MINFI, QINFI, XREFI, YREFI

      N = NIN
      ALFA = ALFAI
      MINF = MINFI
      QINF = QINFI
      XCMREF = XREFI
      YCMREF = YREFI

      DO I=1, N
        READ(*,*) X(I)
      END DO
      DO I=1, N
        READ(*,*) Y(I)
      END DO
      DO I=1, N
        READ(*,*) GAM(I)
      END DO
      DO I=1, N
        READ(*,*) GAM_A(I)
      END DO
      DO I=1, N
        READ(*,*) QLOC(I)
      END DO

      CALL SCALC(X,Y,S,N)
      CALL SEGSPL(X,XP,S,N)
      CALL SEGSPL(Y,YP,S,N)

      CHORD = 1.0
      CALL TECALC

      CALL CPCALC(N,QLOC,QINF,MINF,CPLOC)
      CALL CLCALC(N,X,Y,GAM,GAM_A,ALFA,MINF,QINF, XREFI, YREFI,
     &            CL,CM,CDP, CL_ALF, CL_MSQ)

      IF(SHARP) THEN
        SHV = 1.0
      ELSE
        SHV = 0.0
      ENDIF
      WRITE(*,'(A,1X,4(1X,ES24.16))') 'TE',
     &    ANTE, ASTE, DSTE, SHV
      WRITE(*,'(A)', ADVANCE='NO') 'CP'
      DO I=1, N
        WRITE(*,'(1X,ES24.16)', ADVANCE='NO') CPLOC(I)
      END DO
      WRITE(*,*)
      WRITE(*,'(A,1X,5(1X,ES24.16))') 'CL',
     &    CL, CM, CDP, CL_ALF, CL_MSQ

  900 CONTINUE
      END
