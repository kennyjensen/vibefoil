C     Minimal MRCL/COMSET extracted from XFOIL for test harness.
      SUBROUTINE MRCL(CLS,M_CLS,R_CLS)
      INCLUDE 'XFOIL.INC'
      REAL M_CLS
C
      CLA = MAX( CLS , 0.000001 )
C
      IF(RETYP.LT.1 .OR. RETYP.GT.3) THEN
        RETYP = 1
      ENDIF
      IF(MATYP.LT.1 .OR. MATYP.GT.3) THEN
        MATYP = 1
      ENDIF
C
      IF(MATYP.EQ.1) THEN
        MINF  = MINF1
        M_CLS = 0.
      ELSE IF(MATYP.EQ.2) THEN
        MINF  =  MINF1/SQRT(CLA)
        M_CLS = -0.5*MINF/CLA
      ELSE IF(MATYP.EQ.3) THEN
        MINF  = MINF1
        M_CLS = 0.
      ENDIF
C
      IF(RETYP.EQ.1) THEN
        REINF = REINF1
        R_CLS = 0.
      ELSE IF(RETYP.EQ.2) THEN
        REINF =  REINF1/SQRT(CLA)
        R_CLS = -0.5*REINF/CLA
      ELSE IF(RETYP.EQ.3) THEN
        REINF =  REINF1/CLA
        R_CLS = -REINF /CLA
      ENDIF
C
      IF(MINF .GE. 0.99) THEN
        MINF = 0.99
        M_CLS = 0.
      ENDIF
C
      RRAT = 1.0
      IF(REINF1 .GT. 0.0) RRAT = REINF/REINF1
C
      IF(RRAT .GT. 100.0) THEN
        REINF = REINF1*100.0
        R_CLS = 0.
      ENDIF
C
      RETURN
      END


      SUBROUTINE COMSET
      INCLUDE 'XFOIL.INC'
C
      BETA = SQRT(1.0 - MINF**2)
      BETA_MSQ = -0.5/BETA
C
      TKLAM   = MINF**2 / (1.0 + BETA)**2
      TKL_MSQ =     1.0 / (1.0 + BETA)**2
     &    - 2.0*TKLAM/ (1.0 + BETA) * BETA_MSQ
C
      IF(MINF.EQ.0.0) THEN
        CPSTAR = -999.0
        QSTAR = 999.0
      ELSE
        CPSTAR = 2.0 / (GAMMA*MINF**2)
     &        * (((1.0 + 0.5*GAMM1*MINF**2)
     &            /(1.0 + 0.5*GAMM1        ))**(GAMMA/GAMM1) - 1.0)
        QSTAR = QINF/MINF
     &       * SQRT( (1.0 + 0.5*GAMM1*MINF**2)
     &              /(1.0 + 0.5*GAMM1        ) )
      ENDIF
C
      RETURN
      END
