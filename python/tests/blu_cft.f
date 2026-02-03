C***********************************************************************
C  Extracted CFT from BLU for harness use.
C***********************************************************************
      SUBROUTINE CFT( HK, RT, MSQ, CF, CF_HK, CF_RT, CF_MSQ )
      IMPLICIT REAL (A-H,M,O-Z)
      INCLUDE 'BLPAR.INC'
C
      DATA GAM /1.4/
C
C---- Turbulent skin friction function  ( Cf )    (Coles)
      GM1 = GAM - 1.0
      FC = SQRT(1.0 + 0.5*GM1*MSQ)
      GRT = LOG(RT/FC)
      GRT = MAX(GRT,3.0)
C
      GEX = -1.74 - 0.31*HK
C
      ARG = -1.33*HK
      ARG = MAX(-20.0, ARG )
C
      THK = TANH(4.0 - HK/0.875)
C
      CFO =  CFFAC * 0.3*EXP(ARG) * (GRT/2.3026)**GEX
      CF     = ( CFO  +  1.1E-4*(THK-1.0) ) / FC
      CF_HK  = (-1.33*CFO - 0.31*LOG(GRT/2.3026)*CFO
     &         - 1.1E-4*(1.0-THK**2) / 0.875    ) / FC
      CF_RT  = GEX*CFO/(FC*GRT) / RT
      CF_MSQ = GEX*CFO/(FC*GRT) * (-0.25*GM1/FC**2) - 0.25*GM1*CF/FC**2
C
      RETURN
      END ! CFT
