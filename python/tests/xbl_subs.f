C***********************************************************************
C  Extracted XFOIL BL helper routines for harness use.
C***********************************************************************
      SUBROUTINE DSLIM(DSTR,THET,UEDG,MSQ,HKLIM)
      IMPLICIT REAL (A-H,M,O-Z)
C
      H = DSTR/THET
      CALL HKIN(H,MSQ,HK,HK_H,HK_M)
C
      DH = MAX( 0.0 , HKLIM-HK ) / HK_H
      DSTR = DSTR + DH*THET
C
      RETURN
      END
