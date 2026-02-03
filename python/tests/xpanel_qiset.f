C***********************************************************************
C  Extracted QISET for xoper parity.
C***********************************************************************
      SUBROUTINE QISET
      INCLUDE 'XFOIL.INC'
C
      COSA = COS(ALFA)
      SINA = SIN(ALFA)
C
      DO 5 I=1, N+NW
        QINV  (I) =  COSA*QINVU(I,1) + SINA*QINVU(I,2)
        QINV_A(I) = -SINA*QINVU(I,1) + COSA*QINVU(I,2)
    5 CONTINUE
C
      RETURN
      END
