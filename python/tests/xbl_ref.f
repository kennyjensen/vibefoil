C     Minimal IBLSYS extracted from XFOIL for test harness.
      SUBROUTINE IBLSYS
      INCLUDE 'XFOIL.INC'
      INCLUDE 'XBL.INC'
      INTEGER IV, IS, IBL
      IV = 0
      DO 10 IS=1, 2
        DO 110 IBL=2, NBL(IS)
          IV = IV+1
          ISYS(IBL,IS) = IV
  110   CONTINUE
   10 CONTINUE
      NSYS = IV
      IF(NSYS.GT.2*IVX) STOP '*** IBLSYS overflow. ***'
      RETURN
      END
