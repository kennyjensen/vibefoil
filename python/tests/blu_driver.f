C***********************************************************************
C  BLU driver for CFT.
C***********************************************************************
      PROGRAM BLUDR
      INTEGER NCFT, I
      REAL HK, RT, MSQ, CF, CF_HK, CF_RT, CF_MSQ

      CALL BLPINI

      READ(*,*,END=900) NCFT
      DO I=1, NCFT
        READ(*,*) HK, RT, MSQ
        CALL CFT(HK, RT, MSQ, CF, CF_HK, CF_RT, CF_MSQ)
        WRITE(*,'(A,1X,4(1X,ES24.16))') 'CFT', CF, CF_HK, CF_RT, CF_MSQ
      END DO

  900 CONTINUE
      END
