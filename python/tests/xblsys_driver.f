C***********************************************************************
C  XBLsys function driver for JS/Fortran parity.
C***********************************************************************
      PROGRAM XBLSYSDR
      INTEGER NHKIN, NDIL, NDILW, NHSL, NCFL, NDIT, NHST, NCFT, NHCT
      INTEGER NDSLIM, I
      REAL H, MSQ, HK, HK_H, HK_MSQ
      REAL RT, DI, DI_HK, DI_RT
      REAL HS, HS_HK, HS_RT, HS_MSQ
      REAL CF, CF_HK, CF_RT, CF_MSQ
      REAL US, ST, DI_US, DI_CF, DI_ST, DI_HS
      REAL CFFAC
      REAL HC, HC_HK, HC_MSQ
      REAL DSTR, THET, UEDG, HKLIM

      READ(*,*,END=900) NHKIN, NDIL, NDILW, NHSL, NCFL, NDIT,
     &                 NHST, NCFT, NHCT, NDSLIM

      DO I=1, NHKIN
        READ(*,*) H, MSQ
        CALL HKIN(H, MSQ, HK, HK_H, HK_MSQ)
        WRITE(*,'(A,1X,3(1X,ES24.16))') 'HKIN', HK, HK_H, HK_MSQ
      END DO

      DO I=1, NDIL
        READ(*,*) HK, RT
        CALL DIL(HK, RT, DI, DI_HK, DI_RT)
        WRITE(*,'(A,1X,3(1X,ES24.16))') 'DIL', DI, DI_HK, DI_RT
      END DO

      DO I=1, NDILW
        READ(*,*) HK, RT
        CALL DILW(HK, RT, DI, DI_HK, DI_RT)
        WRITE(*,'(A,1X,3(1X,ES24.16))') 'DILW', DI, DI_HK, DI_RT
      END DO

      DO I=1, NHSL
        READ(*,*) HK, RT, MSQ
        CALL HSL(HK, RT, MSQ, HS, HS_HK, HS_RT, HS_MSQ)
        WRITE(*,'(A,1X,4(1X,ES24.16))') 'HSL', HS, HS_HK, HS_RT, HS_MSQ
      END DO

      DO I=1, NCFL
        READ(*,*) HK, RT, MSQ
        CALL CFL(HK, RT, MSQ, CF, CF_HK, CF_RT, CF_MSQ)
        WRITE(*,'(A,1X,4(1X,ES24.16))') 'CFL', CF, CF_HK, CF_RT, CF_MSQ
      END DO

      DO I=1, NDIT
        READ(*,*) HS, US, CF, ST
        CALL DIT(HS, US, CF, ST, DI, DI_HS, DI_US, DI_CF, DI_ST)
        WRITE(*,'(A,1X,5(1X,ES24.16))') 'DIT', DI, DI_HS, DI_US,
     &        DI_CF, DI_ST
      END DO

      DO I=1, NHST
        READ(*,*) HK, RT, MSQ
        CALL HST(HK, RT, MSQ, HS, HS_HK, HS_RT, HS_MSQ)
        WRITE(*,'(A,1X,4(1X,ES24.16))') 'HST', HS, HS_HK, HS_RT, HS_MSQ
      END DO

      DO I=1, NCFT
        READ(*,*) HK, RT, MSQ, CFFAC
        CALL CFT(HK, RT, MSQ, CFFAC, CF, CF_HK, CF_RT, CF_MSQ)
        WRITE(*,'(A,1X,4(1X,ES24.16))') 'CFT', CF, CF_HK, CF_RT, CF_MSQ
      END DO

      DO I=1, NHCT
        READ(*,*) HK, MSQ
        CALL HCT(HK, MSQ, HC, HC_HK, HC_MSQ)
        WRITE(*,'(A,1X,3(1X,ES24.16))') 'HCT', HC, HC_HK, HC_MSQ
      END DO

      DO I=1, NDSLIM
        READ(*,*) DSTR, THET, UEDG, MSQ, HKLIM
        CALL DSLIM(DSTR, THET, UEDG, MSQ, HKLIM)
        WRITE(*,'(A,1X,ES24.16)') 'DSLIM', DSTR
      END DO

  900 CONTINUE
      END
