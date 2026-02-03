C***********************************************************************
C  Profil driver for PRWALL/UWALL/FS.
C***********************************************************************
      PROGRAM PROFILDR
      INCLUDE 'BLPAR.INC'
      INTEGER NPRW, NUW, NFS, I, NN, J
      REAL DSTAR, THETA, UO, RT, MS, CT, CFFACI
      REAL DO, UI, CF, BB
      REAL DO_DS, DO_TH, DO_UO, DO_RT, DO_MS
      REAL UI_DS, UI_TH, UI_UO, UI_RT, UI_MS
      REAL HS, HS_DS, HS_TH, HS_UO, HS_RT, HS_MS
      REAL CF_DS, CF_TH, CF_UO, CF_RT, CF_MS
      REAL CD, CD_DS, CD_TH, CD_UO, CD_RT, CD_MS, CD_CT
      REAL TH
      REAL Y(2000), U(2000)
      INTEGER INORM, ISPEC
      REAL BSPEC, HSPEC, ETAE, GEO
      REAL ETA(2000), F(2000), S(2000), DELTA

      READ(*,*,END=900) NPRW, NUW, NFS

      DO I=1, NPRW
        READ(*,*) DSTAR, THETA, UO, RT, MS, CT, CFFACI
        CFFAC = CFFACI
        CALL PRWALL(DSTAR, THETA, UO, RT, MS, CT, BB,
     &              DO, DO_DS, DO_TH, DO_UO, DO_RT, DO_MS,
     &              UI, UI_DS, UI_TH, UI_UO, UI_RT, UI_MS,
     &              HS, HS_DS, HS_TH, HS_UO, HS_RT, HS_MS,
     &              CF, CF_DS, CF_TH, CF_UO, CF_RT, CF_MS,
     &              CD, CD_DS, CD_TH, CD_UO, CD_RT, CD_MS, CD_CT)
        WRITE(*,'(A,1X,4(1X,ES24.16))') 'PRWALL', DO, UI, CF, BB
      END DO

      DO I=1, NUW
        READ(*,*) TH, UO, DO, UI, RT, CF, BB, NN
        CALL UWALL(TH, UO, DO, UI, RT, CF, BB, Y, U, NN)
        WRITE(*,'(A)', ADVANCE='NO') 'UWALL'
        DO J=1, NN
          WRITE(*,'(1X,ES24.16)', ADVANCE='NO') Y(J)
          WRITE(*,'(1X,ES24.16)', ADVANCE='NO') U(J)
        END DO
        WRITE(*,*)
      END DO

      DO I=1, NFS
        READ(*,*) INORM, ISPEC, BSPEC, HSPEC, NN, ETAE, GEO
        CALL FS(INORM, ISPEC, BSPEC, HSPEC, NN, ETAE, GEO,
     &          ETA, F, U, S, DELTA)
        WRITE(*,'(A,1X,ES24.16)') 'FSDELTA', DELTA
        WRITE(*,'(A)', ADVANCE='NO') 'FS'
        DO J=1, NN
          WRITE(*,'(1X,ES24.16)', ADVANCE='NO') ETA(J)
          WRITE(*,'(1X,ES24.16)', ADVANCE='NO') F(J)
          WRITE(*,'(1X,ES24.16)', ADVANCE='NO') U(J)
          WRITE(*,'(1X,ES24.16)', ADVANCE='NO') S(J)
        END DO
        WRITE(*,*)
      END DO

  900 CONTINUE
      END
