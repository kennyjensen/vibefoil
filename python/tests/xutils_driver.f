C***********************************************************************
C  Xutils driver for ATANC and SETEXP.
C***********************************************************************
      PROGRAM XUTILSDR
      INTEGER NATANC, NSETEXP, I, NN, J
      REAL Y, X, THOLD, THNEW
      REAL DS1, SMAX
      REAL S(10000)

      READ(*,*,END=900) NATANC, NSETEXP

      DO I=1, NATANC
        READ(*,*) Y, X, THOLD
        THNEW = ATANC(Y, X, THOLD)
        WRITE(*,'(A,1X,ES24.16)') 'ATANC', THNEW
      END DO

      DO I=1, NSETEXP
        READ(*,*) DS1, SMAX, NN
        CALL SETEXP(S, DS1, SMAX, NN)
        WRITE(*,'(A)', ADVANCE='NO') 'SETEXP'
        DO J=1, NN
          WRITE(*,'(1X,ES24.16)', ADVANCE='NO') S(J)
        END DO
        WRITE(*,*)
      END DO

  900 CONTINUE
      END
