C***********************************************************************
C  Xgeom driver for LEFIND.
C***********************************************************************
      PROGRAM XGEOMDR
      INTEGER N, I
      REAL X(10000), Y(10000), S(10000), XP(10000), YP(10000)
      REAL SLE

      READ(*,*,END=900) N
      DO I=1, N
        READ(*,*) X(I)
      END DO
      DO I=1, N
        READ(*,*) Y(I)
      END DO
      DO I=1, N
        READ(*,*) S(I)
      END DO
      DO I=1, N
        READ(*,*) XP(I)
      END DO
      DO I=1, N
        READ(*,*) YP(I)
      END DO

      CALL LEFIND(SLE, X, XP, Y, YP, S, N)
      WRITE(*,'(ES24.16)') SLE

  900 CONTINUE
      END
