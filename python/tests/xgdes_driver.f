C***********************************************************************
C  Xgdes driver for INSIDE, GETXYF, and SSS.
C***********************************************************************
      PROGRAM XGDESDR
      INTEGER N, NINS, NSSS, I, ISIDE
      REAL X(10000), Y(10000), S(10000), XP(10000), YP(10000)
      REAL TOPS, BOTS, XF, YF
      REAL SS, DEL, XBF, YBF, S1, S2
      LOGICAL INSIDE

      READ(*,*,END=900) N, NINS, NSSS
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

      READ(*,*) XF, YF
      CALL GETXYF(X,XP,Y,YP,S,N, TOPS,BOTS,XF,YF)
      WRITE(*,'(A,1X,4(1X,ES24.16))') 'GETXYF', TOPS, BOTS, XF, YF

      DO I=1, NINS
        READ(*,*) XBF, YBF
        IF(INSIDE(X,Y,N,XBF,YBF)) THEN
          WRITE(*,'(A,1X,I1)') 'INSIDE', 1
        ELSE
          WRITE(*,'(A,1X,I1)') 'INSIDE', 0
        ENDIF
      END DO

      DO I=1, NSSS
        READ(*,*) SS, DEL, XBF, YBF, ISIDE
        CALL SSS(SS,S1,S2,DEL,XBF,YBF,X,XP,Y,YP,S,N,ISIDE)
        WRITE(*,'(A,1X,2(1X,ES24.16))') 'SSS', S1, S2
      END DO

  900 CONTINUE
      END
