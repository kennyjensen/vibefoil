C***********************************************************************
C  Extracted XFOIL GDES helper routines for harness use.
C***********************************************************************
      LOGICAL FUNCTION INSIDE(X,Y,N, XF,YF)
      DIMENSION X(N),Y(N)
C-------------------------------------
C     Returns .TRUE. if point XF,YF 
C     is inside contour X(i),Y(i).
C-------------------------------------
      ANGLE = 0.0
      DO 10 I=1, N
        IP = I+1
        IF(I.EQ.N) IP = 1
        XB1 = X(I)  - XF
        YB1 = Y(I)  - YF
        XB2 = X(IP) - XF
        YB2 = Y(IP) - YF
        ANGLE = ANGLE + (XB1*YB2 - YB1*XB2)
     &                   / SQRT((XB1**2 + YB1**2)*(XB2**2 + YB2**2))
 10   CONTINUE
      INSIDE = ABS(ANGLE) .GT. 1.0
      RETURN
      END ! INSIDE


      SUBROUTINE GETXYF(X,XP,Y,YP,S,N, TOPS,BOTS,XF,YF)
      DIMENSION X(N),XP(N),Y(N),YP(N),S(N)
C
      IF(XF .EQ. -999.0)
     &  CALL ASKR('Enter flap hinge x location^',XF)
C
C---- find top and bottom y at hinge x location
      TOPS = S(1) + (X(1) - XF)
      BOTS = S(N) - (X(N) - XF)
      CALL SINVRT(TOPS,XF,X,XP,S,N)      
      CALL SINVRT(BOTS,XF,X,XP,S,N)      
      TOPY = SEVAL(TOPS,Y,YP,S,N)
      BOTY = SEVAL(BOTS,Y,YP,S,N)
C
      WRITE(*,1000) TOPY, BOTY
 1000 FORMAT(/'  Top    surface:  y =', F8.4,'     y/t = 1.0'
     &       /'  Bottom surface:  y =', F8.4,'     y/t = 0.0')
C
      IF(YF .EQ. -999.0)
     & CALL ASKR(
     &  'Enter flap hinge y location (or 999 to specify y/t)^',YF)
C
      IF(YF .EQ. 999.0) THEN
        CALL ASKR('Enter flap hinge relative y/t location^',YREL)
        YF = TOPY*YREL + BOTY*(1.0-YREL)
      ENDIF
C
      RETURN
      END ! GETXYF
