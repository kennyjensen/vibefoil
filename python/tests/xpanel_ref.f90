      SUBROUTINE APCALC
      INCLUDE 'XFOIL.INC'
C
C---- set angles of airfoil panels
      DO 10 I=1, N-1
        SX = X(I+1) - X(I)
        SY = Y(I+1) - Y(I)
        IF(SX.EQ.0.0 .AND. SY.EQ.0.0) THEN
          APANEL(I) = ATAN2( -NY(I) , -NX(I) )
        ELSE
          APANEL(I) = ATAN2( SX , -SY )
        ENDIF
   10 CONTINUE
C
C---- TE panel
      I = N
      IP = 1
      IF(SHARP) THEN
       APANEL(I) = PI
      ELSE
       SX = X(IP) - X(I)
       SY = Y(IP) - Y(I)
       APANEL(I) = ATAN2( -SX , SY ) + PI
      ENDIF
C
      RETURN
      END


      SUBROUTINE NCALC(X,Y,S,N,XN,YN)
C---------------------------------------
C     Calculates normal unit vector
C     components at airfoil panel nodes
C---------------------------------------
      DIMENSION X(N), Y(N), S(N), XN(N), YN(N)
C
      IF(N.LE.1) RETURN
C
      CALL SEGSPL(X,XN,S,N)
      CALL SEGSPL(Y,YN,S,N)
      DO 10 I=1, N
        SX =  YN(I)
        SY = -XN(I)
        SMOD = SQRT(SX*SX + SY*SY)
        IF(SMOD .EQ. 0.0) THEN
         XN(I) = -1.0
         YN(I) = 0.0
        ELSE
         XN(I) = SX/SMOD
         YN(I) = SY/SMOD
        ENDIF
   10 CONTINUE
C
C---- average normal vectors at corner points
      DO 20 I=1, N-1
        IF(S(I) .EQ. S(I+1)) THEN
          SX = 0.5*(XN(I) + XN(I+1))
          SY = 0.5*(YN(I) + YN(I+1))
          SMOD = SQRT(SX*SX + SY*SY)
          IF(SMOD .EQ. 0.0) THEN
           XN(I) = -1.0
           YN(I) = 0.0
           XN(I+1) = -1.0
           YN(I+1) = 0.0
          ELSE
           XN(I)   = SX/SMOD
           YN(I)   = SY/SMOD
           XN(I+1) = SX/SMOD
           YN(I+1) = SY/SMOD
          ENDIF
        ENDIF
 20   CONTINUE
C
      RETURN
      END
