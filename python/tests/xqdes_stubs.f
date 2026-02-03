C***********************************************************************
C  Minimal stubs to satisfy XFOIL/XQDES plotting & UI symbols.
C***********************************************************************
      SUBROUTINE PLOT(X,Y,IPEN)
      RETURN
      END

      SUBROUTINE NEWPEN(ICOL)
      RETURN
      END

      SUBROUTINE PLTINI
      RETURN
      END

      SUBROUTINE PLNUMB(X,Y,H,VAL,ANG,NDEC)
      RETURN
      END

      SUBROUTINE PLCHAR(X,Y,H,STRING,ANGLE,NCHAR)
      CHARACTER*(*) STRING
      RETURN
      END

      SUBROUTINE PLMATH(X,Y,H,STRING,ANGLE,NCHAR)
      CHARACTER*(*) STRING
      RETURN
      END

      SUBROUTINE PLGRID(X1,X2,Y1,Y2,DX,DY)
      RETURN
      END

      SUBROUTINE PLFLUSH
      RETURN
      END

      SUBROUTINE DASH(LDASH)
      RETURN
      END

      SUBROUTINE PLSYMB(X,Y,SIZE,ISYM)
      RETURN
      END

      SUBROUTINE GETCOLOR(ICOL)
      RETURN
      END

      SUBROUTINE NEWCOLORNAME(ICOL,NAME)
      CHARACTER*(*) NAME
      RETURN
      END

      SUBROUTINE NEWCOLOR(ICOL,R,G,B)
      RETURN
      END

      SUBROUTINE GETCURSORXY(X,Y,CH)
      CHARACTER*1 CH
      RETURN
      END

      SUBROUTINE KEYOFF
      RETURN
      END

      SUBROUTINE OFFGET
      RETURN
      END

      SUBROUTINE ANNOT
      RETURN
      END

      SUBROUTINE SPECAL
      RETURN
      END

      SUBROUTINE MODIFY
      RETURN
      END

      SUBROUTINE REPlot
      RETURN
      END

      SUBROUTINE USETZOOM
      RETURN
      END

      SUBROUTINE CLRZOOM
      RETURN
      END

      SUBROUTINE PLEND
      RETURN
      END

      SUBROUTINE PLINITIALIZE
      RETURN
      END

      SUBROUTINE COLORSPECTRUMHUES
      RETURN
      END

      SUBROUTINE PANPLT
      RETURN
      END

      SUBROUTINE CANG
      RETURN
      END

      SUBROUTINE HALF
      RETURN
      END

      SUBROUTINE BENDUMP
      RETURN
      END

      SUBROUTINE BENDUMP2
      RETURN
      END

      SUBROUTINE OPLSET(IDEV,IDEVRP,IPSLU,
     &                  SIZE,PAR,
     &                  XMARG,YMARG,XPAGE,YPAGE,
     &                  CSIZE,SCRNFR,LCURS,LLAND, ICOLS)
      INTEGER ICOLS(2)
      LOGICAL LCURS,LLAND
      RETURN
      END

      SUBROUTINE NAMMOD
      RETURN
      END

      SUBROUTINE BLPINI
      RETURN
      END

      SUBROUTINE IBLSYS
      RETURN
      END
