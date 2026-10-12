C     Canonical operating-point oracle: original XFOIL numerical routines.
C     Only the command interface and plotting initialization are replaced.
C     A fresh process supplies fresh COMMON storage for every case.
      PROGRAM CANONICAL_REFERENCE
      USE, INTRINSIC :: IEEE_EXCEPTIONS
      INCLUDE 'XFOIL.INC'
      INTEGER IDESIN, NITERIN, IVISCIN, ITYPE
      REAL MACHIN, REIN, ALPHAIN, NCRITIN
      CHARACTER*512 COORDINATES
      LOGICAL INVALID_FLAG, OVERFLOW_FLAG, DIVZERO_FLAG
      CALL INIT
      READ(*,*) IDESIN, MACHIN, REIN, ALPHAIN,
     &          NCRITIN, NITERIN, IVISCIN, CTERAT
      READ(*,'(A)') COORDINATES
      NPAN = 240
      MINF1 = MACHIN
      REINF1 = REIN
      ALFA = ALPHAIN*DTOR
      ADEG = ALPHAIN
      WAKLEN = 1.0
      ACRIT(1) = NCRITIN
      ACRIT(2) = NCRITIN
      LVISC = IVISCIN.EQ.1
      LALFA = .TRUE.
      CALL MRCL(1.0,MINF_CL,REINF_CL)
      CALL COMSET
      IF(IDESIN.GT.0) THEN
        CALL NACA(IDESIN,.FALSE.)
      ELSE
C       LOAD in the XFOIL CLI copies the input nodes without repaneling.
C       This is also how the website handles database/custom coordinates.
        CALL LOAD(COORDINATES,ITYPE)
        IF(ITYPE.EQ.0) STOP 2
        CALL ABCOPY(.TRUE.)
      ENDIF
C     Observe exceptional arithmetic without changing traps or rounding.
      CALL IEEE_SET_FLAG(IEEE_INVALID,.FALSE.)
      CALL IEEE_SET_FLAG(IEEE_OVERFLOW,.FALSE.)
      CALL IEEE_SET_FLAG(IEEE_DIVIDE_BY_ZERO,.FALSE.)
      CALL SPECAL
      IF(LVISC) CALL VISCAL(NITERIN)
      IF(.NOT.LVISC) NW = 0
      CALL IEEE_GET_FLAG(IEEE_INVALID,INVALID_FLAG)
      CALL IEEE_GET_FLAG(IEEE_OVERFLOW,OVERFLOW_FLAG)
      CALL IEEE_GET_FLAG(IEEE_DIVIDE_BY_ZERO,DIVZERO_FLAG)
      WRITE(*,*) 'CANONICAL_FLAGS',INVALID_FLAG,OVERFLOW_FLAG,
     &           DIVZERO_FLAG
      WRITE(*,'(A)') 'CANONICAL_RESULT'
      WRITE(*,*) N,NW,LVCONV
      WRITE(*,'(11ES26.17E3)') CL,CM,CD,CDF,CDP,
     &     XOCTR(1),XOCTR(2),MINF,REINF,ACRIT(1),ACRIT(2)
      WRITE(*,'(3ES26.17E3)') RMSBL,RMXBL,RLX
      IF(LVISC) THEN
        WRITE(*,'(6(I0,1X))') NBL(1),NBL(2),IBLTE(1),IBLTE(2),
     &                       ITRAN(1),ITRAN(2)
        WRITE(*,*) TFORCE(1),TFORCE(2)
      ENDIF
      DO I=1,N+NW
        WRITE(*,'(6ES26.17E3)') X(I),Y(I),QINV(I),QVIS(I),
     &                         CPI(I),CPV(I)
      ENDDO
      IF(LVISC) THEN
        DO IS=1,2
          DO IBL=1,NBL(IS)
            WRITE(*,'(I0,1X,11ES26.17E3)') IPAN(IBL,IS),
     &       XSSI(IBL,IS),UEDG(IBL,IS),THET(IBL,IS),DSTR(IBL,IS),
     &       CTAU(IBL,IS),MASS(IBL,IS),TSTR(IBL,IS),TAU(IBL,IS),
     &       DIS(IBL,IS),CTQ(IBL,IS),DELT(IBL,IS)
          ENDDO
        ENDDO
      ENDIF
      WRITE(*,'(A)') 'CANONICAL_END'
      END

      SUBROUTINE PLINITIALIZE
      END
      SUBROUTINE PLPORT
      END
      SUBROUTINE COLORSPECTRUMHUES(NCOLOR,HUES)
      INTEGER NCOLOR
      CHARACTER*(*) HUES
      END
