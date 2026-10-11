C     Headless operating-point oracle. All numerical routines are linked
C     directly from third_party/Xfoil/src; only plotting is stubbed.
      PROGRAM VISCAL_REFERENCE
      INCLUDE 'XFOIL.INC'
      INTEGER IDESIN, NPANIN, NITERIN, ISHARP
      REAL MACHIN, REIN, ALPHAIN, CTERIN, WAKEIN
      CALL INIT
      READ(*,*) IDESIN, NPANIN, CTERIN, MACHIN, REIN,
     &          ALPHAIN, WAKEIN, NITERIN, ISHARP
      NPAN = NPANIN
      CTERAT = CTERIN
      MINF1 = MACHIN
      REINF1 = REIN
      ALFA = ALPHAIN*DTOR
      ADEG = ALPHAIN
      WAKLEN = WAKEIN
      LVISC = .TRUE.
      LALFA = .TRUE.
      CALL MRCL(1.0,MINF_CL,REINF_CL)
      CALL COMSET
      CALL NACA(IDESIN,ISHARP.NE.0)
      CALL SPECAL
      IF(NITERIN.GT.0) CALL VISCAL(NITERIN)
      WRITE(*,'(A)') 'REFERENCE_RESULT'
      WRITE(*,'(9ES26.17E3)') ALFA/DTOR,CL,CM,CD,CDF,CDP,CLF,
     &                      RMSBL,RMXBL
      WRITE(*,'(3ES26.17E3)') RLX,XOCTR(1),XOCTR(2)
      WRITE(*,'(6(I0,1X))') ITRAN(1),ITRAN(2),IMXBL,ISMXBL,
     &                     NBL(1),NBL(2)
      WRITE(*,*) TFORCE(1),TFORCE(2),LVCONV,VMXBL
      WRITE(*,'(2(I0,1X))') N,NW
      DO I=1,N+NW
        WRITE(*,'(5ES26.17E3)') X(I),Y(I),QINV(I),QVIS(I),CPV(I)
      ENDDO
      DO IS=1,2
        DO IBL=1,NBL(IS)
          WRITE(*,'(10ES26.17E3)') UEDG(IBL,IS),THET(IBL,IS),
     &      DSTR(IBL,IS),CTAU(IBL,IS),MASS(IBL,IS),TSTR(IBL,IS),
     &      DSSR(IBL,IS),MDEF(IBL,IS),PDEF(IBL,IS),KDEF(IBL,IS)
        ENDDO
      ENDDO
      END

C     INIT's graphical setup has no effect on numerical solver state.
      SUBROUTINE PLINITIALIZE
      END
      SUBROUTINE PLPORT
      END
      SUBROUTINE COLORSPECTRUMHUES(NCOLOR,HUES)
      INTEGER NCOLOR
      CHARACTER*(*) HUES
      END
