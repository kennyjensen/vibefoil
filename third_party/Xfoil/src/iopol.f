C***********************************************************************
C    Module:  iopol.f
C 
C    Copyright (C) 2000 Mark Drela, Harold Youngren
C 
C    This program is free software; you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation; either version 2 of the License, or
C    (at your option) any later version.
C
C    This program is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with this program; if not, write to the Free Software
C    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
C***********************************************************************

      SUBROUTINE POLREAD(LU,FNPOL,ERROR,
     &                   NAX,NA,CPOL, 
     &                   REYN1,MACH1,ACRIT,XTRIP,
     &                   PTRAT,ETAP,
     &                   NAME, IRETYP,IMATYP,
     &                   NBX,NBL,CPOLSD,CPOLEL,
     &                   CODE,VERSION )
      INCLUDE 'PINDEX.INC'
      CHARACTER*(*) FNPOL, NAME
      LOGICAL ERROR
      CHARACTER*(*) CODE
      REAL CPOL(NAX,IPTOT),
     &     CPOLSD(NAX,2*NBX,JPTOT),
     &     CPOLEL(NAX,NBX,KPTOT)
      REAL MACH1, ACRIT(2*NBX), XTRIP(2*NBX)
C--------------------------------------------------------
C     Reads in polar save file
C
C  Input:
C     LU      logical unit to use for reading
C     FNPOL   name of polar file to be read,
C               if FNPOL(1:1).eq.' ', unit LU will be read
C               if it is already open
C     NAX     polar point array dimension
C     NBX     airfoil element array dimension
C     NBL     number of airfoil elements
C
C  Output:
C     ERROR   T if a READ error occurred
C     NA      number polar points
C     CPOL    polar coefficients and parameters
C     REYN1   Reynolds number for CL=1
C     MACH1   Mach number for CL=1
C     ACRIT   Critical amplification ratio
C     XTRIP   Trip locations
C     PTRAT   Actuator disk total-pressure ratio
C     ETAP    Actuator disk thermal efficiency
C     NAME    airfoil name string
C     IRETYP  flag giving type of Re variation with CL
C     IMATYP  flag giving type of Ma variation with CL
C     NBL     number of airfoil elements
C     CPOLSD  airfoil side-related parameters
C     CPOLEL  airfoil element-related parameters
C     CODE    code used to compute polar
C     VERSION code version
C--------------------------------------------------------
      CHARACTER*256 LINE
      CHARACTER*1 DUMMY

      PARAMETER (NBXL=5)
      REAL RINP(0:IPTOT + 2*NBXL*JPTOT + NBXL*KPTOT)
C
      INTEGER IPOL(IPTOT), ISPOL(2,JPTOT), KPOL(KPTOT)
      INTEGER ITMP(IPTOT+2*JPTOT), ITMP0(IPTOT+2*JPTOT)
      LOGICAL LOPEN, LHEAD, LDLAB
      LOGICAL LIRE, LIMA, LJNC, LJTP
      CHARACTER*20 CPNAME
C
      LLEN = LEN(LINE)
C
      ERROR = .FALSE.
      LHEAD = .TRUE.
C
      NA = 0
      NBL = 1
C
c      KCH = 0
c      KMC = 0
C
      NIPOL = 0
      DO IP = 1, IPTOT
        IPOL(IP) = 0
      ENDDO
      DO JP = 1, JPTOT
        ISPOL(1,JP) = 0
        ISPOL(2,JP) = 0
      ENDDO
      DO KP = 1, KPTOT
        KPOL(KP) = 0
      ENDDO

C---- assume Re,Mach will not be given in header
      IRETYP = 0
      IMATYP = 0
C
C---- do we have to open the file?
      LOPEN = FNPOL .NE. ' '
C
      IF(LOPEN) OPEN(LU,FILE=FNPOL,STATUS='OLD',ERR=90)
C
C=============================================================
C---- start data reading loop
 500  CONTINUE
      READ(LU,1000,END=80) LINE
      IF(LINE.EQ.' ') GO TO 500
C
      IF(LHEAD) THEN
C----- parse to get header info
C
C----- assume this will be the data-label line
       LDLAB = .TRUE.
C
C--------------------------------------------
       K = INDEX(LINE,'Version')
       IF(K.NE.0) THEN
C------ code,version line
        DO K1=1, LLEN
          IF(LINE(K1:K1).NE.' ') GO TO 10
        ENDDO
C
 10     CONTINUE
        IF(K.GT.K1) THEN
         CODE = LINE(K1:K-1)
         READ(LINE(K+7:LLEN),*,ERR=11) VERSION
        ENDIF
 11     CONTINUE
        LDLAB = .FALSE.
       ENDIF
C
C--------------------------------------------
       KF = INDEX(LINE,'for:')
       IF(KF.NE.0) THEN
C------ airfoil name line
        NAME = LINE(KF+5:LLEN)
        LDLAB = .FALSE.
       ENDIF
C
C--------------------------------------------
       KE = INDEX(LINE,'elements')
       IF(KE.GT.0) THEN
C------ element-number line
        READ(LINE(KE-4:KE-1),*,ERR=60) NBL
C------ truncate name line to eliminate elements #
        NAME = LINE(KF+5:KE-4)
C
        IF(NBL .GT. NBX) THEN
         NBL = NBX
         WRITE(*,*) 
     &   'POLREAD: Number of elements set to array limit', NBL
        ENDIF
        LDLAB = .FALSE.
       ENDIF
C
C--------------------------------------------
       KR = INDEX(LINE,'Reynolds number')
       KM = INDEX(LINE,'Mach number')
C
       IF(KR.NE.0) THEN
C------ Re-type line
        IF(KM.GT.KR) THEN
         KEND = KM-1
        ELSE
         KEND = LLEN
        ENDIF
        IF    (INDEX(LINE(KR:KEND),'fixed').NE.0) THEN
         IRETYP = 1
        ELSEIF(INDEX(LINE(KR:KEND),'1/sqrt(CL)').NE.0) THEN
         IRETYP = 2
        ELSEIF(INDEX(LINE(KR:KEND),'1/CL').NE.0) THEN
         IRETYP = 3
        ENDIF
        LDLAB = .FALSE.
       ENDIF
C
       IF(KM.NE.0) THEN
C------ Ma-type line
        IF(KR.GT.KM) THEN
         KEND = KR-1
        ELSE
         KEND = LLEN
        ENDIF
        IF    (INDEX(LINE(KM:KEND),'fixed').NE.0) THEN
         IMATYP = 1
        ELSEIF(INDEX(LINE(KM:KEND),'1/sqrt(CL)').NE.0) THEN
         IMATYP = 2
        ELSEIF(INDEX(LINE(KM:KEND),'1/CL').NE.0) THEN
         IMATYP = 3
        ENDIF
        LDLAB = .FALSE.
       ENDIF
C
C--------------------------------------------
C---- find specified BL trip location
       K = INDEX(LINE,'xtrf')
       IF(K.NE.0) THEN
C------ new style xtrip line
        KT = INDEX(LINE,'(top)')
        KB = INDEX(LINE,'(bottom)')
        KE = INDEX(LINE,'element ')
C--- check for old style trip line 
        KS = INDEX(LINE,'(suc')
        KP = INDEX(LINE,'(pre')
C
        IF(KE.NE.0) THEN
         READ(LINE(KE+7:KE+12),*,ERR=21) N
        ELSE
         N = 1
        ENDIF
        IF(N.LE.NBL) THEN
         IS1 = 2*N-1
         IS2 = 2*N
         XTRIP(IS1) = 1.0
         XTRIP(IS2) = 1.0
         IF(KT.GT.0)  READ(LINE(K+6:KT-1) ,*,ERR=21) XTRIP(IS1)
         IF(KB.GT.KT) READ(LINE(KT+5:KB-1),*,ERR=21) XTRIP(IS2)
         IF(KS.GT.0)  READ(LINE(K+6:KS-1) ,*,ERR=21) XTRIP(IS1)
         IF(KP.GT.KS) READ(LINE(KS+5:KP-1),*,ERR=21) XTRIP(IS2)
        ENDIF
 21     CONTINUE
        LDLAB = .FALSE.
       ENDIF
C
C--------------------------------------------
       K = INDEX(LINE,'Mach =')
       IF(K.NE.0) THEN
        READ(LINE(K+6:LLEN),*,ERR=31) MACH1
 31     CONTINUE
        LDLAB = .FALSE.
       ENDIF
C
C--------------------------------------------
       K = INDEX(LINE,'Re =')
       IF(K.NE.0) THEN
        READ(LINE(K+4:LLEN),*,ERR=32) REYN1
        REYN1 = REYN1 * 1.0E6
 32     CONTINUE
        LDLAB = .FALSE.
       ENDIF
C
C--------------------------------------------
       K = INDEX(LINE,'Ncrit =')
       IF(K.NE.0) THEN
        NINP = 2
        CALL GETFLT(LINE(k+7:LLEN),RINP(1),NINP,ERROR)
        IF(NINP.LE.0 .OR. ERROR) GO TO 33
C
        IF(NINP.EQ.1) THEN
         ACRIT(1) = RINP(1)
         ACRIT(2) = RINP(1)
        ELSE
         ACRIT(1) = RINP(1)
         ACRIT(2) = RINP(2)
        ENDIF
 33     CONTINUE
        LDLAB = .FALSE.
       ENDIF
C
C--------------------------------------------
       K = INDEX(LINE,'pi_p =')
       IF(K.NE.0) THEN
        READ(LINE(K+6:LLEN),*,ERR=34) PTRAT
 34     CONTINUE
        LDLAB = .FALSE.
       ENDIF
C
C--------------------------------------------
       K = INDEX(LINE,'eta_p =')
       IF(K.NE.0) THEN
        READ(LINE(K+7:LLEN),*,ERR=35) ETAP
 35     CONTINUE
        LDLAB = .FALSE.
       ENDIF
C
C--------------------------------------------
       IF(LDLAB .AND. NIPOL.EQ.0) THEN
C------ process line for possible data labels
        DO IP = 1, IPTOT
          CALL STRIP(CPOLNAME(IP),NNAME)
C
C-------- mark this parameter for reading
          K = INDEX(LINE,CPOLNAME(IP)(1:NNAME))
          ITMP0(IP) = K
          ITMP(IP)  = K
        ENDDO
C
        DO JP = 1, JPTOT
          CALL STRIP(CPOLSNAME(JP),NNAME)
C
          CPNAME = 'Top ' // CPOLSNAME(JP)
          K1 = INDEX(LINE,CPNAME(1:NNAME+4))
          CPNAME = 'Top_' // CPOLSNAME(JP)
          K2 = INDEX(LINE,CPNAME(1:NNAME+4))
          ITMP0(IPTOT+JP) = MAX(K1,K2)
          ITMP (IPTOT+JP) = MAX(K1,K2)
C
          CPNAME = 'Bot ' // CPOLSNAME(JP)
          K1 = INDEX(LINE,CPNAME(1:NNAME+4))
          CPNAME = 'Bot_' // CPOLSNAME(JP)
          K2 = INDEX(LINE,CPNAME(1:NNAME+4))
          ITMP0(IPTOT+JP+JPTOT) = MAX(K1,K2)
          ITMP (IPTOT+JP+JPTOT) = MAX(K1,K2)
        ENDDO
C
C------ bubble-sort data label positions in line string
        DO IPASS = 1, IPTOT+2*JPTOT
          DO IP = 1, IPTOT+2*JPTOT-1
            IF(ITMP(IP).GT.ITMP(IP+1)) THEN
             ITMPP1 = ITMP(IP+1)
             ITMP(IP+1) = ITMP(IP)
             ITMP(IP) = ITMPP1
            ENDIF
          ENDDO
        ENDDO
C
C------ assign data position to each parameter
        DO IPT = 1, IPTOT + 2*JPTOT
          IF(ITMP(IPT).GT.0) THEN
           NIPOL = NIPOL + 1
           DO IP = 1, IPTOT
             IF(ITMP(IPT).EQ.ITMP0(IP)) IPOL(IP) = NIPOL
           ENDDO
           DO JP = 1, JPTOT
             IF(ITMP(IPT).EQ.ITMP0(IPTOT+JP      )) ISPOL(1,JP) = NIPOL
             IF(ITMP(IPT).EQ.ITMP0(IPTOT+JPTOT+JP)) ISPOL(2,JP) = NIPOL
           ENDDO
          ENDIF
        ENDDO
C
       ENDIF
C
C--------------------------------------------
       IF(INDEX(LINE,'-----').NE.0) THEN
        LHEAD = .FALSE.
       ENDIF
C
C--------------------------------------------------------------
      ELSE
C----- read polar data lines
       IA = NA + 1
C
       NINP = IPTOT + 2*NBL*JPTOT + NBL*KPTOT
       DO IINP = 1, NINP
         RINP(IINP) = 0.
       ENDDO
       CALL GETFLT(LINE,RINP(1),NINP,ERROR)
       IF(ERROR) GO TO 90
C
       DO IP = 1, IPTOT
         IF(IPOL(IP) .GT. 0) THEN
          CPOL(IA,IP) = RINP(IPOL(IP))
         ELSE
          CPOL(IA,IP) = 0.
         ENDIF
       ENDDO
C
       DO JP = 1, JPTOT
         DO N = 1, NBL
           IS1 = 2*N-1
           IS2 = 2*N
           IF(ISPOL(1,JP) .GT. 0) THEN
            CPOLSD(IA,IS1,JP) = RINP(ISPOL(1,JP)+2*(N-1))
           ELSE
            CPOLSD(IA,IS1,JP) = 0.
           ENDIF
           IF(ISPOL(2,JP) .GT. 0) THEN
            CPOLSD(IA,IS2,JP) = RINP(ISPOL(2,JP)+2*(N-1))
           ELSE
            CPOLSD(IA,IS2,JP) = 0.
           ENDIF
         ENDDO
       ENDDO
C
       DO KP = 1, KPTOT
         DO N = 1, NBL
           IF(KPOL(KP) .GT. 0) THEN
            CPOLEL(IA,N,KP) = RINP(KPOL(KP))
           ELSE
            CPOLEL(IA,N,KP) = 0.
           ENDIF
         ENDDO
       ENDDO

       ACL = MAX( CPOL(IA,ICL) , 0.001 )
C
C
C----- try to find Re, Ma, Ncrit, Xtrip  in polar data
       LIRE = .FALSE.
       LIMA = .FALSE.
       LJNC = .FALSE.
       LJTP = .FALSE.
       DO KP = 1, NIPOL
         IF(IPOL(KP) .EQ. IRE) LIRE = .TRUE.
         IF(IPOL(KP) .EQ. IMA) LIMA = .TRUE.
c         IF(ISPOL(1,KP) .EQ. JNC) LJNC = .TRUE.
c         IF(ISPOL(1,KP) .EQ. JTP) LJTP = .TRUE.
       ENDDO
C
       IF(.NOT. LIRE) THEN
C------ Re was not in polar data... set using header info
        IF    (IRETYP.EQ.1) THEN
         CPOL(IA,IRE) = REYN1
        ELSEIF(IRETYP.EQ.2) THEN
         CPOL(IA,IRE) = REYN1/SQRT(ACL) 
        ELSEIF(IRETYP.EQ.3) THEN
        CPOL(IA,IRE) = REYN1/ACL
        ENDIF
       ENDIF
C
       IF(.NOT. LIMA) THEN
C------ Mach was not in polar data... set using header info
        IF    (IMATYP.EQ.1) THEN
         CPOL(IA,IMA) = MACH1
        ELSEIF(IMATYP.EQ.2) THEN
         CPOL(IA,IMA) = MACH1/SQRT(ACL)
        ELSEIF(IMATYP.EQ.3) THEN
         CPOL(IA,IMA) = MACH1/ACL
        ENDIF
       ENDIF
C
       IF(.NOT. LJNC) THEN
C------ Ncrit was not in polar data... set using header info
        DO IS = 1, 2*NBL
          CPOLSD(IA,IS,JNC) = ACRIT(IS)
        ENDDO
       ENDIF
C
       IF(.NOT. LJTP) THEN
C------ set trip data using header info
        DO IS = 1, 2*NBL
          CPOLSD(IA,IS,JTP) = XTRIP(IS)
        ENDDO
       ENDIF
C
       NA = IA
      ENDIF
C
 60   CONTINUE
C---- go read next line
      GO TO 500
C=============================================================
C
 80   CONTINUE
C---- if file was opened here, then close it
      IF(LOPEN) CLOSE(LU)
      RETURN
C
 90   CONTINUE
      IF(LOPEN) CLOSE(LU)
      ERROR = .TRUE.
      RETURN
C
C..........................................
 1000 FORMAT(A)
      END ! POLREAD


      SUBROUTINE POLWRIT(LU,FNPOL,ERROR, LHEAD,
     &                   NAX, IA1,IA2, CPOL, IPOL,NIPOL,
     &                   REYN1,MACH1,ACRIT,XTRIP,
     &                   PTRAT,ETAP,
     &                   NAME, IRETYP,IMATYP,
     &                   NBX,NBL,CPOLSD,JPOL,NJPOL, 
     &                           CPOLEL,KPOL,NKPOL,
     &                   CODE,VERSION, LQUERY )
      INCLUDE 'PINDEX.INC'
      CHARACTER*(*) FNPOL, NAME
      LOGICAL ERROR, LHEAD,LQUERY
      CHARACTER*(*) CODE
      REAL CPOL(NAX,IPTOT),
     &     CPOLSD(NAX,2*NBX,JPTOT),
     &     CPOLEL(NAX,NBX,KPTOT)
      REAL MACH1, ACRIT(2*NBX), XTRIP(2*NBX)
      INTEGER IPOL(IPTOT), JPOL(JPTOT), KPOL(KPTOT)
C--------------------------------------------------------
C     Writes polar save file
C
C  Input:
C     LU       logical unit to use for writing
C     FNPOL    name of polar file to be read,
C                if FNPOL(1:1).eq.' ', unit LU is assumed 
C                to be already open
C     NAX      polar point array dimension
C     NBX      airfoil element array dimension
C     NBL      number of airfoil elements
C     IA1,IA2  only polar points IA1..IA2 are written
C     CPOL     polar coefficients and parameters
C     IPOL(.)  indices of data quantities to be written
C     NIPOL    number  of data quantities to be written
C     REYN1    Reynolds number for CL=1
C     MACH1    Mach number for CL=1
C     ACRIT    Critical amplification ratio
C     XTRIP    Trip locations
C     PTRAT   Actuator disk total-pressure ratio
C     ETAP    Actuator disk thermal efficiency
C     NAME     airfoil name string
C     IRETYP   flag giving type of Re variation with CL
C     IMATYP   flag giving type of Ma variation with CL
C     NBL      number of airfoil elements
C     CPOLSD   airfoil side-related parameters
C     JPOL(.)  indices of side data quantities to be written
C     NJPOL    number  of side data quantities to be written
C     CPOLEL   airfoil element-related parameters
C     KPOL(.)  indices of element data quantities to be written
C     NKPOL    number  of element data quantities to be written
C     LHEAD    T if header and column label are to be written
C     CODE     code used to compute polar
C     VERSION  code version
C     LQUERY   if T, asks permission to overwrite existing file
C
C  Output:
C     ERROR   T if a OPER or WRITE error occurred
C--------------------------------------------------------
      CHARACTER*29 LINE1, LINE2
      CHARACTER*512 LINEL, LINED, LINEF
      CHARACTER*20 CPNAME
      CHARACTER*1 ANS
      LOGICAL LOPEN
C
      ERROR = .FALSE.
C
C---- do we have to open the file?
      LOPEN = FNPOL .NE. ' '
C
      IF(LOPEN) THEN
       OPEN(LU,FILE=FNPOL,STATUS='OLD',ERR=20)
C
       IF(LQUERY) THEN
        WRITE(*,*)
        WRITE(*,*) 'Output file exists.  Overwrite?  Y'
        READ(*,1000) ANS
C
        IF(INDEX('Nn',ANS).EQ.0) GO TO 22
C
        CLOSE(LU)
        WRITE(*,*) 'Polar file not saved'
        RETURN
       ENDIF
C
 20    OPEN(LU,FILE=FNPOL,STATUS='UNKNOWN',ERR=90)
 22    REWIND(LU)
      ENDIF
C
      IF(LHEAD) THEN
       WRITE(LU,*) ' '
       WRITE(LU,8000) CODE, VERSION
       WRITE(LU,*) ' '
       IF(NBL.EQ.1) THEN
        WRITE(LU,9001) NAME
       ELSE
        WRITE(LU,9002) NAME, NBL
       ENDIF
C
       IFFBC = 0
       ISMOM = 0
C
       IF(IFFBC.NE.0 .AND. ISMOM.NE.0) THEN
        IF(IFFBC.EQ.1)  LINE1 = ' Solid wall far field        '
        IF(IFFBC.EQ.2)  LINE1 = ' Vortex + doublet far field  '
        IF(IFFBC.EQ.3)  LINE1 = ' Constant pressure far field '
        IF(IFFBC.EQ.4)  LINE1 = ' Supersonic wave far field   '
        IF(IFFBC.GE.5)  LINE1 = '                             '
        IF(ISMOM.EQ.1)  LINE2 = '   S-momentum conserved      '
        IF(ISMOM.EQ.2)  LINE2 = '   Entropy conserved         '
        IF(ISMOM.EQ.3)  LINE2 = '   Entropy conserved near LE '
        IF(ISMOM.EQ.4)  LINE2 = '   S-mom conserved at shocks '
        IF(ISMOM.GE.5)  LINE2 = '                             '
        WRITE(LU,9006) LINE1, LINE2
 9006   FORMAT(1X,3X,2A29)
       ENDIF
C
       WRITE(LU,*) ' '
C
       LINE1 = ' '
       LINE2 = ' '
       IF(IRETYP.EQ.1) LINE1 = ' Reynolds number fixed       '
       IF(IRETYP.EQ.2) LINE1 = ' Reynolds number ~ 1/sqrt(CL)'
       IF(IRETYP.EQ.3) LINE1 = ' Reynolds number ~ 1/CL      '
       IF(IMATYP.EQ.1) LINE2 = '   Mach number fixed         '
       IF(IMATYP.EQ.2) LINE2 = '   Mach number ~ 1/sqrt(CL)  '
       IF(IMATYP.EQ.3) LINE2 = '   Mach number ~ 1/CL        '
       WRITE(LU,9005) IRETYP, IMATYP, LINE1, LINE2
C
       WRITE(LU,*) ' '
       DO N = 1, NBL
         IS1 = 2*N-1
         IS2 = 2*N
         IF(NBL.EQ.1) THEN
          WRITE(LU,9011) XTRIP(IS1), XTRIP(IS2)
         ELSE
          WRITE(LU,9012) XTRIP(IS1), XTRIP(IS2), N
         ENDIF
       ENDDO
       WRITE(LU,9015) MACH1, REYN1/1.0E6, ACRIT(1), ACRIT(2)
       IF(PTRAT .NE. 0.0) WRITE(LU,9017) PTRAT, ETAP
       WRITE(LU,*) ' '
C
       LINEL = ' '
       LINED = ' '

       KL = 1
       KD = 1

       DO 30 I = 1, NIPOL
         IP = IPOL(I)
         IF(IP.EQ.0) GO TO 30

         KDOT = INDEX(CPOLFORM(IP),'.')
         IF(KDOT.EQ.0) KDOT = LEN(CPOLFORM(IP))
         READ(CPOLFORM(IP)(2:KDOT-1),*,ERR=95) NFORM

         CALL STRIP(CPOLNAME(IP),NNAME) 
         NBLANK = MAX( (NFORM-NNAME+2)/2 , 0 )

         LINEL(KL+1+NBLANK:KL+NNAME+NBLANK) = CPOLNAME(IP)(1:NNAME)
         KL = KL + NFORM

         LINED(KD+2:KD+NFORM) = '--------------------------------'
         KD = KD + NFORM
 30    CONTINUE

       DO 32 J = 1, NJPOL
         JP = JPOL(J)
         IF(JP.EQ.0) GO TO 32

         KDOT = INDEX(CPOLSFORM(JP),'.')
         IF(KDOT.EQ.0) KDOT = LEN(CPOLSFORM(JP))
         READ(CPOLSFORM(JP)(2:KDOT-1),*,ERR=95) NFORM

         CALL STRIP(CPOLSNAME(JP),NNAME)
         NBLANK = MAX( (NFORM-NNAME-2)/2 , 0 )

         DO N = 1, NBL
           LINEL(KL+1:KL+1) = ' '
           LINED(KD+1:KD+1) = ' '
           KL = KL + 1
           KD = KD + 1

           LINEL(KL+1+NBLANK:KL+4+NNAME+NBLANK) = 
     &      'Top_' // CPOLSNAME(JP)(1:NNAME)
           KL = KL + NFORM

           LINED(KD+2:KD+NFORM) = '--------------------------------'
           KD = KD + NFORM

           LINEL(KL+1+NBLANK:KL+4+NNAME+NBLANK) = 
     &      'Bot_' // CPOLSNAME(JP)(1:NNAME)
           KL = KL + NFORM

           LINED(KD+2:KD+NFORM) = '--------------------------------'
           KD = KD + NFORM
         ENDDO
 32    CONTINUE


       IZERO = ICHAR('0')

       DO 35 N = 1, NBL
         LINEL(KL+1:KL+2) = '  '
         LINED(KD+1:KD+2) = '  '
         KL = KL + 2
         KD = KD + 2

         DO 34 K = 1, NKPOL
           KP = KPOL(K)
           IF(KP.EQ.0) GO TO 34

           KDOT = INDEX(CPOLNFORM(KP),'.')
           IF(KDOT.EQ.0) KDOT = LEN(CPOLNFORM(KP))
           READ(CPOLNFORM(KP)(2:KDOT-1),*,ERR=95) NFORM

           CPNAME = CPOLNNAME(KP)
           CALL STRIP(CPNAME,NNAME)
           NBLANK = MAX( (NFORM-NNAME-1)/2 , 0 )
           
           KN = INDEX(CPNAME(1:NNAME),'_n') + 1

           CPNAME(KN:KN) = CHAR(IZERO+N)

           LINEL(KL+1+NBLANK:KL+NNAME+NBLANK) = CPNAME(1:NNAME)
           KL = KL + NFORM

           LINED(KD+2:KD+NFORM) = '--------------------------------'
           KD = KD + NFORM

 34      CONTINUE

 35    CONTINUE

 

C       LINEL = 
C     & '  alpha     CL        CD       CDp       CM    Top_Xtr Bot_Xtr'
CCC     1234567890123456789012345678901234567890123456789012345678901234567890
C       K = 62
C
C
C       LINEL =
C     & ' ------- -------- --------- --------- -------- ------- -------'
CCC       3.453   1.3750   0.00921   0.00512  -0.1450  0.9231  0.5382
CCC       3.453   1.3750   0.00921   0.00213  -0.1450  0.9231  0.5382
C       K = 62

       WRITE(LU,1000) LINEL(1:KL)
       WRITE(LU,1000) LINED(1:KD)
C
      ENDIF
C

      LINEF = '(1X'
      KF = 3
      DO I = 1, NIPOL
        IP = IPOL(I)
        NF = LEN(CPOLFORM(IP))
        LINEF(KF+1:KF+NF+1) = ',' // CPOLFORM(IP)
        KF = KF + NF + 1
      ENDDO
      DO N = 1, NBL
        LINEF(KF+1:KF+3) = ',1X'
        KF = KF + 3
        DO J = 1, NJPOL
          JP = JPOL(J)
          NF = LEN(CPOLSFORM(JP))
          LINEF(KF+1:KF+NF+1) = ',' // CPOLSFORM(JP)
          KF = KF + NF + 1
          LINEF(KF+1:KF+NF+1) = ',' // CPOLSFORM(JP)
          KF = KF + NF + 1
        ENDDO
      ENDDO
      DO N = 1, NBL
        LINEF(KF+1:KF+3) = ',2X'
        KF = KF + 3
        DO K = 1, NKPOL
          KP = KPOL(K)
          NF = LEN(CPOLNFORM(KP))
          LINEF(KF+1:KF+NF+1) = ',' // CPOLNFORM(KP)
          KF = KF + NF + 1
        ENDDO
      ENDDO
      LINEF(KF+1:KF+1) = ')'
      KF = KF + 1

      DO 40 IA = IA1, IA2
        WRITE(LU,LINEF)
     &         (CPOL(IA,IPOL(I)), I=1, NIPOL),
     &        ((CPOLSD(IA,IS,JPOL(J)), J=1, NJPOL), IS=1, 2*NBL),
     &        ((CPOLEL(IA,N ,KPOL(K)), K=1, NKPOL), N=1, NBL)
   40 CONTINUE
C
C
 80   CONTINUE
C---- if file was opened here, then close it
      IF(LOPEN) CLOSE(LU)
      RETURN
C
 90   CONTINUE
      ERROR = .TRUE.
      RETURN
C
 95   CONTINUE
      WRITE(*,*) '? Bad CPOLFORM set up in PINDEX.INC'
      STOP
C
C......................................................................
 1000 FORMAT(A)
 8000 FORMAT(7X,A,9X,'Version', F5.2)
 9001 FORMAT(1X,'Calculated polar for: ', A)
 9002 FORMAT(1X,'Calculated polar for: ', A, I4,' elements')
 9005 FORMAT(1X,I1,I2,2A29)
 9011 FORMAT(1X,
     &'xtrf = ',F7.3,' (top)    ',F9.3,' (bottom)  ')
 9012 FORMAT(1X,
     &'xtrf = ',F7.3,' (top)    ',F9.3,' (bottom)     element', I3)
 9015 FORMAT(1X,
     &'Mach = ',F7.3,5X,'Re = ',F9.3,' e 6',5X,'Ncrit = ',20F7.3)
 9017 FORMAT(1X,
     &'pi_p = ',F7.4,5X,'eta_p = ',F9.4)
 9100 FORMAT(1X,F7.3,F9.4,2F10.5,F9.4,2F8.4  ,  F9.5)
CCC      3.453   1.3750   0.00921     0.500  -0.1450  0.9231  0.5382 -0.00942
CCC      3.453   1.3750   0.00921     0.500  -0.1450  0.9231  0.5382
      END



      SUBROUTINE POLREF(LU,FNREF,ERROR,
     &            NFX,NF,XYREF,LABREF )
      INCLUDE 'PINDEX.INC'
      CHARACTER*(*) FNREF,LABREF
      LOGICAL ERROR
      DIMENSION NF(4)
      DIMENSION XYREF(NFX,2,4)
C--------------------------------------------------------
C     Reads in polar reference data file
C
C  Input:
C     LU      logical unit to use for reading
C     FNREF   name of polar file to be read,
C               if FNREF(1:1).eq.' ', unit LU is assumed 
C               to be already open
C     NFX     polar point array dimension
C
C  Output:
C     ERROR      T if a READ error occurred
C     NF(.)      number of points in each data block
C     XYREF(...) reference polar data
C     LABREF(.)  reference polar label
C--------------------------------------------------------
      LOGICAL LOPEN
      CHARACTER*80 LINE
C
      ERROR = .FALSE.
      LOPEN = FNREF(1:1) .NE. ' '
      IF(LOPEN) OPEN(LU,FILE=FNREF,STATUS='OLD',ERR=900)
C
C---- try to read data label
      READ(LU,1000,END=900) LINE
 1000 FORMAT(A)
C
C---- set data label if present
      IF(LINE(1:1).EQ.'#') THEN
        LABREF = LINE(2:80)
      ELSE
        LABREF = ' '
        REWIND(LU)
      ENDIF
C
      DO 100 K=1, 4
        DO 10 I=1, NFX
          READ(LU,*,END=11,ERR=900) XYREF(I,1,K), XYREF(I,2,K)
          IF(XYREF(I,1,K) .EQ. 999.0) GO TO 11
   10   CONTINUE
   11   NF(K) = I-1
  100 CONTINUE
      IF(LOPEN) CLOSE(LU)
      RETURN
C
  900 CONTINUE
      ERROR = .TRUE.
C
      RETURN
      END ! POLREF
