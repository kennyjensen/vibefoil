C     Fixed-form driver for XFOIL xbl routines.
      PROGRAM XBL_DRIVER
      INCLUDE 'XFOIL.INC'
      INCLUDE 'XBL.INC'
      INTEGER MODE
      INTEGER I, J, IS, NIN, IV
      INTEGER LALFAI
      INTEGER NBL1, NBL2
      INTEGER IBLTE1, IBLTE2
      REAL DSTR_IN, THET_IN, UEDG_IN, MSQ_IN, HKLIM_IN
      REAL CLS_IN, MCLS, RCLS
      REAL MINF1_IN, REINF1_IN
      REAL XLE_IN, YLE_IN, XTE_IN, YTE_IN, SLE_IN, SST_IN
      REAL XSTRIP1, XSTRIP2
      REAL BETA, BETA_MSQ, HERAT, HERAT_MS
C
      READ(*,*,END=900) MODE
C
      IF (MODE .EQ. 1) THEN
C----- DSLIM
        READ(*,*) DSTR_IN, THET_IN, UEDG_IN, MSQ_IN, HKLIM_IN
        DSTR = DSTR_IN
        THET = THET_IN
        UEDG = UEDG_IN
        MSQ = MSQ_IN
        HKLIM = HKLIM_IN
        CALL DSLIM(DSTR,THET,UEDG,MSQ,HKLIM)
        WRITE(*,'(ES24.15)') DSTR
C
      ELSE IF (MODE .EQ. 2) THEN
C----- BLPINI
        CALL BLPINI
        WRITE(*,'(10ES24.15)') SCCON, GACON, GBCON, GCCON, DLCON,
     &                          CTRCON, CTRCEX, DUXCON, CTCON, CFFAC
C
      ELSE IF (MODE .EQ. 3) THEN
C----- MRCL
        READ(*,*) CLS_IN
        READ(*,*) MATYP, RETYP
        READ(*,*) MINF1_IN, REINF1_IN
        MINF1 = MINF1_IN
        REINF1 = REINF1_IN
        CALL MRCL(CLS_IN, MCLS, RCLS)
        WRITE(*,'(4ES24.15)') MINF, REINF, MCLS, RCLS
C
      ELSE IF (MODE .EQ. 4) THEN
C----- COMSET
        READ(*,*) MINF1_IN
        MINF = MINF1_IN
        CALL COMSET
        WRITE(*,'(2ES24.15)') TKLAM, TKL_MSQ
C
      ELSE IF (MODE .EQ. 5) THEN
C----- IBLSYS
        READ(*,*) NBL1, NBL2
        NBL(1) = NBL1
        NBL(2) = NBL2
        CALL IBLSYS
        WRITE(*,'(I0)') NSYS
        DO 50 IS=1,2
          DO 40 I=2,NBL(IS)
            WRITE(*,'(I0,1X,I0,1X,I0)') IS, I, ISYS(I,IS)
   40     CONTINUE
   50   CONTINUE
C
      ELSE IF (MODE .EQ. 6) THEN
C----- XIFSET
        READ(*,*) NIN
        N = NIN
        READ(*,*) IBLTE1, IBLTE2
        IBLTE(1) = IBLTE1
        IBLTE(2) = IBLTE2
        READ(*,*) NBL1, NBL2
        NBL(1) = NBL1
        NBL(2) = NBL2
        DO 152 I=1, NBL1-IBLTE1
          WGAP(I) = 0.0
  152   CONTINUE
        DO 153 I=1, NBL2-IBLTE2
          WGAP(I) = 0.0
  153   CONTINUE
        READ(*,*) XLE_IN, YLE_IN, XTE_IN, YTE_IN
        READ(*,*) SLE_IN, SST_IN
        READ(*,*) XSTRIP1, XSTRIP2
        XLE = XLE_IN
        YLE = YLE_IN
        XTE = XTE_IN
        YTE = YTE_IN
        SLE = SLE_IN
        SST = SST_IN
        XSTRIP(1) = XSTRIP1
        XSTRIP(2) = XSTRIP2
        DO 60 I=1,N
          READ(*,*) X(I), Y(I), S(I)
   60   CONTINUE
        DO 70 I=1,NBL1
          READ(*,*) XSSI(I,1)
   70   CONTINUE
        DO 80 I=1,NBL2
          READ(*,*) XSSI(I,2)
   80   CONTINUE
        READ(*,*) IS
        CALL XIFSET(IS)
        WRITE(*,'(ES24.15)') XIFORC
C
      ELSE IF (MODE .EQ. 7) THEN
C----- UESET
        READ(*,*) NIN
        N = NIN
        READ(*,*) NBL1, NBL2
        NBL(1) = NBL1
        NBL(2) = NBL2
        DO 90 IS=1,2
          DO 100 I=2,NBL(IS)
            READ(*,*) IPAN(I,IS), VTI(I,IS), UINV(I,IS),
     &               MASS(I,IS)
  100     CONTINUE
   90   CONTINUE
        DO 120 I=1,N
          DO 110 J=1,N
            READ(*,*) DIJ(I,J)
  110     CONTINUE
  120   CONTINUE
        CALL UESET
        DO 140 IS=1,2
          DO 130 I=2,NBL(IS)
            WRITE(*,'(I0,1X,I0,1X,ES24.15)') IS, I, UEDG(I,IS)
  130     CONTINUE
  140   CONTINUE
      ELSE IF (MODE .EQ. 8) THEN
C----- MRCHUE (smoke: relies on XIFSET only for small NBL)
        READ(*,*) NIN
        N = NIN
        READ(*,*) MINF, QINF, GAMMA, GAMM1, REINF, HVRAT, ANTE
        GAMBL = GAMMA
        GM1BL = GAMM1
        QINFBL = QINF
        BETA = SQRT(1.0 - MINF**2)
        BETA_MSQ = -0.5/BETA
        TKLAM   = MINF**2 / (1.0 + BETA)**2
        TKL_MSQ =     1.0 / (1.0 + BETA)**2
     &    - 2.0*TKLAM/ (1.0 + BETA) * BETA_MSQ
        TKBL = TKLAM
        TKBL_MS = TKL_MSQ
        RSTBL = (1.0 + 0.5*GM1BL*MINF**2) ** (1.0/GM1BL)
        RSTBL_MS = 0.5*RSTBL/(1.0 + 0.5*GM1BL*MINF**2)
        HSTINV = GM1BL*(MINF/QINFBL)**2 / (1.0 + 0.5*GM1BL*MINF**2)
        HSTINV_MS = GM1BL*(1.0/QINFBL)**2 / (1.0 + 0.5*GM1BL*MINF**2)
     &                - 0.5*GM1BL*HSTINV / (1.0 + 0.5*GM1BL*MINF**2)
        HERAT = 1.0 - 0.5*QINFBL**2*HSTINV
        HERAT_MS = -0.5*QINFBL**2*HSTINV_MS
        REYBL    = REINF * SQRT(HERAT**3) * (1.0+HVRAT)/(HERAT+HVRAT)
        REYBL_RE =         SQRT(HERAT**3) * (1.0+HVRAT)/(HERAT+HVRAT)
        REYBL_MS = REYBL * (1.5/HERAT - 1.0/(HERAT+HVRAT))*HERAT_MS
        READ(*,*) IBLTE1, IBLTE2
        IBLTE(1) = IBLTE1
        IBLTE(2) = IBLTE2
        READ(*,*) NBL1, NBL2
        NBL(1) = NBL1
        NBL(2) = NBL2
        DO 212 I=1, NBL1-IBLTE1
          WGAP(I) = 0.0
  212   CONTINUE
        DO 213 I=1, NBL2-IBLTE2
          WGAP(I) = 0.0
  213   CONTINUE
        READ(*,*) XLE_IN, YLE_IN, XTE_IN, YTE_IN
        READ(*,*) SLE_IN, SST_IN
        READ(*,*) XSTRIP1, XSTRIP2
        XLE = XLE_IN
        YLE = YLE_IN
        XTE = XTE_IN
        YTE = YTE_IN
        SLE = SLE_IN
        SST = SST_IN
        XSTRIP(1) = XSTRIP1
        XSTRIP(2) = XSTRIP2
        DO 160 I=1,N
          READ(*,*) X(I), Y(I), S(I)
 160   CONTINUE
        DO 170 I=1,NBL1
          READ(*,*) XSSI(I,1)
 170   CONTINUE
        DO 180 I=1,NBL2
          READ(*,*) XSSI(I,2)
 180   CONTINUE
        DO 190 I=2,NBL1
          READ(*,*) UEDG(I,1)
  190   CONTINUE
        DO 195 I=2,NBL2
          READ(*,*) UEDG(I,2)
  195   CONTINUE
        CALL BLPINI
        CALL MRCHUE
        WRITE(*,'(I0,1X,I0)') ITRAN(1), ITRAN(2)
        DO 196 I=2,NBL1
          WRITE(*,'(ES24.15)') THET(I,1)
  196   CONTINUE
        DO 197 I=2,NBL2
          WRITE(*,'(ES24.15)') THET(I,2)
  197   CONTINUE
        DO 198 I=2,NBL1
          WRITE(*,'(ES24.15)') DSTR(I,1)
  198   CONTINUE
        DO 199 I=2,NBL2
          WRITE(*,'(ES24.15)') DSTR(I,2)
  199   CONTINUE
        DO 201 I=2,NBL1
          WRITE(*,'(ES24.15)') CTAU(I,1)
  201   CONTINUE
        DO 202 I=2,NBL2
          WRITE(*,'(ES24.15)') CTAU(I,2)
  202   CONTINUE
        DO 203 I=2,NBL1
          WRITE(*,'(ES24.15)') MASS(I,1)
  203   CONTINUE
        DO 204 I=2,NBL2
          WRITE(*,'(ES24.15)') MASS(I,2)
  204   CONTINUE
C
      ELSE IF (MODE .EQ. 9) THEN
C----- MRCHDU (smoke: relies on XIFSET only for small NBL)
        READ(*,*) NIN
        N = NIN
        READ(*,*) MINF, QINF, GAMMA, GAMM1, REINF, HVRAT, ANTE
        GAMBL = GAMMA
        GM1BL = GAMM1
        QINFBL = QINF
        BETA = SQRT(1.0 - MINF**2)
        BETA_MSQ = -0.5/BETA
        TKLAM   = MINF**2 / (1.0 + BETA)**2
        TKL_MSQ =     1.0 / (1.0 + BETA)**2
     &    - 2.0*TKLAM/ (1.0 + BETA) * BETA_MSQ
        TKBL = TKLAM
        TKBL_MS = TKL_MSQ
        RSTBL = (1.0 + 0.5*GM1BL*MINF**2) ** (1.0/GM1BL)
        RSTBL_MS = 0.5*RSTBL/(1.0 + 0.5*GM1BL*MINF**2)
        HSTINV = GM1BL*(MINF/QINFBL)**2 / (1.0 + 0.5*GM1BL*MINF**2)
        HSTINV_MS = GM1BL*(1.0/QINFBL)**2 / (1.0 + 0.5*GM1BL*MINF**2)
     &                - 0.5*GM1BL*HSTINV / (1.0 + 0.5*GM1BL*MINF**2)
        HERAT = 1.0 - 0.5*QINFBL**2*HSTINV
        HERAT_MS = -0.5*QINFBL**2*HSTINV_MS
        REYBL    = REINF * SQRT(HERAT**3) * (1.0+HVRAT)/(HERAT+HVRAT)
        REYBL_RE =         SQRT(HERAT**3) * (1.0+HVRAT)/(HERAT+HVRAT)
        REYBL_MS = REYBL * (1.5/HERAT - 1.0/(HERAT+HVRAT))*HERAT_MS
        READ(*,*) IBLTE1, IBLTE2
        IBLTE(1) = IBLTE1
        IBLTE(2) = IBLTE2
        READ(*,*) NBL1, NBL2
        NBL(1) = NBL1
        NBL(2) = NBL2
        READ(*,*) XLE_IN, YLE_IN, XTE_IN, YTE_IN
        READ(*,*) SLE_IN, SST_IN
        READ(*,*) XSTRIP1, XSTRIP2
        XLE = XLE_IN
        YLE = YLE_IN
        XTE = XTE_IN
        YTE = YTE_IN
        SLE = SLE_IN
        SST = SST_IN
        XSTRIP(1) = XSTRIP1
        XSTRIP(2) = XSTRIP2
        DO 200 I=1,N
          READ(*,*) X(I), Y(I), S(I)
 200   CONTINUE
        DO 210 I=1,NBL1
          READ(*,*) XSSI(I,1)
 210   CONTINUE
        DO 220 I=1,NBL2
          READ(*,*) XSSI(I,2)
 220   CONTINUE
        DO 225 I=2,NBL1
          READ(*,*) THET(I,1), DSTR(I,1), CTAU(I,1),
     &               MASS(I,1), UEDG(I,1)
  225   CONTINUE
        DO 226 I=2,NBL2
          READ(*,*) THET(I,2), DSTR(I,2), CTAU(I,2),
     &               MASS(I,2), UEDG(I,2)
  226   CONTINUE
        CALL BLPINI
        CALL MRCHDU
        WRITE(*,'(I0,1X,I0)') ITRAN(1), ITRAN(2)
        DO 227 I=2,NBL1
          WRITE(*,'(ES24.15)') THET(I,1)
  227   CONTINUE
        DO 228 I=2,NBL2
          WRITE(*,'(ES24.15)') THET(I,2)
  228   CONTINUE
        DO 229 I=2,NBL1
          WRITE(*,'(ES24.15)') DSTR(I,1)
  229   CONTINUE
        DO 230 I=2,NBL2
          WRITE(*,'(ES24.15)') DSTR(I,2)
  230   CONTINUE
        DO 231 I=2,NBL1
          WRITE(*,'(ES24.15)') CTAU(I,1)
  231   CONTINUE
        DO 232 I=2,NBL2
          WRITE(*,'(ES24.15)') CTAU(I,2)
  232   CONTINUE
        DO 233 I=2,NBL1
          WRITE(*,'(ES24.15)') MASS(I,1)
  233   CONTINUE
        DO 234 I=2,NBL2
          WRITE(*,'(ES24.15)') MASS(I,2)
  234   CONTINUE
C
      ELSE IF (MODE .EQ. 10) THEN
C----- UPDATE
        READ(*,*) NIN
        N = NIN
        READ(*,*) NBL1, NBL2
        NBL(1) = NBL1
        NBL(2) = NBL2
        READ(*,*) IBLTE1, IBLTE2
        IBLTE(1) = IBLTE1
        IBLTE(2) = IBLTE2
        READ(*,*) NSYS
        READ(*,*) MATYP, LALFAI
        LALFA = (LALFAI .NE. 0)
        READ(*,*) CL, CLSPEC, ALFA, DTOR
        READ(*,*) MINF, QINF, GAMMA, GAMM1, MINF_CL
        READ(*,*) HVRAT
        DO 240 I=1,N
          READ(*,*) X(I), Y(I)
  240   CONTINUE
        IV = 0
        DO 250 IS=1,2
          DO 245 I=2,NBL(IS)
            IV = IV + 1
            ISYS(I,IS) = IV
            READ(*,*) IPAN(I,IS), VTI(I,IS), UINV(I,IS),
     &               UINV_A(I,IS), DSTR(I,IS), THET(I,IS),
     &               CTAU(I,IS), MASS(I,IS), UEDG(I,IS)
  245     CONTINUE
  250   CONTINUE
        DO 260 IV=1,NSYS
          DO 255 J=1,2
            DO 254 I=1,3
              READ(*,*) VDEL(I,J,IV)
  254       CONTINUE
  255     CONTINUE
  260   CONTINUE
        DO 280 I=1,N
          DO 270 J=1,N
            READ(*,*) DIJ(I,J)
  270     CONTINUE
  280   CONTINUE
        CALL UPDATE
        WRITE(*,'(4ES24.15)') CL, ALFA, RMSBL, RLX
        DO 300 IS=1,2
          DO 290 I=2,NBL(IS)
            WRITE(*,'(I0,1X,I0,5ES24.15)') IS, I, UEDG(I,IS),
     &               THET(I,IS), DSTR(I,IS), CTAU(I,IS), MASS(I,IS)
  290     CONTINUE
  300   CONTINUE
      ENDIF
C
  900 CONTINUE
      END
