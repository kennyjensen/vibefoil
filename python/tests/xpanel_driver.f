C     Fixed-form driver for XFOIL xpanel routines.
      PROGRAM XPANEL_DRIVER
      INCLUDE 'XFOIL.INC'
      INCLUDE 'XBL.INC'
      INTEGER MODE
      INTEGER I, J, NIN, NWIN
      INTEGER GEOLIN, SIGLIN, SHARPF, LIMAGEF
      INTEGER IST_IN
      REAL ALFA_IN, QINF_IN, ANTE_IN, ASTE_IN, DSTE_IN
      REAL QOPI_IN, HOPI_IN, PI_IN, YIMAGE_IN
      REAL XI, YI, NXI, NYI
      REAL PSI, PSI_NI
      REAL WAKLEN_IN, CHORD_IN
C
      READ(*,*,END=900) MODE
C
      IF (MODE .EQ. 1) THEN
C----- NCALC
        READ(*,*) NIN
        DO 10 I=1,NIN
          READ(*,*) X(I), Y(I), S(I)
   10   CONTINUE
        CALL NCALC(X,Y,S,NIN,NX,NY)
        WRITE(*,'(I0)') NIN
        DO 20 I=1,NIN
          WRITE(*,'(2ES24.15)') NX(I), NY(I)
   20   CONTINUE
C
      ELSE IF (MODE .EQ. 2) THEN
C----- APCALC
        READ(*,*) NIN, SHARPF
        READ(*,*) ANTE_IN, ASTE_IN, DSTE_IN
        N = NIN
        SHARP = (SHARPF .NE. 0)
        ANTE = ANTE_IN
        ASTE = ASTE_IN
        DSTE = DSTE_IN
        PI = 4.0*ATAN(1.0)
        DO 30 I=1,NIN
          READ(*,*) X(I), Y(I), NX(I), NY(I)
   30   CONTINUE
        CALL APCALC
        WRITE(*,'(I0)') NIN
        DO 40 I=1,NIN
          WRITE(*,'(ES24.15)') APANEL(I)
   40   CONTINUE
C
      ELSE IF (MODE .EQ. 3) THEN
C----- PSILIN
        READ(*,*) NIN, NWIN, I, GEOLIN, SIGLIN
        READ(*,*) SHARPF, LIMAGEF
        READ(*,*) ALFA_IN, QINF_IN
        READ(*,*) ANTE_IN, ASTE_IN, DSTE_IN
        READ(*,*) QOPI_IN, HOPI_IN, PI_IN
        READ(*,*) YIMAGE_IN
        N = NIN
        NW = NWIN
        SHARP = (SHARPF .NE. 0)
        LIMAGE = (LIMAGEF .NE. 0)
        ALFA = ALFA_IN
        QINF = QINF_IN
        ANTE = ANTE_IN
        ASTE = ASTE_IN
        DSTE = DSTE_IN
        QOPI = QOPI_IN
        HOPI = HOPI_IN
        PI = PI_IN
        YIMAGE = YIMAGE_IN
        DO 50 J=1,N
          READ(*,*) X(J), Y(J), S(J), XP(J), YP(J), NX(J),
     &             NY(J), GAM(J)
   50   CONTINUE
        DO 60 J=1,N+NW
          READ(*,*) SIG(J)
   60   CONTINUE
        DO 70 J=1,N
          READ(*,*) QF0(J), QF1(J), QF2(J), QF3(J)
   70   CONTINUE
        READ(*,*) XI, YI, NXI, NYI
        CALL PSILIN(I,XI,YI,NXI,NYI,PSI,PSI_NI,
     &              GEOLIN.NE.0, SIGLIN.NE.0)
        WRITE(*,'(9ES24.15)') PSI, PSI_NI, QTAN1, QTAN2, QTANM,
     &                        Z_QINF, Z_ALFA, Z_QDOF0, Z_QDOF1
        WRITE(*,'(2ES24.15)') Z_QDOF2, Z_QDOF3
        WRITE(*,'(I0)') N
        DO 80 J=1,N
          WRITE(*,'(4ES24.15)') DZDG(J), DZDN(J), DQDG(J), GAM(J)
   80   CONTINUE
        WRITE(*,'(I0)') N+NW
        DO 90 J=1,N+NW
          WRITE(*,'(2ES24.15)') DZDM(J), DQDM(J)
   90   CONTINUE
C
      ELSE IF (MODE .EQ. 4) THEN
C----- PSWLIN
        READ(*,*) NIN, NWIN, I
        READ(*,*) QOPI_IN, PI_IN
        N = NIN
        NW = NWIN
        QOPI = QOPI_IN
        PI = PI_IN
        DO 100 J=1,N+NW
          READ(*,*) X(J), Y(J)
  100   CONTINUE
        DO 110 J=1,N+NW
          READ(*,*) APANEL(J)
  110   CONTINUE
        DO 120 J=1,N+NW
          READ(*,*) SIG(J)
  120   CONTINUE
        READ(*,*) XI, YI, NXI, NYI
        CALL PSWLIN(I,XI,YI,NXI,NYI,PSI,PSI_NI)
        WRITE(*,'(2ES24.15)') PSI, PSI_NI
        WRITE(*,'(I0)') N+NW
        DO 130 J=1,N+NW
          WRITE(*,'(2ES24.15)') DZDM(J), DQDM(J)
  130   CONTINUE
C
      ELSE IF (MODE .EQ. 5) THEN
C----- GGCALC
        READ(*,*) NIN, SHARPF
        READ(*,*) ALFA_IN, QINF_IN
        READ(*,*) ANTE_IN, ASTE_IN, DSTE_IN
        READ(*,*) QOPI_IN, HOPI_IN, PI_IN
        N = NIN
        SHARP = (SHARPF .NE. 0)
        ALFA = ALFA_IN
        QINF = QINF_IN
        ANTE = ANTE_IN
        ASTE = ASTE_IN
        DSTE = DSTE_IN
        QOPI = QOPI_IN
        HOPI = HOPI_IN
        PI = PI_IN
        DO 140 J=1,N
          READ(*,*) X(J), Y(J), S(J), XP(J), YP(J), NX(J),
     &             NY(J)
  140   CONTINUE
        DO 150 J=1,N
          READ(*,*) SIG(J)
  150   CONTINUE
        DO 160 J=1,N
          READ(*,*) QF0(J), QF1(J), QF2(J), QF3(J)
  160   CONTINUE
        XTE = 0.5*(X(1)+X(N))
        YTE = 0.5*(Y(1)+Y(N))
        CALL GGCALC
        WRITE(*,'(I0)') N
        DO 170 J=1,N
          WRITE(*,'(2ES24.15)') QINVU(J,1), QINVU(J,2)
  170   CONTINUE
C
      ELSE IF (MODE .EQ. 6) THEN
C----- QWCALC
        READ(*,*) NIN, NWIN
        READ(*,*) SHARPF
        READ(*,*) ALFA_IN, QINF_IN
        READ(*,*) ANTE_IN, ASTE_IN, DSTE_IN
        READ(*,*) QOPI_IN, HOPI_IN, PI_IN
        N = NIN
        NW = NWIN
        SHARP = (SHARPF .NE. 0)
        ALFA = ALFA_IN
        QINF = QINF_IN
        ANTE = ANTE_IN
        ASTE = ASTE_IN
        DSTE = DSTE_IN
        QOPI = QOPI_IN
        HOPI = HOPI_IN
        PI = PI_IN
        DO 180 J=1,N+NW
          READ(*,*) X(J), Y(J), S(J), XP(J), YP(J), NX(J),
     &             NY(J)
  180   CONTINUE
        DO 190 J=1,N
          READ(*,*) GAM(J)
  190   CONTINUE
        DO 200 J=1,N+NW
          READ(*,*) SIG(J)
  200   CONTINUE
        DO 210 J=1,N+NW
          READ(*,*) QINVU(J,1), QINVU(J,2)
  210   CONTINUE
        CALL QWCALC
        WRITE(*,'(I0)') N+NW
        DO 220 J=1,N+NW
          WRITE(*,'(2ES24.15)') QINVU(J,1), QINVU(J,2)
  220   CONTINUE
C
      ELSE IF (MODE .EQ. 7) THEN
C----- QDCALC (runs GGCALC)
        READ(*,*) NIN, NWIN, SHARPF
        READ(*,*) ALFA_IN, QINF_IN
        READ(*,*) ANTE_IN, ASTE_IN, DSTE_IN
        READ(*,*) QOPI_IN, HOPI_IN, PI_IN
        N = NIN
        NW = NWIN
        SHARP = (SHARPF .NE. 0)
        ALFA = ALFA_IN
        QINF = QINF_IN
        ANTE = ANTE_IN
        ASTE = ASTE_IN
        DSTE = DSTE_IN
        QOPI = QOPI_IN
        HOPI = HOPI_IN
        PI = PI_IN
        DO 230 J=1,N+NW
          READ(*,*) X(J), Y(J), S(J), XP(J), YP(J), NX(J),
     &             NY(J)
  230   CONTINUE
        DO 240 J=1,N+NW
          READ(*,*) SIG(J)
  240   CONTINUE
        DO 250 J=1,N
          READ(*,*) QF0(J), QF1(J), QF2(J), QF3(J)
  250   CONTINUE
        XTE = 0.5*(X(1)+X(N))
        YTE = 0.5*(Y(1)+Y(N))
        CALL GGCALC
        LADIJ = .FALSE.
        CALL QDCALC
        WRITE(*,'(I0)') N+NW
        DO 260 I=1,N
          WRITE(*,'(6ES24.15)') DIJ(I,1), DIJ(I,2), DIJ(I,N),
     &                          DIJ(I,N+1), DIJ(I,N+NW-1),
     &                          DIJ(I,N+NW)
  260   CONTINUE
C
      ELSE IF (MODE .EQ. 8) THEN
C----- XYWAKE
        READ(*,*) NIN, NWIN, SHARPF
        READ(*,*) ALFA_IN, QINF_IN, WAKLEN_IN, CHORD_IN
        READ(*,*) ANTE_IN, ASTE_IN, DSTE_IN
        READ(*,*) QOPI_IN, HOPI_IN, PI_IN
        N = NIN
        NW = NWIN
        SHARP = (SHARPF .NE. 0)
        ALFA = ALFA_IN
        QINF = QINF_IN
        WAKLEN = WAKLEN_IN
        CHORD = CHORD_IN
        ANTE = ANTE_IN
        ASTE = ASTE_IN
        DSTE = DSTE_IN
        QOPI = QOPI_IN
        HOPI = HOPI_IN
        PI = PI_IN
        DO 270 J=1,N
          READ(*,*) X(J), Y(J), S(J), XP(J), YP(J), NX(J),
     &             NY(J)
  270   CONTINUE
        DO 280 J=1,N
          READ(*,*) GAM(J)
  280   CONTINUE
        DO 290 J=1,N+NW
          READ(*,*) SIG(J)
  290   CONTINUE
        CALL XYWAKE
        WRITE(*,'(I0)') N+NW
        DO 300 J=1,N+NW
          WRITE(*,'(5ES24.15)') X(J), Y(J), S(J), NX(J), NY(J)
  300   CONTINUE
C
      ELSE IF (MODE .EQ. 9) THEN
C----- STFIND
        READ(*,*) NIN
        N = NIN
        DO 310 J=1,N
          READ(*,*) S(J), GAM(J)
  310   CONTINUE
        CALL STFIND
        WRITE(*,'(I0)') IST
        WRITE(*,'(3ES24.15)') SST, SST_GO, SST_GP
C
      ELSE IF (MODE .EQ. 10) THEN
C----- IBLPAN
        READ(*,*) NIN, NWIN, IST_IN
        N = NIN
        NW = NWIN
        IST = IST_IN
        CALL IBLPAN
        WRITE(*,'(I0,1X,I0)') IBLTE(1), IBLTE(2)
        WRITE(*,'(I0,1X,I0)') NBL(1), NBL(2)
        DO 320 I=1,NBL(1)
          WRITE(*,'(I0,1X,I0)') IPAN(I,1), NINT(VTI(I,1))
  320   CONTINUE
        DO 330 I=1,NBL(2)
          WRITE(*,'(I0,1X,I0)') IPAN(I,2), NINT(VTI(I,2))
  330   CONTINUE
C
      ELSE IF (MODE .EQ. 11) THEN
C----- XICALC
        READ(*,*) NIN, NWIN, IST_IN, SHARPF
        READ(*,*) ANTE_IN, ASTE_IN, DSTE_IN
        N = NIN
        NW = NWIN
        IST = IST_IN
        SHARP = (SHARPF .NE. 0)
        ANTE = ANTE_IN
        ASTE = ASTE_IN
        DSTE = DSTE_IN
        DO 340 J=1,N+NW
          READ(*,*) X(J), Y(J), S(J), XP(J), YP(J)
  340   CONTINUE
        READ(*,*) IBLTE(1), IBLTE(2)
        READ(*,*) NBL(1), NBL(2)
        DO 350 I=1,NBL(1)
          READ(*,*) IPAN(I,1)
  350   CONTINUE
        DO 360 I=1,NBL(2)
          READ(*,*) IPAN(I,2)
  360   CONTINUE
        READ(*,*) SST
        CALL XICALC
        WRITE(*,'(I0,1X,I0)') NBL(1), NBL(2)
        DO 370 I=1,NBL(1)
          WRITE(*,'(ES24.15)') XSSI(I,1)
  370   CONTINUE
        DO 380 I=1,NBL(2)
          WRITE(*,'(ES24.15)') XSSI(I,2)
  380   CONTINUE
        WRITE(*,'(I0)') NW
        DO 390 I=1,NW
          WRITE(*,'(ES24.15)') WGAP(I)
  390   CONTINUE
C
      ELSE IF (MODE .EQ. 12) THEN
C----- UICALC
        READ(*,*) NIN, NWIN
        N = NIN
        NW = NWIN
        READ(*,*) NBL(1), NBL(2)
        DO 400 I=1,NBL(1)
          READ(*,*) IPAN(I,1), VTI(I,1)
  400   CONTINUE
        DO 410 I=1,NBL(2)
          READ(*,*) IPAN(I,2), VTI(I,2)
  410   CONTINUE
        DO 420 J=1,N+NW
          READ(*,*) QINV(J), QINV_A(J)
  420   CONTINUE
        CALL UICALC
        DO 430 I=1,NBL(1)
          WRITE(*,'(2ES24.15)') UINV(I,1), UINV_A(I,1)
  430   CONTINUE
        DO 440 I=1,NBL(2)
          WRITE(*,'(2ES24.15)') UINV(I,2), UINV_A(I,2)
  440   CONTINUE
C
      ELSE IF (MODE .EQ. 13) THEN
C----- QVFUE
        READ(*,*) NIN, NWIN
        N = NIN
        NW = NWIN
        READ(*,*) NBL(1), NBL(2)
        DO 450 I=1,NBL(1)
          READ(*,*) IPAN(I,1), VTI(I,1), UEDG(I,1)
  450   CONTINUE
        DO 460 I=1,NBL(2)
          READ(*,*) IPAN(I,2), VTI(I,2), UEDG(I,2)
  460   CONTINUE
        CALL QVFUE
        WRITE(*,'(I0)') N+NW
        DO 470 J=1,N+NW
          WRITE(*,'(ES24.15)') QVIS(J)
  470   CONTINUE
C
      ELSE IF (MODE .EQ. 14) THEN
C----- QISET
        READ(*,*) NIN, NWIN
        READ(*,*) ALFA_IN
        N = NIN
        NW = NWIN
        ALFA = ALFA_IN
        DO 480 J=1,N+NW
          READ(*,*) QINVU(J,1), QINVU(J,2)
  480   CONTINUE
        CALL QISET
        WRITE(*,'(I0)') N+NW
        DO 490 J=1,N+NW
          WRITE(*,'(2ES24.15)') QINV(J), QINV_A(J)
  490   CONTINUE
C
      ELSE IF (MODE .EQ. 15) THEN
C----- GAMQV
        READ(*,*) NIN
        N = NIN
        DO 500 J=1,N
          READ(*,*) QVIS(J), QINV_A(J)
  500   CONTINUE
        CALL GAMQV
        WRITE(*,'(I0)') N
        DO 510 J=1,N
          WRITE(*,'(2ES24.15)') GAM(J), GAM_A(J)
  510   CONTINUE
C
      ELSE IF (MODE .EQ. 16) THEN
C----- STMOVE
        READ(*,*) NIN, NWIN, IST_IN
        READ(*,*) SHARPF
        READ(*,*) ANTE_IN, ASTE_IN, DSTE_IN
        N = NIN
        NW = NWIN
        IST = IST_IN
        SHARP = (SHARPF .NE. 0)
        ANTE = ANTE_IN
        ASTE = ASTE_IN
        DSTE = DSTE_IN
        DO 520 J=1,N
          READ(*,*) X(J), Y(J), S(J), XP(J), YP(J), GAM(J)
  520   CONTINUE
        READ(*,*) IBLTE(1), IBLTE(2)
        READ(*,*) NBL(1), NBL(2)
        DO 530 I=1,NBL(1)
          READ(*,*) IPAN(I,1), VTI(I,1), UEDG(I,1), CTAU(I,1),
     &             THET(I,1), DSTR(I,1), XSSI(I,1)
  530   CONTINUE
        DO 540 I=1,NBL(2)
          READ(*,*) IPAN(I,2), VTI(I,2), UEDG(I,2), CTAU(I,2),
     &             THET(I,2), DSTR(I,2), XSSI(I,2)
  540   CONTINUE
        DO 550 J=1,N+NW
          READ(*,*) QINV(J), QINV_A(J)
  550   CONTINUE
        DO 560 J=1,N+NW
          READ(*,*) QVIS(J)
  560   CONTINUE
        READ(*,*) SST
        CALL STMOVE
        WRITE(*,'(I0)') IST
        WRITE(*,'(3ES24.15)') SST, SST_GO, SST_GP
        WRITE(*,'(I0,1X,I0)') IBLTE(1), IBLTE(2)
        WRITE(*,'(I0,1X,I0)') NBL(1), NBL(2)
        DO 570 I=1,NBL(1)
          WRITE(*,'(6ES24.15)') UEDG(I,1), CTAU(I,1), THET(I,1),
     &                          DSTR(I,1), XSSI(I,1), MASS(I,1)
  570   CONTINUE
        DO 580 I=1,NBL(2)
          WRITE(*,'(6ES24.15)') UEDG(I,2), CTAU(I,2), THET(I,2),
     &                          DSTR(I,2), XSSI(I,2), MASS(I,2)
  580   CONTINUE
      ENDIF
C
  900 CONTINUE
      END
