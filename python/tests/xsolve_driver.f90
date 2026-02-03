      PROGRAM XSOLVDR
      INCLUDE 'XFOIL.INC'
      INTEGER MODE
      INTEGER NSIZ, NN, NRHS
      INTEGER N, NSYS
      INTEGER I, J, K, L
      INTEGER IBL1, IBL2, IVTE1, IVZ
      REAL VACC
      REAL S1, S2
      INTEGER, PARAMETER :: MAXN = 64
      REAL Z(MAXN,MAXN), R(MAXN,MAXN)
      REAL A(MAXN,MAXN), B(MAXN)
      INTEGER INDX(MAXN)
      INTEGER NVA, NVM, NVDEL

      READ(*,*,END=900) MODE

      IF (MODE.EQ.1) THEN
        READ(*,*) NSIZ, NN, NRHS
        DO I=1, NN
          DO J=1, NN
            READ(*,*) Z(I,J)
          END DO
        END DO
        DO I=1, NN
          DO J=1, NRHS
            READ(*,*) R(I,J)
          END DO
        END DO
        CALL GAUSS(MAXN, NN, Z, R, NRHS)
        WRITE(*,*) NN, NRHS
        DO I=1, NN
          WRITE(*,'(100(1X,ES24.16))') (Z(I,J), J=1, NN)
        END DO
        DO I=1, NN
          WRITE(*,'(100(1X,ES24.16))') (R(I,J), J=1, NRHS)
        END DO
        GO TO 900
      END IF

      IF (MODE.EQ.2) THEN
        READ(*,*) NSIZ, N
        DO I=1, N
          DO J=1, N
            READ(*,*) A(I,J)
          END DO
        END DO
        DO I=1, N
          READ(*,*) B(I)
        END DO
        CALL LUDCMP(MAXN, N, A, INDX)
        CALL BAKSUB(MAXN, N, A, INDX, B)
        WRITE(*,*) N
        DO I=1, N
          WRITE(*,'(100(1X,ES24.16))') (A(I,J), J=1, N)
        END DO
        WRITE(*,'(100(1X,I8))') (INDX(I), I=1, N)
        WRITE(*,'(100(1X,ES24.16))') (B(I), I=1, N)
        GO TO 900
      END IF

      IF (MODE.EQ.3) THEN
        READ(*,*) NSYS
        READ(*,*) IBL1, IBL2, IVTE1, IVZ
        READ(*,*) VACC
        READ(*,*) S1, S2

        NSYS = NSYS
        N = 2
        VACCEL = VACC
        S(1) = S1
        S(N) = S2

        IBLTE(1) = IBL1
        IBLTE(2) = IBL2
        DO I=1, NSYS
          ISYS(I,1) = 0
          ISYS(I,2) = 0
        END DO
        IF (IBL1.GE.1 .AND. IBL1.LE.NSYS) THEN
          ISYS(IBL1,1) = IVTE1
        END IF
        IF (IBL2+1.GE.1 .AND. IBL2+1.LE.NSYS) THEN
          ISYS(IBL2+1,2) = IVZ
        END IF

        DO K=1, 3
          DO J=1, 2
            DO I=1, NSYS
              READ(*,*) VA(K,J,I)
            END DO
          END DO
        END DO
        DO K=1, 3
          DO J=1, 2
            DO I=1, NSYS
              READ(*,*) VB(K,J,I)
            END DO
          END DO
        END DO
        DO K=1, 3
          DO L=1, NSYS
            DO I=1, NSYS
              READ(*,*) VM(K,L,I)
            END DO
          END DO
        END DO
        DO K=1, 3
          DO J=1, 2
            DO I=1, NSYS
              READ(*,*) VDEL(K,J,I)
            END DO
          END DO
        END DO
        DO K=1, 3
          DO J=1, 2
            READ(*,*) VZ(K,J)
          END DO
        END DO

        CALL BLSOLV

        NVA = 3 * 2 * NSYS
        NVM = 3 * NSYS * NSYS
        NVDEL = 3 * 2 * NSYS

        WRITE(*,*) NSYS
        WRITE(*,*) NVA
        DO K=1, 3
          DO J=1, 2
            DO I=1, NSYS
              WRITE(*,'(ES24.16)') VA(K,J,I)
            END DO
          END DO
        END DO
        WRITE(*,*) NVM
        DO K=1, 3
          DO L=1, NSYS
            DO I=1, NSYS
              WRITE(*,'(ES24.16)') VM(K,L,I)
            END DO
          END DO
        END DO
        WRITE(*,*) NVDEL
        DO K=1, 3
          DO J=1, 2
            DO I=1, NSYS
              WRITE(*,'(ES24.16)') VDEL(K,J,I)
            END DO
          END DO
        END DO
        GO TO 900
      END IF

  900 CONTINUE
      END
