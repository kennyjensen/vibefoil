C***********************************************************************
C  Xoper driver for SPECAL (QISET) parity.
C***********************************************************************
      PROGRAM XOPERDR
      INCLUDE 'XFOIL.INC'
      INTEGER NIN, NWIN, I
      REAL ALFAI

      READ(*,*,END=900) NIN, NWIN
      READ(*,*) ALFAI

      N = NIN
      NW = NWIN
      ALFA = ALFAI

      DO I=1, N+NW
        READ(*,*) QINVU(I,1), QINVU(I,2)
      END DO

      CALL QISET

      WRITE(*,'(A)', ADVANCE='NO') 'QINV'
      DO I=1, N+NW
        WRITE(*,'(1X,ES24.16)', ADVANCE='NO') QINV(I)
      END DO
      WRITE(*,*)
      WRITE(*,'(A)', ADVANCE='NO') 'QINVA'
      DO I=1, N+NW
        WRITE(*,'(1X,ES24.16)', ADVANCE='NO') QINV_A(I)
      END DO
      WRITE(*,*)

  900 CONTINUE
      END
