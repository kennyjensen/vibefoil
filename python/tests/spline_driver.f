C***********************************************************************
C  Spline driver for SCALC/SEGSPL/SEVAL/DEVAL/D2VAL/SINVRT.
C***********************************************************************
      PROGRAM SPLINEDR
      INTEGER N, NQ, NSINV, I
      REAL X(10000), Y(10000), S(10000), XS(10000)
      REAL SQ, SI, XI
      REAL V1, V2, V3

      READ(*,*,END=900) N, NQ, NSINV
      DO I=1, N
        READ(*,*) X(I)
      END DO
      DO I=1, N
        READ(*,*) Y(I)
      END DO

      CALL SCALC(X,Y,S,N)
      CALL SEGSPL(X,XS,S,N)

      WRITE(*,'(A)', ADVANCE='NO') 'SCALC'
      DO I=1, N
        WRITE(*,'(1X,ES24.16)', ADVANCE='NO') S(I)
      END DO
      WRITE(*,*)

      WRITE(*,'(A)', ADVANCE='NO') 'SEGSPL'
      DO I=1, N
        WRITE(*,'(1X,ES24.16)', ADVANCE='NO') XS(I)
      END DO
      WRITE(*,*)

      DO I=1, NQ
        READ(*,*) SQ
        V1 = SEVAL(SQ, X, XS, S, N)
        V2 = DEVAL(SQ, X, XS, S, N)
        V3 = D2VAL(SQ, X, XS, S, N)
        WRITE(*,'(A,1X,3(1X,ES24.16))') 'EVAL', V1, V2, V3
      END DO

      DO I=1, NSINV
        READ(*,*) SI, XI
        CALL SINVRT(SI, XI, X, XS, S, N)
        WRITE(*,'(A,1X,ES24.16)') 'SINVRT', SI
      END DO

  900 CONTINUE
      END
