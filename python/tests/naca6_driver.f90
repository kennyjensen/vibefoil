PROGRAM NACA6_DRIVER
  USE NACAauxilary, ONLY: Thickness6, MeanLine6, MeanLine6M, InterpolateCombinedAirfoil
  IMPLICIT NONE

  INTEGER :: nside, family
  REAL :: a, cl, toc
  INTEGER :: i
  REAL :: an, anp, frac
  REAL, ALLOCATABLE, DIMENSION(:) :: x, yt, ymean, ymeanp, yu, yl

  READ(*,*) family, a, cl, toc, nside

  ALLOCATE(x(nside), yt(nside), ymean(nside), ymeanp(nside), yu(nside), yl(nside))

  an = 1.5
  anp = an + 1.0
  DO i = 1, nside
    frac = REAL(i - 1) / REAL(nside - 1)
    IF (i == nside) THEN
      x(i) = 1.0
    ELSE
      x(i) = 1.0 - anp * frac * (1.0 - frac)**an - (1.0 - frac)**anp
    END IF
  END DO

  CALL Thickness6(family, toc, x, yt)
  IF (family < 6) THEN
    CALL MeanLine6(a, cl, x, ymean, ymeanp)
  ELSE
    CALL MeanLine6M(cl, x, ymean, ymeanp)
  END IF
  CALL InterpolateCombinedAirfoil(x, yt, ymean, ymeanp, yu, yl)

  WRITE(*,'(I0)') nside
  DO i = 1, nside
    WRITE(*,'(3ES24.15)') x(i), yu(i), yl(i)
  END DO

  DEALLOCATE(x, yt, ymean, ymeanp, yu, yl)
END PROGRAM NACA6_DRIVER
