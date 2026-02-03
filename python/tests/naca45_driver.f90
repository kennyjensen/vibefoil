program naca45_driver
  implicit none
  integer :: ides, nside, nb, i
  real, allocatable :: xx(:), yt(:), yc(:), xb(:), yb(:)
  character(len=32) :: name

  read (*, *, end=100) ides, nside

  allocate(xx(nside), yt(nside), yc(nside))
  allocate(xb(2 * nside), yb(2 * nside))

  if (ides .le. 9999) then
     call naca4(ides, xx, yt, yc, nside, xb, yb, nb, name)
  else
     call naca5(ides, xx, yt, yc, nside, xb, yb, nb, name)
  endif

  write (*, '(I0)') nb
  do i = 1, nb
     write (*, '(F20.12,1X,F20.12)') xb(i), yb(i)
  end do

100 continue
end program naca45_driver
