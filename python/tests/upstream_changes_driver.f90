program upstream_changes_driver
  implicit none
  integer :: mode, n, i, ides, nb
  logical :: ltes
  real :: hk, rt, msq, cf, cf_hk, cf_rt, cf_msq, first, last, sle, factor
  real :: x(5000), y(5000), s(5000), xp(5000), yp(5000), work(5000)
  character(len=32) :: name
  read (*, *) mode
  select case(mode)
  case(1)
    read (*, *) hk, rt, msq, factor
    call test_cffac(factor)
    call cft(hk, rt, msq, cf, cf_hk, cf_rt, cf_msq)
    write (*, '(4ES26.17)') cf, cf_hk, cf_rt, cf_msq
  case(2)
    read (*, *) ides, n, ltes
    if (ides <= 9999) then
      call naca4(ides, s, xp, yp, n, x, y, nb, name, ltes)
    else
      call naca5(ides, s, xp, yp, n, x, y, nb, name, ltes)
    end if
    do i = 1, nb
      write (*, '(2ES26.17)') x(i), y(i)
    end do
  case(3)
    read (*, *) n, first, last
    do i = 1, n
      read (*, *) s(i), x(i)
    end do
    call segspld(x, xp, s, n, first, last)
    write (*, '(ES26.17)') xp(1:n)
  case(4)
    read (*, *) n
    do i = 1, n
      read (*, *) x(i), y(i)
    end do
    call scalc(x, y, s, n)
    call segspl(x, xp, s, n)
    call segspl(y, yp, s, n)
    call lefind(sle, x, xp, y, yp, s, n)
    write (*, '(ES26.17)') sle
  case(5)
    call test_psilin
  case(6)
    call test_mrcl
  case(7)
    call test_tecalc
  end select
end program
