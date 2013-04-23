subroutine knotopen(k, kpm, d)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) k,kpm
  !f2py intent(out) d
  !f2py depend(kpm) d

  !Input
  integer, intent(in) ::  k, kpm

  !Output
  double precision, intent(out) ::  d(kpm)

  !Working
  integer i, m
  double precision den

  m = kpm - k
  den = m - k + 1

  do i = 1,k
     d(i) = 0
  end do
  do i = k+1,m
     d(i) = (i-k)/den
  end do
  do i = m+1,m+k
     d(i) = 1
  end do

end subroutine knotopen
