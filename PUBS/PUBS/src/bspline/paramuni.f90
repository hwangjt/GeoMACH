subroutine paramuni(kpm, m, n, d, t)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) kpm,m,n,d
  !f2py intent(out) t
  !f2py depend(kpm) d
  !f2py depend(n) t

  !Input
  integer, intent(in) ::  kpm, m, n
  double precision, intent(in) ::  d(kpm)

  !Output
  double precision, intent(out) ::  t(n)

  !Working
  double precision P(n), C(m), den
  double precision xnew,x0,x,f0,f
  double precision B(kpm-m)
  integer i, j, k, l, i0

  k = kpm - m

  den = m-1
  do l=1,m
     C(l) = (l-1)/den
  end do

  den = n-1
  do l=1,n
     P(l) = (l-1)/den
     x0 = 0
     x = 1

     call basis(k, kpm, x0, d, B, i0)
     f0 = -P(l)
     do i=1,k
        f0 = f0 + B(i)*C(i0+i)
     end do

     call basis(k, kpm, x, d, B, i0)
     f = -P(l)
     do i=1,k
        f = f + B(i)*C(i0+i)
     end do

     do j=1,100
        !print *,j,abs(f)
        if (abs(f) .lt. 1e-15) then
           exit
        end if
        xnew = x - f*(x-x0)/(f-f0)
        if (xnew .lt. 0) then
           xnew = 0 - xnew
        else if (xnew .gt. 1) then
           xnew = 2 - xnew
        end if
        x0 = x
        x = xnew
        f0 = f

        call basis(k, kpm, x, d, B, i0)
        f = -P(l)
        do i=1,k
           f = f + B(i)*C(i0+i)
        end do
     end do
     t(l) = x
  end do

end subroutine paramuni
