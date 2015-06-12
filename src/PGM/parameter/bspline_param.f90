subroutine param_custom(kpm, m, n, knot, pos, param)

  implicit none

  !Input
  integer, intent(in) ::  kpm, m, n
  double precision, intent(in) ::  knot(kpm), pos(m)

  !Output
  double precision, intent(out) ::  param(n)

  !Working
  double precision cp(m), pt
  double precision xnew, x0, x, f0, f
  double precision B(kpm-m)
  integer i, j, k, l, i0

  k = kpm - m

  cp(:) = (pos(:) - pos(1)) / (pos(m) - pos(1))

  do l = 1, n
     pt = 1.0 * (l - 1) / (n - 1)
     pt = (pt - pos(1)) / (pos(m) - pos(1))
     if ((0 .le. pt) .and. (pt .le. 1)) then
        x0 = 0
        x = 1

        call basis(k, kpm, x0, knot, B, i0)
        f0 = -pt
        do i = 1, k
           f0 = f0 + B(i) * cp(i0+i)
        end do

        call basis(k, kpm, x, knot, B, i0)
        f = -pt
        do i = 1, k
           f = f + B(i) * cp(i0+i)
        end do

        do j=1,100
           if (abs(f) .lt. 1e-15) then
              exit
           end if
           xnew = x - f * (x-x0) / (f-f0)
           if (xnew .lt. 0) then
              xnew = 0 - xnew
           else if (xnew .gt. 1) then
              xnew = 2 - xnew
           end if
           x0 = x
           x = xnew
           f0 = f

           call basis(k, kpm, x, knot, B, i0)
           f = -pt
           do i = 1, k
              f = f + B(i) * cp(i0+i)
           end do
        end do
        param(l) = x
     else
        param(l) = -1.0
     end if
  end do

end subroutine param_custom
