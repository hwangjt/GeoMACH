!        k  i0       \t/   m      m+k
!  D D D D...D D D D D D...D D D D D
!              0 0 0 1
!              0 0 B B
!              0 B B B
!               j1  j2
!                 .
!                 .
!                 .
!              B B B B
!

subroutine basis(k, kpm, t, d, B, i0)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) k,kpm,t,d
  !f2py intent(out) B,i0
  !f2py depend(kpm) d
  !f2py depend(k) B

  !Input
  integer, intent(in) ::  k, kpm
  double precision, intent(in) ::  t, d(kpm)

  !Output
  double precision, intent(out) ::  B(k)
  integer, intent(out) ::  i0

  !Working
  integer i, j, j1, j2, l, m, n

  m = kpm - k

  i0 = -1
  do i=k,m
     if ((d(i) .le. t) .and. (t .lt. d(i+1))) then
        i0 = i-k
     end if
  end do

  B(:) = 0
  B(k) = 1

  if (t .eq. d(m+k)) then
     i0 = m-k
  end if

  do i=2,k
     l = i-1
     j1 = k-l
     j2 = k
     n = i0 + j1
     if (d(n+l+1) .ne. d(n+1)) then
        B(j1) = (d(n+l+1)-t)/(d(n+l+1)-d(n+1))*B(j1+1)
     else
        B(j1) = 0
     end if
     do j=j1+1,j2-1
        n = i0 + j
        if (d(n+l) .ne. d(n)) then
           B(j) = (t-d(n))/(d(n+l)-d(n))*B(j) 
        else
           B(j) = 0
        end if
        if (d(n+l+1) .ne. d(n+1)) then
           B(j) = B(j) + (d(n+l+1)-t)/(d(n+l+1)-d(n+1))*B(j+1)
        end if
     end do
     n = i0 + j2
     if (d(n+l) .ne. d(n)) then
        B(j2) = (t-d(n))/(d(n+l)-d(n))*B(j2) 
     else
        B(j2) = 0
     end if
  end do

end subroutine basis
