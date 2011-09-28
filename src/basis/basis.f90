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
subroutine basisold(k, kpm, t, d, B, i0)

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

  if (t .eq. d(1)) then
     i0 = 0
     B(:) = 0
     B(1) = 1
  else if (t .eq. d(m+k)) then
     i0 = m-k
     B(:) = 0
     B(k) = 1
  else
     i0 = -1
     do i=k,m
        if ((d(i) .le. t) .and. (t .lt. d(i+1))) then
           i0 = i-k
        end if
     end do

     B(:) = 0
     B(k) = 1

     do i=2,k
        l = i-1
        j1 = k-l
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
        j2 = k
        n = i0 + j2
        if (d(n+l) .ne. d(n)) then
           B(j2) = (t-d(n))/(d(n+l)-d(n))*B(j2) 
        else
           B(j2) = 0
        end if
     end do
  end if

end subroutine basisold


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
  double precision b1, b2, f1, f2, s1, s2, den, F(k), S(k)
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
     j2 = k
     n = i0 + j2
     if (d(n+l) .ne. d(n)) then
        B(j2) = (t-d(n))/(d(n+l)-d(n))*B(j2) 
     else
        B(j2) = 0
     end if
  end do

end subroutine basis


subroutine basis1(k, kpm, t, d, F, i0)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) k,kpm,t,d
  !f2py intent(out) F,i0
  !f2py depend(kpm) d
  !f2py depend(k) F

  !Input
  integer, intent(in) ::  k, kpm
  double precision, intent(in) ::  t, d(kpm)

  !Output
  double precision, intent(out) ::  F(k)
  integer, intent(out) ::  i0

  !Working
  double precision B(k), b1, b2, f1, f2, den
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

  F(:) = 0.0
  do i=2,k
     l = i-1
     j1 = k-l
     j2 = k
     do j=j1,j2
        n = i0 + j
        if (d(n+l) .ne. d(n)) then
           den = d(n+l)-d(n)
           b1 = (t-d(n))/den*B(j) 
           f1 = (B(j)+(t-d(n))*F(j))/den
        else
           b1 = 0
           f1 = 0
        end if
        if ((j .ne. j2) .and. (d(n+l+1) .ne. d(n+1))) then
           den = d(n+l+1)-d(n+1)
           b2 = (d(n+l+1)-t)/den*B(j+1)
           f2 = ((d(n+l+1)-t)*F(j+1)-B(j+1))/den
        else
           b2 = 0
           f2 = 0
        end if
        B(j) = b1 + b2
        F(j) = f1 + f2
     end do
  end do

end subroutine basis1


subroutine basis2(k, kpm, t, d, S, i0)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) k,kpm,t,d
  !f2py intent(out) S,i0
  !f2py depend(kpm) d
  !f2py depend(k) S

  !Input
  integer, intent(in) ::  k, kpm
  double precision, intent(in) ::  t, d(kpm)

  !Output
  double precision, intent(out) ::  S(k)
  integer, intent(out) ::  i0

  !Working
  double precision B(k), F(k), b1, b2, f1, f2, s1, s2, den
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

  F(:) = 0.0
  S(:) = 0.0
  do i=2,k
     l = i-1
     j1 = k-l
     j2 = k
     do j=j1,j2
        n = i0 + j
        if (d(n+l) .ne. d(n)) then
           den = d(n+l)-d(n)
           b1 = (t-d(n))/den*B(j) 
           f1 = (B(j)+(t-d(n))*F(j))/den
           s1 = (2*F(j)+(t-d(n))*S(j))/den
        else
           b1 = 0
           f1 = 0
           s1 = 0
        end if
        if ((j .ne. j2) .and. (d(n+l+1) .ne. d(n+1))) then
           den = d(n+l+1)-d(n+1)
           b2 = (d(n+l+1)-t)/den*B(j+1)
           f2 = ((d(n+l+1)-t)*F(j+1)-B(j+1))/den
           s2 = ((d(n+l+1)-t)*S(j+1)-2*F(j+1))/den
        else
           b2 = 0
           f2 = 0
           s2 = 0
        end if
        B(j) = b1 + b2
        F(j) = f1 + f2
        if (i .gt. 2) then
           S(j) = s1 + s2
        end if
     end do
  end do

end subroutine basis2
