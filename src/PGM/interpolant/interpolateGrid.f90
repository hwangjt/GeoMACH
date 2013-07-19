subroutine interpolateFrames(i, nu0, nv0, iu, iv, Q, dQdw)

  implicit none

  !Input
  integer, intent(in) ::  i, nu0, nv0, iu(4), iv(4)

  !Output
  double precision, intent(inout) ::  Q(nu0,nv0,3)
  double precision, intent(out) ::  dQdw(nu0,nv0,3)

  !Working
  integer j, nu, nv
  
  dQdw(:,:,:) = 0.0
  do j=1,3
     nv = iv(j+1)-iv(j)
     nu = iu(i+1)-iu(i)
     if ((nu.gt.0) .and. (nv.gt.0)) then
        call coonsPatch(nu+1, nv+1, &
             Q(iu(i):iu(i+1),iv(j),:), &
             Q(iu(i):iu(i+1),iv(j+1),:), &
             Q(iu(i),iv(j):iv(j+1),:), &
             Q(iu(i+1),iv(j):iv(j+1),:), &
             Q(iu(i):iu(i+1),iv(j):iv(j+1),:), &
             dQdw(iu(i):iu(i+1),iv(j):iv(j+1),:))
     end if
  end do  

end subroutine interpolateFrames



subroutine bezierCurve(n, P0, D0, P1, D1, P)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, P0, D0, P1, D1
  !f2py intent(out) P
  !f2py depend(n) P

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  P0(3), D0(3), P1(3), D1(3)

  !Output
  double precision, intent(out) ::  P(n,3)

  !Working
  double precision den, u
  double precision A(3), B(3), C(3), D(3)
  logical given0, given1
  integer i

  den = 1.0/(n-1)

  if (dot_product(D0,D0)**0.5 .lt. 1e-10) then
     given0 = .False.
  else
     given0 = .True.
  end if

  if (dot_product(D1,D1)**0.5 .lt. 1e-10) then
     given1 = .False.
  else
     given1 = .True.
  end if

  if (given0 .and. given1) then
     A = P0
     B = D0/3.0 + P0
     C = D1/3.0 + P1
     D = P1
  elseif (given0 .and. (.not. given1)) then
     A = P0
     B = D0/3.0 + P0
     C = D0/3.0 + 2.0/3.0*P0 + P1/3.0
     D = P1
  elseif ((.not. given0) .and. given1) then
     A = P0
     B = D1/3.0 + 2.0/3.0*P1 + P0/3.0
     C = D1/3.0 + P1
     D = P1
  else
     A = P0
     B = 2.0/3.0*P0 + P1/3.0
     C = 2.0/3.0*P1 + P0/3.0
     D = P1
  end if

  do i=1,n
     u = (i-1)*den
     P(i,:) = (1-u)**3*A + 3*u*(1-u)**2*B + 3*u**2*(1-u)*C + u**3*D
  end do

end subroutine bezierCurve
