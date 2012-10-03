subroutine bezierCurve(n, A, B, S, T, P)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, A, B, S, T
  !f2py intent(out) P
  !f2py depend(n) P

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  A(3), B(3), S(3), T(3)

  !Output
  double precision, intent(out) ::  P(n,3)

  !Working
  double precision den, u
  double precision P0(3), P1(3), P2(3), P3(3)
  logical Sgiven, Tgiven
  integer i

  den = 1.0/(n-1)

  if (dot_product(S,S) .lt. 1e-14) then
     Sgiven = .false.
  else
     Sgiven = .true.
  end if

  if (dot_product(T,T) .lt. 1e-14) then
     Tgiven = .false.
  else
     Tgiven = .true.
  end if

  if ((.not. Sgiven) .and. (.not. Tgiven)) then
     P0 = A
     P1 = B
     P2 = 0
     P3 = 0
  elseif ((Sgiven) .and. (.not. Tgiven)) then
     P0 = A
     P1 = S/2.0 + A
     P2 = B
     P3 = 0
  elseif ((.not. Sgiven) .and. (Tgiven)) then
     P0 = A
     P1 = T/2.0 + B
     P2 = B
     P3 = 0
  else
     P0 = A
     P1 = S/3.0 + A
     P2 = T/3.0 + B
     P3 = B
  end if

  do i=1,n
     u = (i-1)*den
     if ((.not. Sgiven) .and. (.not. Tgiven)) then
        P(i,:) = (1-u)*P0 + u*P1
     elseif (Sgiven .and. Tgiven) then
        P(i,:) = (1-u)**3*P0 + 3*u*(1-u)**2*P1 + 3*u**2*(1-u)*P2 + u**3*P3
     else
        P(i,:) = (1-u)**2*P0 + 2*u*(1-u)*P1 + u**2*P2
     end if
  end do

end subroutine bezierCurve



subroutine coonsPatch(nu, nv, P0v, P1v, Pu0, Pu1, P)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nu, nv, P0v, P1v, Pu0, Pu1
  !f2py intent(out) P
  !f2py depend(nu) Pu0, Pu1
  !f2py depend(nv) P0v, P1v
  !f2py depend(nu,nv) P

  !Input
  integer, intent(in) ::  nu, nv
  double precision, intent(in) ::  P0v(nv,3), P1v(nv,3), Pu0(nu,3), Pu1(nu,3)

  !Output
  double precision, intent(out) ::  P(nu,nv,3)

  !Working
  double precision P00(3), P01(3), P10(3), P11(3)
  double precision denu, denv
  double precision u, v
  integer i, j

  P00 = P0v(1,:)
  P10 = P1v(1,:)
  P01 = P0v(nv,:)
  P11 = P1v(nv,:)

  denu = 1.0/(nu-1)
  denv = 1.0/(nv-1)

  do i=1,nu
     u = (i-1)*denu
     do j=1,nv
        v = (j-1)*denv
        P(i,j,:) = (1-u)*P0v(j,:) + u*P1v(j,:) + (1-v)*Pu0(i,:) + v*Pu1(i,:)
        P(i,j,:) = P(i,j,:) - (1-u)*(1-v)*P00 - u*(1-v)*P10 - (1-u)*v*P01 - u*v*P11
     end do
  end do

end subroutine coonsPatch
