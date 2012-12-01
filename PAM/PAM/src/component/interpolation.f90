subroutine computeJunction(nQu, nQv, nu1, nu2, nu3, nv1, nv2, nv3, &
     f0, m0, mQT, mQB, mQL, mQR, mA, mB, fQ, Q)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nQu, nQv, nu1, nu2, nu3, nv1, nv2, nv3, f0, m0, mQT, mQB, mQL, mQR, mA, mB, fQ
  !f2py intent(out) Q
  !f2py depend(nv2) mQT, mQB
  !f2py depend(nu2) mQL, mQR
  !f2py depend(nQu, nQv) fQ, Q

  !Input
  integer, intent(in) ::  nQu, nQv, nu1, nu2, nu3, nv1, nv2, nv3
  double precision, intent(in) ::  f0, m0
  double precision, intent(in) ::  mQT(nv2,3), mQB(nv2,3)
  double precision, intent(in) ::  mQL(nu2,3), mQR(nu2,3)
  double precision, intent(in) ::  mA(2,2,3), mB(2,2,3)
  double precision, intent(in) ::  fQ(nQu,nQv,3)

  !Output
  double precision, intent(out) ::  Q(nQu,nQv,3)

  !Working
  double precision vNW(nu1,3), vNE(nu1,3), vSW(nu3,3), vSE(nu3,3)
  double precision hNW(nv1,3), hSW(nv1,3), hNE(nv3,3), hSE(nv3,3)
  double precision A(3), B(3)
  integer u1, u2, v1, v2

  Q(:,:,:) = 0.0

  u1 = nu1
  u2 = nu1 + nu2 - 1
  v1 = nv1
  v2 = nv1 + nv2 - 1

  A = fQ(1,nv1,:)
  B = fQ(2,nv1,:)
  call bezierCurve(nu1, A, f0*(B-A), mA(1,1,:), m0*(mB(1,1,:)-mA(1,1,:)), vNW)
  A = fQ(1,nv1+nv2-1,:)
  B = fQ(2,nv1+nv2-1,:)
  call bezierCurve(nu1, A, f0*(B-A), mA(1,2,:), m0*(mB(1,2,:)-mA(1,2,:)), vNE)

  A = fQ(nQu,nv1,:)
  B = fQ(nQu-1,nv1,:)
  call bezierCurve(nu3, mA(2,1,:), m0*(mB(2,1,:)-mA(2,1,:)), A, f0*(B-A), vSW)
  A = fQ(nQu,nv1+nv2-1,:)
  B = fQ(nQu-1,nv1+nv2-1,:)
  call bezierCurve(nu3, mA(2,2,:), m0*(mB(2,2,:)-mA(2,2,:)), A, f0*(B-A), vSE)

  A = fQ(nu1,1,:)
  B = fQ(nu1,2,:)
  call bezierCurve(nv1, A, f0*(B-A), mA(1,1,:), m0*(mB(1,1,:)-mA(1,1,:)), hNW)
  A = fQ(nu1+nu2-1,1,:)
  B = fQ(nu1+nu2-1,2,:)
  call bezierCurve(nv1, A, f0*(B-A), mA(2,1,:), m0*(mB(2,1,:)-mA(2,1,:)), hSW)

  A = fQ(nu1,nQv,:)
  B = fQ(nu1,nQv-1,:)
  call bezierCurve(nv3, mA(1,2,:), m0*(mB(1,2,:)-mA(1,2,:)), A, f0*(B-A), hNE)
  A = fQ(nu1+nu2-1,nQv,:)
  B = fQ(nu1+nu2-1,nQv-1,:)
  call bezierCurve(nv3, mA(2,2,:), m0*(mB(2,2,:)-mA(2,2,:)), A, f0*(B-A), hSE)

  call coonsPatch(nu1, nv1, fQ(:u1,1,:), vNW, fQ(1,:v1,:), hNW, Q(:u1,:v1,:))
  call coonsPatch(nu1, nv2, vNW, vNE, fQ(1,v1:v2,:), mQT, Q(:u1,v1:v2,:))
  call coonsPatch(nu1, nv3, vNE, fQ(:u1,nQv,:), fQ(1,v2:,:), hNE, Q(:u1,v2:,:))

  call coonsPatch(nu2, nv1, fQ(u1:u2,1,:), mQL, hNW, hSW, Q(u1:u2,:v1,:))
  call coonsPatch(nu2, nv3, mQR, fQ(u1:u2,nQv,:), hNE, hSE, Q(u1:u2,v2:,:))

  call coonsPatch(nu3, nv1, fQ(u2:,1,:), vSW, hSW, fQ(nQu,:v1,:), Q(u2:,:v1,:))
  call coonsPatch(nu3, nv2, vSW, vSE, mQB, fQ(nQu,v1:v2,:), Q(u2:,v1:v2,:))
  call coonsPatch(nu3, nv3, vSE, fQ(u2:,nQv,:), hSE, fQ(nQu,v2:,:), Q(u2:,v2:,:))

end subroutine computeJunction

  

subroutine quad2Dcurve(k, n, P1, P2, n1, n2, P)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) k, n, P1, P2, n1, n2
  !f2py intent(out) P
  !f2py depend(n) P

  !Input
  integer, intent(in) ::  k, n
  double precision, intent(in) ::  P1(3), P2(3), n1(3), n2(3)

  !Output
  double precision, intent(out) ::  P(n,3)

  !Working
  integer i, j, iP
  double precision det, R1, R2, B(3), t, den, cross_ij

  i = k + 1
  j = k + 2
  if (i .gt. 3) then
     i = i - 3
  end if
  if (j .gt. 3) then
     j = j - 3
  end if

  det = cross_ij(i,j,n1,n2)
  R1 = cross_ij(i,j,P1,n1)
  R2 = cross_ij(i,j,P2,n2)

  if (abs(det) .gt. 1e-14) then
     B(i) = (-n2(i)*R1 + n1(i)*R2)/det
     B(j) = (-n2(j)*R1 + n1(j)*R2)/det
     B(k) = 0.5*P1(k) + 0.5*P2(k)
  else
     B(:) = 0.0
  end if

  den = 1.0/(n-1)
  do iP=1,n
     t = (iP-1)*den
     P(iP,:) = (1-t)**2*P1 + 2*t*(1-t)*B + t**2*P2
  end do

end subroutine quad2Dcurve



function cross_ij(i, j, u, v)
  
  implicit none

  integer, intent(in) ::  i, j
  double precision, intent(in) ::  u(3), v(3)
  double precision ::  cross_ij

  cross_ij = u(i)*v(j) - u(j)*v(i)

end function cross_ij



subroutine bezierCurve(n, A, S, B, T, P)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, A, S, B, T
  !f2py intent(out) P
  !f2py depend(n) P

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  A(3), S(3), B(3), T(3)

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



subroutine coonsPatch(nu, nv, Pu0, Pu1, P0v, P1v, P)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nu, nv, Pu0, Pu1, P0v, P1v
  !f2py intent(out) P
  !f2py depend(nu) Pu0, Pu1
  !f2py depend(nv) P0v, P1v
  !f2py depend(nu,nv) P

  !Input
  integer, intent(in) ::  nu, nv
  double precision, intent(in) ::  Pu0(nu,3), Pu1(nu,3), P0v(nv,3), P1v(nv,3)

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
