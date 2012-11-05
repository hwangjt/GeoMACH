subroutine computeWingRotations(nj, pos, rot0)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nj, pos
  !f2py intent(out) rot0
  !f2py depend(nj) pos, rot0

  !Input
  integer, intent(in) ::  nj
  double precision, intent(in) ::  pos(nj,3)

  !Output
  double precision, intent(out) ::  rot0(nj,3)

  !Working
  integer j
  double precision pi, t(3), t1(3), t2(3), z, one, p, q, v(2)

  pi = 2*acos(0.0)
  z = 0.0
  one = 1.0
  do j=1,nj
     if (j .eq. 1) then
        t = pos(j+1,:) - pos(j,:)
     else if (j .eq. nj) then
        t = pos(j,:) - pos(j-1,:)
     else
        t1 = pos(j,:) - pos(j-1,:)
        t2 = pos(j+1,:) - pos(j,:)
        t = t1/dot_product(t1,t1)**0.5 + t2/dot_product(t2,t2)**0.5
     end if
     v = (/t(3),t(2)/)
     call arc_tan(v, one, one, p)
     v = (/(t(2)**2+t(3)**2)**0.5,t(1)/)
     call arc_tan(v, one, one, q)
     rot0(j,1) = p
     rot0(j,2) = q
     rot0(j,3) = z
  end do

end subroutine computeWingRotations  



subroutine computeWingSections(ni, nj, r, offset, chord, &
     pos, rot, shape0, Q)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) ni, nj, r, offset, chord, pos, rot, shape0
  !f2py intent(out) Q
  !f2py depend(nj) chord, pos, rot
  !f2py depend(ni,nj) shape0, Q

  !Input
  integer, intent(in) ::  ni, nj
  double precision, intent(in) ::  r(3), offset(3)
  double precision, intent(in) ::  chord(nj), pos(nj,3), rot(nj,3)
  double precision, intent(in) ::  shape0(ni,nj,3)

  !Output
  double precision, intent(out) ::  Q(ni,nj,3)

  !Working
  integer i, j
  double precision T(3,3), dT_drot(3,3,3)

  do j=1,nj
     call computeRtnMtx(rot(j,:), T, dT_drot)
     do i=1,ni
        Q(i,j,:) = (matmul(T,shape0(i,j,:)-r) + r)*chord(j) + pos(j,:) + offset
     end do
  end do

end subroutine computeWingSections



subroutine computeRtnMtx(rot, T, dT_drot)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) rot
  !f2py intent(out) T, dT_drot

  !Input
  double precision, intent(in) ::  rot(3)

  !Output
  double precision, intent(out) ::  T(3,3), dT_drot(3,3,3)

  !Working
  double precision p, q, r
  double precision Tp(3,3), Tq(3,3), Tr(3,3)
  double precision dTp(3,3), dTq(3,3), dTr(3,3)

  dTp(:,:) = 0.0
  dTq(:,:) = 0.0
  dTr(:,:) = 0.0

  p = rot(1)
  q = rot(2)
  r = rot(3)

  Tp(:,:) = 0.0
  Tp(1,1) = 1.0
  Tp(2,2) = cos(p)
  Tp(2,3) = sin(p)
  Tp(3,2) = -sin(p)
  Tp(3,3) = cos(p)

  Tq(:,:) = 0.0
  Tq(2,2) = 1.0
  Tq(1,1) = cos(q)
  Tq(1,3) = sin(q)
  Tq(3,1) = -sin(q)
  Tq(3,3) = cos(q)

  Tr(:,:) = 0.0
  Tr(3,3) = 1.0
  Tr(1,1) = cos(r)
  Tr(1,2) = sin(r)
  Tr(2,1) = -sin(r)
  Tr(2,2) = cos(r)

  dTp(2,2) = -sin(p)
  dTp(2,3) = cos(p)
  dTp(3,2) = -cos(p)
  dTp(3,3) = -sin(p)

  dTq(1,1) = -sin(q)
  dTq(1,3) = cos(q)
  dTq(3,1) = -cos(q)
  dTq(3,3) = -sin(q)

  dTr(1,1) = -sin(r)
  dTr(1,2) = cos(r)
  dTr(2,1) = -cos(r)
  dTr(2,2) = -sin(r)

  T = matmul(matmul(Tr,Tq),Tp)
  dT_drot(:,:,1) = matmul(matmul(Tr,Tq),dTp)
  dT_drot(:,:,2) = matmul(matmul(Tr,dTq),Tp)
  dT_drot(:,:,3) = matmul(matmul(dTr,Tq),Tp)

end subroutine computeRtnMtx



subroutine arctan2pi(P, t, dt_dP)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) P
  !f2py intent(out) t, dt_dP

  !Input
  double precision, intent(in) ::  P(2)

  !Output
  double precision, intent(out) ::  t, dt_dP(2)

  !Working
  double precision x, y, pi

  pi = 2*acos(0.0)

  x = P(1)
  y = P(2)

  if (x .eq. 0) then
     if (y .gt. 0) then
        t = pi/2.0
     else if (y .lt. 0) then
        t = 3*pi/2.0
     end if
  else if (y .eq. 0) then
     if (x .gt. 0) then
        t = 0
     else if (x .lt. 0) then
        t = pi
     end if
  else if (x .lt. 0) then
     t = atan(y/x) + pi   
  else if (y .lt. 0) then
     t = atan(y/x) + 2*pi    
  else if (y .gt. 0) then
     t = atan(y/x)
  else
     t = 0
  end if

  dt_dP(1) = -y/(x**2+y**2)
  dt_dP(2) = x/(x**2+y**2)

end subroutine arctan2pi
