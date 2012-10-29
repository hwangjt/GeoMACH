subroutine computeWingRotations(nj, pos, rot0)

  implicit none

  !Input
  integer, intent(in) ::  nj
  double precision, intent(in) ::  pos(nj,3)

  !Output
  double precision, intent(out) ::  rot0(nj,3)

  !Working
  integer j
  double precision t(3), t1(3), t2(3), z, one, p, q, v(2), pi

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
     call arc_tan(dot_product(v,v)**0.5, one, one, q)
     rot0(j,1) = p*180.0/pi
     rot0(j,2) = q*180.0/pi
     rot0(j,3) = z
  end do

end subroutine computeWingRotations  



subroutine computeWingSections(ni, nj, r, offset, chord, &
     pos, rot, nor, shape0, Q)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) ni, nj, r, offset, chord, pos, rot, nor, shape0
  !f2py intent(out) Q
  !f2py depend(nj) chord, pos, rot, nor
  !f2py depend(ni,nj) shape0, Q

  !Input
  integer, intent(in) ::  ni, nj
  double precision, intent(in) ::  r(3), offset(3)
  double precision, intent(in) ::  chord(nj), pos(nj,3), rot(nj,3), nor(nj,3)
  double precision, intent(in) ::  shape0(ni,nj,3)

  !Output
  double precision, intent(out) ::  Q(ni,nj,3)

  !Working
  integer i, j
  double precision T(3,3), rot0(nj,3)

  call computeWingRotations(nj, pos, rot0)

  do j=1,nj
     call computeRtnMtx(rot(j,:) + rot0(j,:)*nor(j,:), T)
     do i=1,ni
        Q(i,j,:) = (matmul(T,shape0(i,j,:)-r) + r)*chord(j) + pos(j,:) + offset
     end do
  end do

end subroutine computeWingSections



subroutine computeRtnMtx(rot, T)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) rot
  !f2py intent(out) T

  !Input
  double precision, intent(in) ::  rot(3)

  !Output
  double precision, intent(out) ::  T(3,3)

  !Working
  double precision p, q, r, T0(3,3), pi
  integer k

  pi = 2*acos(0.0)

  p = rot(1)*pi/180.0
  q = rot(2)*pi/180.0
  r = rot(3)*pi/180.0

  T(:,:) = 0.0
  do k=1,3
     T(k,k) = 1.0
  end do

  T0(:,:) = 0.0
  T0(1,1) = 1.0
  T0(2,2) = cos(p)
  T0(2,3) = sin(p)
  T0(3,2) = -sin(p)
  T0(3,3) = cos(p)
  T = matmul(T0, T)

  T0(:,:) = 0.0
  T0(2,2) = 1.0
  T0(1,1) = cos(q)
  T0(1,3) = sin(q)
  T0(3,1) = -sin(q)
  T0(3,3) = cos(q)
  T = matmul(T0, T)

  T0(:,:) = 0.0
  T0(3,3) = 1.0
  T0(1,1) = cos(r)
  T0(1,2) = sin(r)
  T0(2,1) = -sin(r)
  T0(2,2) = cos(r)
  T = matmul(T0, T)

end subroutine computeRtnMtx
