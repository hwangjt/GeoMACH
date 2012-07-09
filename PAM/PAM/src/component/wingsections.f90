subroutine computeWingRotations(nj, posx, posy, posz, rot0)

  implicit none

  !Input
  integer, intent(in) ::  nj
  double precision, intent(in) ::  posx(nj), posy(nj), posz(nj)

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
        t(1) = posx(j+1) - posx(j)
        t(2) = posy(j+1) - posy(j)
        t(3) = posz(j+1) - posz(j)
     else if (j .eq. nj) then
        t(1) = posx(j) - posx(j-1)
        t(2) = posy(j) - posy(j-1)
        t(3) = posz(j) - posz(j-1)
     else
        t1(1) = posx(j) - posx(j-1)
        t1(2) = posy(j) - posy(j-1)
        t1(3) = posz(j) - posz(j-1)
        t2(1) = posx(j+1) - posx(j)
        t2(2) = posy(j+1) - posy(j)
        t2(3) = posz(j+1) - posz(j)
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



subroutine computeWingSections(ni, nj, rx, ry, posx, posy, posz, &
     rotx, roty, rotz, prpx, prpy, chord, shape0, Q)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) ni, nj, rx, ry, posx, posy, posz, rotx, roty, rotz, prpx, prpy, chord, shape0
  !f2py intent(out) Q
  !f2py depend(nj) posx, posy, posz, rotx, roty, rotz, prpx, prpy, chord
  !f2py depend(ni,nj) shape0
  !f2py depend(ni,nj) Q

  !Input
  integer, intent(in) ::  ni, nj
  double precision, intent(in) ::  rx, ry
  double precision, intent(in) ::  posx(nj), posy(nj), posz(nj)
  double precision, intent(in) ::  rotx(nj), roty(nj), rotz(nj)
  double precision, intent(in) ::  prpx(nj), prpy(nj)
  double precision, intent(in) ::  chord(nj)
  double precision, intent(in) ::  shape0(ni,nj,3)

  !Output
  double precision, intent(out) ::  Q(ni,nj,3)

  !Working
  integer i, j
  double precision T(3,3), rotj(3), shapej(ni,3), rot0(nj,3)

  call computeWingRotations(nj, posx, posy, posz, rot0)

  do j=1,nj
     rotj(1) = rotx(j) + rot0(j,1)*prpx(j)
     rotj(2) = roty(j) + rot0(j,2)*prpy(j)
     rotj(3) = rotz(j)
     call computeRtnMtx(rotj, T)
     shapej(:,1) = shape0(:,j,1) - rx
     shapej(:,2) = shape0(:,j,2) - ry
     shapej(:,3) = shape0(:,j,3)
     do i=1,ni        
        Q(i,j,:) = matmul(T,shapej(i,:))
        Q(i,j,1) = Q(i,j,1) + rx
        Q(i,j,2) = Q(i,j,2) + ry
        Q(i,j,:) = Q(i,j,:)*chord(j)
     end do
     Q(:,j,1) = Q(:,j,1) + posx(j)
     Q(:,j,2) = Q(:,j,2) + posy(j)
     Q(:,j,3) = Q(:,j,3) + posz(j)
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
