subroutine computeTipTangents(ax1, ax2, rot, hT, vT, dhT_drot, dvT_drot)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) ax1, ax2, rot
  !f2py intent(out) hT, vT, dhT_drot, dvT_drot
  
  !Input
  integer, intent(in) ::  ax1, ax2
  double precision, intent(in) ::  rot(3)

  !Output
  double precision, intent(out) ::  hT(3), vT(3)
  double precision, intent(out) ::  dhT_drot(3,3), dvT_drot(3,3)

  !Working
  integer k, ax3
  double precision e2(3), e3(3)
  double precision T(3,3), dT_drot(3,3,3)

  ax3 = 6 - ax1 - ax2
  call computeRtnMtx(ax1, ax2, ax3, rot, T, dT_drot)

  e2(:) = 0.0
  e3(:) = 0.0
  e2(2) = 1.0
  e3(3) = 1.0

  hT = matmul(T,e3)
  vT = matmul(T,e2)
  do k=1,3
     dhT_drot(:,k) = matmul(dT_drot(:,:,k),e3)
     dvT_drot(:,k) = matmul(dT_drot(:,:,k),e2)
  end do

end subroutine computeTipTangents



subroutine computeCone(front, bot, nu, nv, nz, ny, f0, m0, pC, hT, vT, &
     shapeR, shapeT, shapeL, shapeB, shape0, Q)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) front, bot, nu, nv, nz, ny, pC, hT, vT, shapeR, shapeT, shapeL, shapeB, shape0
  !f2py intent(out) Q
  !f2py depend(ny) shapeR, shapeL
  !f2py depend(nz) shapeT, shapeB
  !f2py depend(ny,nz) shape0, Q

  !Input
  logical, intent(in) ::  front, bot
  integer, intent(in) ::  nu, nv, nz, ny
  double precision, intent(in) ::  f0, m0, pC(3), hT(3), vT(3)
  double precision, intent(in) ::  shapeR(ny,2,3), shapeT(nz,2,3)
  double precision, intent(in) ::  shapeL(ny,2,3), shapeB(nz,2,3)
  double precision, intent(in) ::  shape0(ny,nz)

  !Output
  double precision, intent(out) ::  Q(ny,nz,3)

  !Working
  double precision left(ny,2,3), right(ny,2,3), top(nz,2,3), bottom(nz,2,3)
  double precision hCurves(3,2,nv,3), vCurves(2,3,nu,3)
  double precision pL(2,3), pR(2,3), pT(2,3), pB(2,3)
  double precision nL(3), nR(3), nT(3), nB(3)
  double precision dQdw(ny,nz,3)
  integer k

  Q(:,:,:) = 0.0
  dQdw(:,:,:) = 0.0

  if (front) then
     left = shapeR(ny:1:-1,:,:)
     top = shapeT
     right = shapeL
     bottom = shapeB(nz:1:-1,:,:)
  else
     left = shapeL
     top = shapeT(nz:1:-1,:,:)
     right = shapeR(ny:1:-1,:,:)
     bottom = shapeB
  end if

  hCurves(1,1,:,:) = top(1:nv,1,:)
  hCurves(1,2,:,:) = top(nz-nv+1:nz,1,:)
  hCurves(3,1,:,:) = bottom(1:nv,1,:)
  hCurves(3,2,:,:) = bottom(nz-nv+1:nz,1,:)
  vCurves(1,1,:,:) = left(1:nu,1,:)
  vCurves(2,1,:,:) = left(ny-nu+1:ny,1,:)
  vCurves(1,3,:,:) = right(1:nu,1,:)
  vCurves(2,3,:,:) = right(ny-nu+1:ny,1,:)

  call midValues(ny, left, pL)
  call midValues(ny, right, pR)
  call midValues(nz, top, pT)
  call midValues(nz, bottom, pB)

  nL = pL(1,:) - pL(2,:)
  nR = pR(2,:) - pR(1,:)
  nT = pT(1,:) - pT(2,:)
  nB = pB(2,:) - pB(1,:)

  call bezierCurve(nv, pL(1,:), nL*f0, pC, -hT*m0, hCurves(2,1,:,:))
  call bezierCurve(nv, pC, hT*m0, pR(1,:), -nR*f0, hCurves(2,2,:,:))
  call bezierCurve(nu, pT(1,:), nT*f0, pC, -vT*m0, vCurves(1,2,:,:))
  call bezierCurve(nu, pC, vT*m0, pB(1,:), -nB*f0, vCurves(2,2,:,:))

  call coonsPatch(nu, nv, vCurves(1,1,:,:), vCurves(1,2,:,:), &
       hCurves(1,1,:,:), hCurves(2,1,:,:), Q(1:nu,1:nv,:), dQdw(1:nu,1:nv,:))
  call coonsPatch(nu, nv, vCurves(1,2,:,:), vCurves(1,3,:,:), &
       hCurves(1,2,:,:), hCurves(2,2,:,:), Q(1:nu,nz-nv+1:nz,:), dQdw(1:nu,nz-nv+1:nz,:))
  if (bot) then
     call coonsPatch(nu, nv, vCurves(2,1,:,:), vCurves(2,2,:,:), &
          hCurves(2,1,:,:), hCurves(3,1,:,:), Q(ny-nu+1:ny,1:nv,:), dQdw(ny-nu+1:ny,1:nv,:))
     call coonsPatch(nu, nv, vCurves(2,2,:,:), vCurves(2,3,:,:), &
          hCurves(2,2,:,:), hCurves(3,2,:,:), Q(ny-nu+1:ny,nz-nv+1:nz,:), dQdw(ny-nu+1:ny,nz-nv+1:nz,:))
  end if

  do k=1,3
     Q(:,:,k) = Q(:,:,k) + shape0(:,:)*dQdw(:,:,k)
  end do

end subroutine computeCone



subroutine midValues(n, P, val)

  implicit none

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  P(n,2,3)

  !Output
  double precision, intent(out) ::  val(2,3)

  !Working
  integer k

  do k=1,2
     val(k,:) = 0.5*P(ceiling((n+1)/2.0),k,:) + 0.5*P(floor((n+1)/2.0),k,:)
  end do

end subroutine midValues



subroutine computeShape(ni, nj, t1, t2, fillet, shape0, Q)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) ni, nj, t1, t2, fillet, shape0
  !f2py intent(out) Q
  !f2py depend(nj) fillet
  !f2py depend(ni,nj) shape0, Q

  !Input
  integer, intent(in) ::  ni, nj
  double precision, intent(in) ::  t1, t2
  double precision, intent(in) ::  fillet(nj,4)
  double precision, intent(in) ::  shape0(ni,nj)

  !Output
  double precision, intent(out) ::  Q(ni,nj,3)

  !Working
  integer j
  double precision pi, taU, tbU, taL, tbL, one

  pi = 2*acos(0.0)
  one = 1.0
  Q(:,:,:) = 0.0
  do j=1,nj
     taU = atan(fillet(j,1))/pi
     tbU = atan(1.0/fillet(j,2))/pi
     taL = atan(fillet(j,3))/pi
     tbL = atan(1.0/fillet(j,4))/pi
     call computeRoundedSection(ni, one, one, taU, tbU, taL, tbL, &
          t1, t2, shape0(:,j), Q(:,j,1), Q(:,j,2))
  end do
     
end subroutine computeShape



subroutine computeRoundedSection(n, rx, ry, taU, tbU, taL, tbL, t1, t2, shape0, x, y)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, rx, ry, taU, tbU, taL, tbL, t1, t2, shape0
  !f2py intent(out) x, y
  !f2py depend(n) shape0, x, y

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  rx, ry, taU, tbU, taL, tbL, t1, t2, shape0(n)

  !Output
  double precision, intent(out) ::  x(n), y(n)

  !Working
  integer i
  double precision xU, yU, xL, yL, sxU, syU, sxL, syL
  double precision pi, t, tt, nx, ny, norm, nMap

  pi = 2*acos(0.0)
  tt = (t2 - t1)/(n - 1)

  xU = ry*tan((0.5-tbU)*pi)
  yU = rx*tan(taU*pi)     
  xL = ry*tan((0.5-tbL)*pi)
  yL = rx*tan(taL*pi)
  
  sxU = rx - xU
  syU = ry - yU
  sxL = rx - xL
  syL = ry - yL
  
  do i=1,n
     t = t1 + (i-1)*tt
     if (t .lt. 0) then
        t = t + 2.0
     end if
     if (t .le. taU) then
        x(i) = rx
        y(i) = rx*tan(t*pi)
        nx = 1.0
        ny = 0.0
     else if (t .le. tbU) then
        t = nMap(t, taU, tbU)/2.0
        x(i) = xU + sxU*cos(t*pi)
        y(i) = yU + syU*sin(t*pi)
        nx = syU*cos(t*pi)
        ny = sxU*sin(t*pi)
     else if (t .le. 1-tbU) then
        x(i) = ry*tan((0.5-t)*pi)
        y(i) = ry
        nx = 0.0
        ny = 1.0
     else if (t .le. 1-taU) then
        t = nMap(t, 1-tbU, 1-taU)/2.0 + 0.5
        x(i) = -xU + sxU*cos(t*pi)
        y(i) =  yU + syU*sin(t*pi)
        nx = syU*cos(t*pi)
        ny = sxU*sin(t*pi)
     else if (t .le. 1+taL) then
        x(i) = -rx
        y(i) = rx*tan((1-t)*pi)
        nx = -1.0
        ny =  0.0
     else if (t .le. 1+tbL) then
        t = nMap(t, 1+taL, 1+tbL)/2.0 + 1.0
        x(i) = -xL + sxL*cos(t*pi)
        y(i) = -yL + syL*sin(t*pi)
        nx = syL*cos(t*pi)
        ny = sxL*sin(t*pi)
     else if (t .le. 2-tbL) then
        x(i) = -ry*tan((1.5-t)*pi)
        y(i) = -ry
        nx =  0.0
        ny = -1.0
     else if (t .le. 2-taL) then
        t = nMap(t, 2-tbL, 2-taL)/2.0 + 1.5
        x(i) =  xL + sxL*cos(t*pi)
        y(i) = -yL + syL*sin(t*pi)
        nx = syL*cos(t*pi)
        ny = sxL*sin(t*pi)
     else
        x(i) = rx
        y(i) = rx*tan(t*pi)  
        nx = 1.0
        ny = 0.0
     end if
     norm = (nx**2 + ny**2)**0.5
     x(i) = x(i) + nx/norm*shape0(i)
     y(i) = y(i) + ny/norm*shape0(i)
  end do

end subroutine computeRoundedSection



function nMap(t, t1, t2)

  implicit none
  double precision, intent(in) ::  t, t1, t2
  double precision ::  nMap

  nMap = (t - t1)/(t2 - t1)

end function nMap
