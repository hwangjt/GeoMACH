subroutine computeCone1(front, bot, nu, nv, nz, ny, L, dz, &
     shapeR, shapeT, shapeL, shapeB, Q)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) front, bot, nu, nv, nz, ny, L, dz, shapeR, shapeT, shapeL, shapeB
  !f2py intent(out) Q
  !f2py depend(ny) shapeR, shapeL
  !f2py depend(nz) shapeT, shapeB
  !f2py depend(ny,nz) Q

  !Input
  logical, intent(in) ::  front, bot
  integer, intent(in) ::  nu, nv, nz, ny
  double precision, intent(in) ::  L, dz
  double precision, intent(in) ::  shapeR(ny,2,3), shapeT(nz,2,3)
  double precision, intent(in) ::  shapeL(ny,2,3), shapeB(nz,2,3)

  !Output
  double precision, intent(out) ::  Q(ny,nz,3)

  !Working
  double precision left(ny,2,3), right(ny,2,3), top(nz,2,3), bottom(nz,2,3)
  double precision hCurves(3,2,nv,3), vCurves(2,3,nu,3)
  double precision tempu(nu,3), tempv(nv,3)
  double precision pC(3), pL(2,3), pR(2,3), pT(2,3), pB(2,3)
  double precision nL(3), nR(3), nT(3), nB(3)
  double precision e1(3), e2(3), e3(3), pz(3)
  double precision tempQ(nu,nv,3)

  e1(:) = 0.0
  e2(:) = 0.0
  e3(:) = 0.0
  e1(1) = 1.0
  e2(2) = 1.0
  e3(3) = 1.0

  pC(1) = 0.0
  pC(2) = 0.0
  pC(3) = L

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

  pz = -e3*L/abs(L)*dz
  nL = pL(2,:) - pL(1,:) + pz
  nR = pR(2,:) - pR(1,:) + pz
  nT = pT(2,:) - pT(1,:) + pz
  nB = pB(2,:) - pB(1,:) + pz

  call quad2Dcurve(2, nv, pL(1,:), pC, nL, e1, tempv)
  hCurves(2,1,:,:) = tempv
  call quad2Dcurve(2, nv, pC, pR(1,:), e1, nR, tempv)
  hCurves(2,2,:,:) = tempv
  call quad2Dcurve(1, nu, pT(1,:), pC, nT, e2, tempu)
  vCurves(1,2,:,:) = tempu
  call quad2Dcurve(1, nu, pC, pB(1,:), e2, nB, tempu)
  vCurves(2,2,:,:) = tempu

  call coonsPatch(nu, nv, hCurves(1,1,:,:), hCurves(2,1,:,:), &
       vCurves(1,1,:,:), vCurves(1,2,:,:), tempQ)
  Q(1:nu,1:nv,:) = tempQ
  call coonsPatch(nu, nv, hCurves(1,2,:,:), hCurves(2,2,:,:), &
       vCurves(1,2,:,:), vCurves(1,3,:,:), tempQ)
  Q(1:nu,nz-nv+1:nz,:) = tempQ
  if (bot) then
     call coonsPatch(nu, nv, hCurves(2,1,:,:), hCurves(3,1,:,:), &
          vCurves(2,1,:,:), vCurves(2,2,:,:), tempQ)
     Q(ny-nu+1:ny,1:nv,:) = tempQ
     call coonsPatch(nu, nv, hCurves(2,2,:,:), hCurves(3,2,:,:), &
          vCurves(2,2,:,:), vCurves(2,3,:,:), tempQ)
     Q(ny-nu+1:ny,nz-nv+1:nz,:) = tempQ
  end if

end subroutine computeCone1



subroutine computeCone2(nu, nv, nQ, r, offset, pos, rot, Q0, Q, dQ_drot)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nu, nv, nQ, r, offset, pos, rot, Q0
  !f2py intent(out) Q, dQ_drot
  !f2py depend(nu,nv) Q0, Q
  !f2py depend(nQ) dQ_drot

  !Input
  integer, intent(in) ::  nu, nv, nQ
  double precision, intent(in) ::  r(3), offset(3), pos(3), rot(3)
  double precision, intent(in) ::  Q0(nu,nv,3)

  !Output
  double precision, intent(out) ::  Q(nu,nv,3), dQ_drot(nQ,3)

  !Working
  integer u, v, k
  double precision T(3,3), dT_drot(3,3,3)

  call computeRtnMtx(rot, T, dT_drot)
  
  do u=1,nu
     do v=1,nv
        Q(u,v,:) = matmul(T,Q0(u,v,:)-r) + r + offset + pos
        do k=1,3
           dQ_drot((k-1)*nu*nv+(v-1)*nu+u,k) = r(k) + offset(k) + pos(k) + &
                dot_product(T(k,:),Q0(u,v,:)-r) 
        end do
     end do
  end do

end subroutine computeCone2



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



subroutine computeShape(ni, nj, t1, t2, radii, fillet, shape0, Q)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) ni, nj, t1, t2, radii, fillet, shape0
  !f2py intent(out) Q
  !f2py depend(nj) radii, fillet
  !f2py depend(ni,nj) shape0, Q

  !Input
  integer, intent(in) ::  ni, nj
  double precision, intent(in) ::  t1, t2
  double precision, intent(in) ::  radii(nj,3), fillet(nj,4)
  double precision, intent(in) ::  shape0(ni,nj)

  !Output
  double precision, intent(out) ::  Q(ni,nj,3)

  !Working
  integer j
  double precision rx, ry, taU, tbU, taL, tbL
  double precision x(ni), y(ni)

  Q(:,:,:) = 0.0
  do j=1,nj
     rx = radii(j,1)
     ry = radii(j,2)

     taU = fillet(j,1)/2.0
     tbU = fillet(j,2)/2.0
     taL = fillet(j,3)/2.0
     tbL = fillet(j,4)/2.0

     call computeRoundedSection(ni, rx, ry, taU, tbU, taL, tbL, t1, t2, shape0(:,j), x, y)
     Q(:,j,1) = x
     Q(:,j,2) = y
  end do
     
end subroutine computeShape



subroutine nMap(t, t1, t2, val)
  
  !Input
  double precision, intent(in) ::  t, t1, t2

  !Output
  double precision, intent(out) ::  val

  val = (t - t1)/(t2 - t1)

end subroutine nMap



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
  double precision pi, t, tt, val, nx, ny, norm

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
        call nMap(t, taU, tbU, val)
        t = val/2.0
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
        call nMap(t, 1-tbU, 1-taU, val)
        t = val/2.0 + 0.5
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
        call nMap(t, 1+taL, 1+tbL, val)
        t = val/2.0 + 1.0
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
        call nMap(t, 2-tbL, 2-taL, val)
        t = val/2.0 + 1.5
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
