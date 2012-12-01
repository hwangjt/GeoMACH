subroutine computeCone1(front, bot, nu, nv, nz, ny, L, dx, &
     shapeR, shapeT, shapeL, shapeB, shape0, Q)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) front, bot, nu, nv, nz, ny, L, dx, shapeR, shapeT, shapeL, shapeB
  !f2py intent(out) Q
  !f2py depend(ny) shapeR, shapeL
  !f2py depend(nz) shapeT, shapeB
  !f2py depend(ny,nz) shape0, Q

  !Input
  logical, intent(in) ::  front, bot
  integer, intent(in) ::  nu, nv, nz, ny
  double precision, intent(in) ::  L, dx
  double precision, intent(in) ::  shapeR(ny,2,3), shapeT(nz,2,3)
  double precision, intent(in) ::  shapeL(ny,2,3), shapeB(nz,2,3)
  double precision, intent(in) ::  shape0(ny,nz)

  !Output
  double precision, intent(out) ::  Q(ny,nz,3)

  !Working
  double precision left(ny,2,3), right(ny,2,3), top(nz,2,3), bottom(nz,2,3)
  double precision hCurves(3,2,nv,3), vCurves(2,3,nu,3)
  double precision pC(3), pL(2,3), pR(2,3), pT(2,3), pB(2,3)
  double precision nL(3), nR(3), nT(3), nB(3)
  double precision e1(3), e2(3), e3(3), px(3)
  double precision dQdw(ny,nz,3)
  integer k

  Q(:,:,:) = 0.0
  dQdw(:,:,:) = 0.0

  e1(:) = 0.0
  e2(:) = 0.0
  e3(:) = 0.0
  e1(1) = 1.0
  e2(2) = 1.0
  e3(3) = 1.0

  pC(1) = L
  pC(2) = 0.0
  pC(3) = 0.0

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

  px = -e1*L/abs(L)*dx
  nL = pL(2,:) - pL(1,:) + px
  nR = pR(2,:) - pR(1,:) + px
  nT = pT(2,:) - pT(1,:) + px
  nB = pB(2,:) - pB(1,:) + px

  call quad2Dcurve(2, nv, pL(1,:), pC, nL, e3, hCurves(2,1,:,:))
  call quad2Dcurve(2, nv, pC, pR(1,:), e3, nR, hCurves(2,2,:,:))
  call quad2Dcurve(3, nu, pT(1,:), pC, nT, e2, vCurves(1,2,:,:))
  call quad2Dcurve(3, nu, pC, pB(1,:), e2, nB, vCurves(2,2,:,:))

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

end subroutine computeCone1



subroutine computeCone2(ax1, ax2, nu, nv, nQ, r, offset, pos, rot, Q0, Q, dQ_drot)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) ax1, ax2, nu, nv, nQ, r, offset, pos, rot, Q0
  !f2py intent(out) Q, dQ_drot
  !f2py depend(nu,nv) Q0, Q
  !f2py depend(nQ) dQ_drot

  !Input
  integer, intent(in) ::  ax1, ax2, nu, nv, nQ
  double precision, intent(in) ::  r(3), offset(3), pos(3), rot(3)
  double precision, intent(in) ::  Q0(nu,nv,3)

  !Output
  double precision, intent(out) ::  Q(nu,nv,3), dQ_drot(nQ,3,3)

  !Working
  integer u, v, k, ax3
  double precision T(3,3), dT_drot(3,3,3)

  ax3 = 6 - ax1 - ax2
  call computeRtnMtx(ax1, ax2, ax3, rot, T, dT_drot)
  
  do u=1,nu
     do v=1,nv
        Q(u,v,:) = matmul(T,Q0(u,v,:)-r) + r + offset + pos
        do k=1,3
           dQ_drot((k-1)*nu*nv+(v-1)*nu+u,:,k) = matmul(dT_drot(:,:,k),Q0(u,v,:)-r)
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
  double precision pi, rz, ry, taU, tbU, taL, tbL

  pi = 2*acos(0.0)
  Q(:,:,:) = 0.0
  do j=1,nj
     rz = radii(j,1)
     ry = radii(j,2)

     taU = atan(ry/rz*fillet(j,1))/pi
     tbU = atan(ry/rz/fillet(j,2))/pi
     taL = atan(ry/rz*fillet(j,3))/pi
     tbL = atan(ry/rz/fillet(j,4))/pi

     call computeRoundedSection(ni, rz, ry, taU, tbU, taL, tbL, &
          t1, t2, shape0(:,j), Q(:,j,3), Q(:,j,2))
  end do
     
end subroutine computeShape



function nMap(t, t1, t2)

  implicit none
  double precision, intent(in) ::  t, t1, t2
  double precision ::  nMap

  nMap = (t - t1)/(t2 - t1)

end function nMap



subroutine computeRoundedSection(n, rz, ry, taU, tbU, taL, tbL, t1, t2, shape0, z, y)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, rz, ry, taU, tbU, taL, tbL, t1, t2, shape0
  !f2py intent(out) z, y
  !f2py depend(n) shape0, z, y

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  rz, ry, taU, tbU, taL, tbL, t1, t2, shape0(n)

  !Output
  double precision, intent(out) ::  z(n), y(n)

  !Working
  integer i
  double precision zU, yU, zL, yL, szU, syU, szL, syL
  double precision pi, t, tt, nz, ny, norm, nMap

  pi = 2*acos(0.0)
  tt = (t2 - t1)/(n - 1)

  zU = ry*tan((0.5-tbU)*pi)
  yU = rz*tan(taU*pi)     
  zL = ry*tan((0.5-tbL)*pi)
  yL = rz*tan(taL*pi)
  
  szU = rz - zU
  syU = ry - yU
  szL = rz - zL
  syL = ry - yL
  
  do i=1,n
     t = t1 + (i-1)*tt
     if (t .lt. 0) then
        t = t + 2.0
     end if
     if (t .le. taU) then
        z(i) = rz
        y(i) = rz*tan(t*pi)
        nz = 1.0
        ny = 0.0
     else if (t .le. tbU) then
        t = nMap(t, taU, tbU)/2.0
        z(i) = zU + szU*cos(t*pi)
        y(i) = yU + syU*sin(t*pi)
        nz = syU*cos(t*pi)
        ny = szU*sin(t*pi)
     else if (t .le. 1-tbU) then
        z(i) = ry*tan((0.5-t)*pi)
        y(i) = ry
        nz = 0.0
        ny = 1.0
     else if (t .le. 1-taU) then
        t = nMap(t, 1-tbU, 1-taU)/2.0 + 0.5
        z(i) = -zU + szU*cos(t*pi)
        y(i) =  yU + syU*sin(t*pi)
        nz = syU*cos(t*pi)
        ny = szU*sin(t*pi)
     else if (t .le. 1+taL) then
        z(i) = -rz
        y(i) = rz*tan((1-t)*pi)
        nz = -1.0
        ny =  0.0
     else if (t .le. 1+tbL) then
        t = nMap(t, 1+taL, 1+tbL)/2.0 + 1.0
        z(i) = -zL + szL*cos(t*pi)
        y(i) = -yL + syL*sin(t*pi)
        nz = syL*cos(t*pi)
        ny = szL*sin(t*pi)
     else if (t .le. 2-tbL) then
        z(i) = -ry*tan((1.5-t)*pi)
        y(i) = -ry
        nz =  0.0
        ny = -1.0
     else if (t .le. 2-taL) then
        t = nMap(t, 2-tbL, 2-taL)/2.0 + 1.5
        z(i) =  zL + szL*cos(t*pi)
        y(i) = -yL + syL*sin(t*pi)
        nz = syL*cos(t*pi)
        ny = szL*sin(t*pi)
     else
        z(i) = rz
        y(i) = rz*tan(t*pi)  
        nz = 1.0
        ny = 0.0
     end if
     norm = (nz**2 + ny**2)**0.5
     z(i) = z(i) + nz/norm*shape0(i)
     y(i) = y(i) + ny/norm*shape0(i)
  end do

end subroutine computeRoundedSection
