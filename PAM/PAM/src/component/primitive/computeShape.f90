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
