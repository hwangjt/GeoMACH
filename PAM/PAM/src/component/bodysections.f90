subroutine computeCone(full, ni, nj, L, y0, y1, y2, ry1, ry2, rz1, rz2, dx, Q)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) full, ni, nj, L, y0, y1, y2, ry1, ry2, rz1, rz2, dx
  !f2py intent(out) Q
  !f2py depend(ni,nj) Q

  !Input
  logical, intent(in) ::  full
  integer, intent(in) ::  ni, nj
  double precision, intent(in) ::  L, y0, y1, y2, ry1, ry2, rz1, rz2, dx
 
  !Output
  double precision, intent(out) ::  Q(ni, nj, 3)

  !Working
  integer i, j
  double precision ii, jj
  double precision C(3,3,3), ir2, pi, z, u, v, one
  double precision Ov(3), lv(3), uO(3), ul(3)
  double precision OO(3), Ol(3), lO(3), ll(3)

  pi = 2*acos(0.0)
  ir2 = 1.0/2**0.5

  z = 0.0
  one = 1.0

  C(1,1,:) = (/ L , y1 - ry1*ir2 , -rz1*ir2 /)
  C(2,1,:) = (/ L , y1 - ry1     , z /)
  C(3,1,:) = (/ L , y1 - ry1*ir2 , rz1*ir2 /)
  C(1,2,:) = (/ L , y1           , -rz1 /)
  C(2,2,:) = (/ z , y0           , z /)
  C(3,2,:) = (/ L , y1           , rz1 /)
  C(1,3,:) = (/ L , y1 + ry1*ir2 ,-rz1*ir2 /)
  C(2,3,:) = (/ L , y1 + ry1     , z /)
  C(3,3,:) = (/ L , y1 + ry1*ir2 , rz1*ir2 /)

  do i=1,ni
     ii = -1.0 + 2.0*(i-1)/(ni-1)
     do j=1,nj
        jj = 1.0*(j-1)/(nj-1)
        if (full .eqv. .True.) then
           jj = -1.0 + 2.0*jj
        end if
        if ((ii .lt. 0) .and. (jj .lt. 0)) then
           OO = C(1,3,:)
           Ol = C(2,3,:)
           lO = C(1,2,:)
           ll = C(2,2,:)
           call computeBoundaryValue(-one, jj, L, y0, y1, y2, ry1, ry2, rz1, rz2, dx, C, Ov)
           call computeBoundaryValue( z, jj, L, y0, y1, y2, ry1, ry2, rz1, rz2, dx, C, lv)
           call computeBoundaryValue( ii,-one, L, y0, y1, y2, ry1, ry2, rz1, rz2, dx, C, uO)
           call computeBoundaryValue( ii, z, L, y0, y1, y2, ry1, ry2, rz1, rz2, dx, C, ul)
           u = ii + 1.0
           v = jj + 1.0
        else if ((ii .lt. 0) .and. (jj .ge. 0)) then
           OO = C(2,3,:)
           Ol = C(3,3,:)
           lO = C(2,2,:)
           ll = C(3,2,:)
           call computeBoundaryValue(-one, jj, L, y0, y1, y2, ry1, ry2, rz1, rz2, dx, C, Ov)
           call computeBoundaryValue( z, jj, L, y0, y1, y2, ry1, ry2, rz1, rz2, dx, C, lv)
           call computeBoundaryValue( ii, z, L, y0, y1, y2, ry1, ry2, rz1, rz2, dx, C, uO)
           call computeBoundaryValue( ii, one, L, y0, y1, y2, ry1, ry2, rz1, rz2, dx, C, ul)
           u = ii + 1.0
           v = jj
        else if ((ii .ge. 0) .and. (jj .lt. 0)) then
           OO = C(1,2,:)
           Ol = C(2,2,:)
           lO = C(1,1,:)
           ll = C(2,1,:)
           call computeBoundaryValue( z, jj, L, y0, y1, y2, ry1, ry2, rz1, rz2, dx, C, Ov)
           call computeBoundaryValue( one, jj, L, y0, y1, y2, ry1, ry2, rz1, rz2, dx, C, lv)
           call computeBoundaryValue( ii,-one, L, y0, y1, y2, ry1, ry2, rz1, rz2, dx, C, uO)
           call computeBoundaryValue( ii, z, L, y0, y1, y2, ry1, ry2, rz1, rz2, dx, C, ul)
           u = ii
           v = jj + 1.0
        else if ((ii .ge. 0) .and. (jj .ge. 0)) then
           OO = C(2,2,:)
           Ol = C(3,2,:)
           lO = C(2,1,:)
           ll = C(3,1,:)
           call computeBoundaryValue( z, jj, L, y0, y1, y2, ry1, ry2, rz1, rz2, dx, C, Ov)
           call computeBoundaryValue( one, jj, L, y0, y1, y2, ry1, ry2, rz1, rz2, dx, C, lv)
           call computeBoundaryValue( ii, z, L, y0, y1, y2, ry1, ry2, rz1, rz2, dx, C, uO)
           call computeBoundaryValue( ii, one, L, y0, y1, y2, ry1, ry2, rz1, rz2, dx, C, ul)
           u = ii
           v = jj
        end if
        Q(i,j,:) = Ov*(1-u) + lv*u + uO*(1-v) + ul*v
        Q(i,j,:) = Q(i,j,:) - OO*(1-u)*(1-v) - Ol*(1-u)*v - lO*u*(1-v) - ll*u*v
     end do
  end do
  
end subroutine computeCone



subroutine computeBoundaryValue(ii, jj, L, y0, y1, y2, ry1, ry2, rz1, rz2, dx, C, val)

  implicit none

  !Input
  double precision, intent(in) ::  ii, jj, L, y0, y1, y2, ry1, ry2, rz1, rz2, dx, C(3,3,3)

  !Output
  double precision, intent(out) ::  val(3)

  !Working
  integer iType, jType
  double precision pi, ir2, t, z, one
  double precision P1(2), P2(2), n1(2), n2(2), B0(2), B(3)

  pi = 2*acos(0.0)
  ir2 = 1.0/2**0.5

  z = 0.0
  one = 1.0

  if (ii .eq. -1.0) then
     iType = 1
  else if (ii .lt. 0.0) then
     iType = 2
  else if (ii .eq. 0.0) then
     iType = 3
  else if (ii .lt. 1.0) then
     iType = 4
  else if (ii .eq. 1.0) then
     iType = 5
  else
     print *, 'Error: iType'
  end if

  if (jj .eq. -1.0) then
     jType = 1
  else if (jj .lt. 0.0) then
     jType = 2
  else if (jj .eq. 0.0) then
     jType = 3
  else if (jj .lt. 1.0) then
     jType = 4
  else if (jj .eq. 1.0) then
     jType = 5
  else
     print *, 'Error: jType'
  end if
     
  if (iType .eq. 1) then
     t = (0.5 - jj*0.25)*pi
     val = (/ L , y1 + ry1*sin(t) , rz1*cos(t) /)
  else if (iType .eq. 5) then
     t = (1.5 + jj*0.25)*pi
     val = (/ L , y1 + ry1*sin(t) , rz1*cos(t) /)
  else if (jType .eq. 1) then
     t = (1.0 + ii*0.25)*pi
     val = (/ L , y1 + ry1*sin(t) , rz1*cos(t) /)
  else if (jType .eq. 5) then
     t = (0.0 - ii*0.25)*pi
     val = (/ L , y1 + ry1*sin(t) , rz1*cos(t) /)
  else if (iType .eq. 2) then
     P1(:) = (/ L , y1 + ry1 /)
     P2(:) = (/ z , y0 /)
     n1(:) = (/ dx, y2 + ry2 - y1 - ry1 /)
     n2(:) = (/ z , one /)
     call computeInterpolantB(P1, P2, n1, n2, B0)
     B = (/ B0(1) , B0(2) , z /)
     t = ii + 1
     val = C(2,3,:)*(1-t)**2 + B*2*t*(1-t) + C(2,2,:)*t**2
     val(3) = C(2,3,3)*(1-t) + C(2,2,3)*t
  else if (iType .eq. 4) then
     P1(:) = (/ z , y0 /)
     P2(:) = (/ L , y1 - ry1 /)
     n1(:) = (/ z , one /)
     n2(:) = (/ dx, y2 - ry2 - y1 + ry1 /)
     call computeInterpolantB(P1, P2, n1, n2, B0)
     B = (/ B0(1) , B0(2) , z /)
     t = ii
     val = C(2,2,:)*(1-t)**2 + B*2*t*(1-t) + C(2,1,:)*t**2
     val(3) = C(2,2,3)*(1-t) + C(2,1,3)*t
  else if (jType .eq. 2) then
     P1(:) = (/ L , -rz1 /)
     P2(:) = (/ z , z /)
     n1(:) = (/ dx, rz1 - rz2 /)
     n2(:) = (/ z , one /)
     call computeInterpolantB(P1, P2, n1, n2, B0)
     B = (/ B0(1) , z , B0(2) /)
     t = jj + 1 
     val = C(1,2,:)*(1-t)**2 + B*2*t*(1-t) + C(2,2,:)*t**2
     val(2) = C(1,2,2)*(1-t) + C(2,2,2)*t
  else if (jType .eq. 4) then
     P1(:) = (/ z , z /)
     P2(:) = (/ L , rz1 /)
     n1(:) = (/ z , one /)
     n2(:) = (/ dx, rz2 - rz1 /)
     call computeInterpolantB(P1, P2, n1, n2, B0)
     B = (/ B0(1) , z , B0(2) /)
     t = jj
     val = C(2,2,:)*(1-t)**2 + B*2*t*(1-t) + C(3,2,:)*t**2
     val(2) = C(2,2,2)*(1-t) + C(3,2,2)*t
  else if ((iType .eq. 3) .and. (jType .eq. 3)) then
     val = C(2,2,:)
  else
     print *, 'Error: boundary value not assigned', iType, jType
  end if     

end subroutine computeBoundaryValue



subroutine computeInterpolantB(P1, P2, n1, n2, B)

  implicit none

  !Input
  double precision, intent(in) ::  P1(2), P2(2), n1(2), n2(2)

  !Output
  double precision, intent(out) ::  B(2)

  !Working
  double precision det, R1, R2

  call crossproduct(n1,n2,det)
  call crossproduct(P1,n1,R1)
  call crossproduct(P2,n2,R2)

  if (abs(det) .gt. 1e-14) then
     B(1) = (-n2(1)*R1 + n1(1)*R2)/det
     B(2) = (-n2(2)*R1 + n1(2)*R2)/det
  else
     B(:) = 0.0
  end if

end subroutine computeInterpolantB



subroutine computeRoundedSection(n, rz, ry, ta1, tb1, ta2, tb2, t1, t2, z, y)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n, rz, ry, ta1, tb1, ta2, tb2, t1, t2
  !f2py intent(out) z, y
  !f2py depend(n) z, y

  !Input
  integer, intent(in) ::  n
  double precision, intent(in) ::  rz, ry, ta1, tb1, ta2, tb2, t1, t2

  !Output
  double precision, intent(out) ::  z(n), y(n)

  !Working
  double precision pi, ta, tb, t, tt, val
  double precision z0, y0, sz, sy, z1, y1, z2, y2
  integer i

  pi = 2*acos(0.0)

  z1 = ry*tan((0.5-tb1/2.0)*pi)
  y1 = rz*tan(ta1*pi/2.0)

  z2 = ry*tan((0.5-tb2/2.0)*pi)
  y2 = rz*tan(ta2*pi/2.0)

  tt = (t2 - t1)/(n - 1)
  do i=1,n
     t = t1 + (i-1)*tt
     if (t .lt. 0) then
        t = t + 2.0
     end if
     if (t .le. 1) then
        ta = ta1/2.0
        tb = tb1/2.0
        z0 = z1
        y0 = y1
     else
        ta = ta2/2.0
        tb = tb2/2.0
        z0 = z2
        y0 = y2
     end if
     sz = rz - z0
     sy = ry - y0
     if ((0 .le. t) .and. (t .le. ta)) then
        z(i) = rz
        y(i) = rz*tan(t*pi)
     else if ((ta .le. t) .and. (t .le. tb)) then
        call nMap(t, ta, tb, val)
        t = val/2.0
        z(i) = z0 + sz*cos(t*pi)
        y(i) = y0 + sy*sin(t*pi)
     else if ((tb .le. t) .and. (t .le. 1-tb)) then
        z(i) = ry*tan((0.5-t)*pi)
        y(i) = ry
     else if ((1-tb .le. t) .and. (t .le. 1-ta)) then
        call nMap(t, 1-tb, 1-ta, val)
        t = val/2.0 + 0.5
        z(i) = -z0 + sz*cos(t*pi)
        y(i) = y0 + sy*sin(t*pi)
     else if ((1-ta .le. t) .and. (t .le. 1+ta)) then
        z(i) = -rz
        y(i) = rz*tan((1-t)*pi)
     else if ((1+ta .le. t) .and. (t .le. 1+tb)) then
        call nMap(t, 1+ta, 1+tb, val)
        t = val/2.0 + 1.0
        z(i) = -z0 + sz*cos(t*pi)
        y(i) = -y0 + sy*sin(t*pi)
     else if ((1+tb .le. t) .and. (t .le. 2-tb)) then
        z(i) = -ry*tan((1.5-t)*pi)
        y(i) = -ry
     else if ((2-tb .le. t) .and. (t .le. 2-ta)) then
        call nMap(t, 2-tb, 2-ta, val)
        t = val/2.0 + 1.5
        z(i) = z0 + sz*cos(t*pi)
        y(i) = -y0 + sy*sin(t*pi)
     else if ((2-ta .le. t) .and. (t .le. 2)) then
        z(i) = rz
        y(i) = rz*tan(t*pi)  
     end if
  end do

end subroutine computeRoundedSection



subroutine nMap(t, t1, t2, val)
  
  !Input
  double precision, intent(in) ::  t, t1, t2

  !Output
  double precision, intent(out) ::  val

  val = (t - t1)/(t2 - t1)

end subroutine nMap
