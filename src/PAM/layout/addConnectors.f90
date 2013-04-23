subroutine addConnectors(nvert, nedge, nvert0, nedge0, verts0, edges0, &
     quadrants, verts, edges)
  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvert, nedge, nvert0, nedge0, verts0, edges0, quadrants
  !f2py intent(out) verts, edges
  !f2py depend(nvert0) verts0
  !f2py depend(nedge0) edges0
  !f2py depend(nvert0) quadrants
  !f2py depend(nvert) verts
  !f2py depend(nedge) edges

  !Input
  integer, intent(in) ::  nvert, nedge, nvert0, nedge0
  double precision, intent(in) ::  verts0(nvert0,2), edges0(nedge0,5)
  logical, intent(in) ::  quadrants(nvert0,4)

  !Output
  double precision, intent(out) ::  verts(nvert,2), edges(nedge,5)

  !Working
  integer v, e, v0
  double precision v1, v2

  do v=1,nvert0
     verts(v,:) = verts0(v,:)
  end do

  do e=1,nedge0
     edges(e,:) = edges0(e,:)
  end do

  v = nvert0 + 1
  e = nedge0 + 1
  do v0=1,nvert0
     v1 = verts0(v0,1)
     v2 = verts0(v0,2)
     if (.not. ((v1.eq.0) .or. (v1.eq.1) .or. (v2.eq.0) .or. (v2.eq.1))) then
        if (.not. quadrants(v0,1)) then
           verts(v,1) = 1.0
           verts(v,2) = v2
           edges(e,1) = v
           edges(e,2) = v0
           v = v + 1
           e = e + 1
        end if
        if (.not. quadrants(v0,2)) then
           verts(v,1) = v1
           verts(v,2) = 1.0
           edges(e,1) = v
           edges(e,2) = v0
           v = v + 1
           e = e + 1
        end if
        if (.not. quadrants(v0,3)) then
           verts(v,1) = 0.0
           verts(v,2) = v2
           edges(e,1) = v
           edges(e,2) = v0
           v = v + 1
           e = e + 1
        end if
        if (.not. quadrants(v0,4)) then
           verts(v,1) = v1
           verts(v,2) = 0.0
           edges(e,1) = v
           edges(e,2) = v0
           v = v + 1
           e = e + 1
        end if
     end if
  end do

end subroutine addConnectors



subroutine countConnectors(nvert, nedge, Lx, Ly, verts, edges, count, quadrants)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvert, nedge, Lx, Ly, verts, edges
  !f2py intent(out) count, quadrants
  !f2py depend(nvert) verts
  !f2py depend(nedge) edges
  !f2py depend(nvert) quadrants

  !Input
  integer, intent(in) ::  nvert, nedge
  double precision, intent(in) ::  Lx, Ly
  double precision, intent(in) ::  verts(nvert,2), edges(nedge,5)

  !Output
  integer, intent(out) ::  count
  logical, intent(out) ::  quadrants(nvert,4)

  !Working
  integer v, e, d
  double precision P1(2), P2(2), t, pi
  double precision v1, v2

  pi = 2*acos(0.0)

  quadrants(:,:) = .False.
  do e=1,nedge
     do d=1,2
        if (d .eq. 1) then
           v = int(edges(e,1))
           P1 = verts(int(edges(e,1)),:)
           P2 = verts(int(edges(e,2)),:)
        else
           v = int(edges(e,2))
           P1 = verts(int(edges(e,2)),:)
           P2 = verts(int(edges(e,1)),:)
        end if
        call arc_tan(P2-P1, Lx, Ly, t)
        t = t/pi
        if ((0 .le. t) .and. (t .le. 0.25)) then
           quadrants(v,1) = .True.
        else if ((0.25 .le. t) .and. (t .le. 0.75)) then
           quadrants(v,2) = .True.
        else if ((0.75 .le. t) .and. (t .le. 1.25)) then
           quadrants(v,3) = .True.
        else if ((1.25 .le. t) .and. (t .le. 1.75)) then
           quadrants(v,4) = .True.
        else if ((1.75 .le. t) .and. (t .le. 2.00)) then
           quadrants(v,1) = .True.
        end if
     end do
  end do

  count = 0
  do v=1,nvert
     v1 = verts(v,1)
     v2 = verts(v,2)
     if (.not. ((v1.eq.0) .or. (v1.eq.1) .or. (v2.eq.0) .or. (v2.eq.1))) then
        if (.not. quadrants(v,1)) then
           count = count + 1
        end if
        if (.not. quadrants(v,2)) then
           count = count + 1
        end if
        if (.not. quadrants(v,3)) then
           count = count + 1
        end if
        if (.not. quadrants(v,4)) then
           count = count + 1
        end if
     end if
  end do

end subroutine countConnectors



subroutine arc_tan(P, Lx, Ly, t)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) P, Lx, Ly
  !f2py intent(out) t

  !Input
  double precision, intent(in) ::  P(2), Lx, Ly

  !Output
  double precision, intent(out) ::  t

  !Working
  double precision x, y, pi

  pi = 2*acos(0.0)

  x = Lx*P(1)
  y = Ly*P(2)

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

end subroutine arc_tan
