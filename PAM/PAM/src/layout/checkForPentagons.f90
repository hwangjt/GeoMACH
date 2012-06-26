subroutine hasPentagons(nvert, nedge, verts, edges, val)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvert, nedge, verts, edges
  !f2py intent(out) val
  !f2py depend(nvert) verts
  !f2py depend(nedge) edges

  !Input
  integer, intent(in) ::  nvert, nedge
  double precision, intent(in) ::  verts(nvert,2), edges(nedge,5)

  !Output
  logical, intent(out) ::  val

  !Working
  integer v0, v1, v2, e0, e1, e2
  double precision pi
  integer count, counter

  pi = 2*acos(0.0)

  counter = 0
  val = .False.
  do v0=1,nvert
     do e0=1,nedge
        e1 = e0
        call getOtherV(v0, edges(e1,1:2), v1)
        if (v1 .ne. 0) then
           turns: do count=1,3
              call turnRight(nvert, nedge, v1, e1, verts, edges, v2, e2)
              if ((v2 .eq. 0) .or. (v2 .eq. v0)) then
                 exit turns
              else if (count .eq. 3) then
                 val = .True.
                 counter = counter + 1
              else
                 v1 = v2
                 e1 = e2
              end if
           end do turns
        end if
     end do
  end do

end subroutine hasPentagons



subroutine turnRight(nvert, nedge, v1, e1, verts, edges, v2, e2)

  implicit none

  !Input
  integer, intent(in) ::  nvert, nedge
  integer, intent(in) ::  v1, e1
  double precision, intent(in) ::  verts(nvert,2), edges(nedge,5)

  !Output
  integer, intent(out) ::  v2, e2

  !Working
  integer v0, v, e, mine
  double precision t1, t, mint, pi

  pi = 2*acos(0.0)

  call getOtherV(v1, edges(e1,1:2), v0)
  call arc_tan(verts(v0,:)-verts(v1,:), 1.0, 1.0, t1)
  t1 = t1/pi
  
  mint = 4.0
  mine = 0
  do e=1,nedge
     call getOtherV(v1, edges(e,1:2), v)
     if ((v .ne. 0) .and. (e .ne. e1)) then
        call arc_tan(verts(v,:)-verts(v1,:), 1.0, 1.0, t)
        t = t/pi
        if (t .le. t1) then
           t = t + 2.0
        end if
        if (t .lt. mint) then
           mint = t
           mine = e
        end if
     end if
  end do
  if (mint - t1 .lt. 1) then
     e2 = mine
     call getOtherV(v1, edges(e2,1:2), v2)
  else
     v2 = 0
     e2 = 0
  end if

end subroutine turnRight



subroutine getOtherV(v0, edge, v)

  implicit none

  !Input
  integer, intent(in) ::  v0
  double precision, intent(in) ::  edge(2)

  !Output
  integer, intent(out) ::  v

  if (int(edge(1)) .eq. v0) then
     v = int(edge(2))
  else if (int(edge(2)) .eq. v0) then
     v = int(edge(1))
  else
     v = 0
  end if

end subroutine getOtherV
