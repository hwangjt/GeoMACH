subroutine deleteDuplicatePolygons(npoly, npoly0, poly_vert0, poly_edge0, poly_vert, poly_edge)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) npoly, npoly0, poly_vert0, poly_edge0
  !f2py intent(out) poly_vert, poly_edge
  !f2py depend(npoly0) poly_vert0, poly_edge0
  !f2py depend(npoly) poly_vert, poly_edge

  !Input
  integer, intent(in) ::  npoly, npoly0
  integer, intent(in) ::  poly_vert0(npoly0,5), poly_edge0(npoly0,5)
  
  !Output
  integer, intent(out) ::  poly_vert(npoly,5), poly_edge(npoly,5)

  !Working
  integer p1, p2
  integer polyID(npoly0), index
  logical same

  polyID(:) = 0
  
  index = 0
  do p1=1,npoly0
     if (polyID(p1) .eq. 0) then
        index = index + 1
        polyID(p1) = index
        poly_vert(index,:) = poly_vert0(p1,:)
        poly_edge(index,:) = poly_edge0(p1,:)
        do p2=p1+1,npoly0
           if (polyID(p2) .eq. 0) then
              call comparePolygons(poly_vert0(p1,:), poly_vert0(p2,:), same)
              if (same) then
                 polyID(p2) = index
              end if
           end if
        end do
     end if
  end do

end subroutine deleteDuplicatePolygons



subroutine comparePolygons(A, B, same)

  implicit none

  !Input
  integer, intent(in) ::  A(5), B(5)

  !Output
  logical, intent(out) ::  same

  !Working
  integer sizeA, sizeB
  integer k

  sizeA = 0
  sizeB = 0
  do k=1,5
     if (A(k) .ne. 0) then
        sizeA = sizeA + 1
     end if
     if (B(k) .ne. 0) then
        sizeB = sizeB + 1
     end if
  end do

  same = .False.
  if ((sizeA.eq.3) .and. (sizeB.eq.3)) then        
     if ((A(1).eq.B(1)).and.(A(2).eq.B(2)).and.(A(3).eq.B(3))) then
        same = .True.
     else if ((A(1).eq.B(2)).and.(A(2).eq.B(3)).and.(A(3).eq.B(1))) then
        same = .True.
     else if ((A(1).eq.B(3)).and.(A(2).eq.B(1)).and.(A(3).eq.B(2))) then
        same = .True.
     end if
  else if ((sizeA.eq.4) .and. (sizeB.eq.4)) then
     if ((A(1).eq.B(1)).and.(A(2).eq.B(2)) & 
          .and.(A(3).eq.B(3)).and.(A(4).eq.B(4))) then
        same = .True.
     else if ((A(1).eq.B(2)).and.(A(2).eq.B(3)) & 
          .and.(A(3).eq.B(4)).and.(A(4).eq.B(1))) then
        same = .True.
     else if ((A(1).eq.B(3)).and.(A(2).eq.B(4)) & 
          .and.(A(3).eq.B(1)).and.(A(4).eq.B(2))) then
        same = .True.
     else if ((A(1).eq.B(4)).and.(A(2).eq.B(1)) & 
          .and.(A(3).eq.B(2)).and.(A(4).eq.B(3))) then
        same = .True.
     end if
  else if ((sizeA.eq.5) .and. (sizeB.eq.5)) then
     if ((A(1).eq.B(1)).and.(A(2).eq.B(2)).and.(A(3).eq.B(3)) & 
          .and.(A(4).eq.B(4)).and.(A(5).eq.B(5))) then
        same = .True.
     else if ((A(1).eq.B(2)).and.(A(2).eq.B(3)).and.(A(3).eq.B(4)) & 
          .and.(A(4).eq.B(5)).and.(A(5).eq.B(1))) then
        same = .True.
     else if ((A(1).eq.B(3)).and.(A(2).eq.B(4)).and.(A(3).eq.B(5)) & 
          .and.(A(4).eq.B(1)).and.(A(5).eq.B(2))) then
        same = .True.
     else if ((A(1).eq.B(4)).and.(A(2).eq.B(5)).and.(A(3).eq.B(1)) & 
          .and.(A(4).eq.B(2)).and.(A(5).eq.B(3))) then
        same = .True.
     else if ((A(1).eq.B(5)).and.(A(2).eq.B(1)).and.(A(3).eq.B(2)) & 
          .and.(A(4).eq.B(3)).and.(A(5).eq.B(4))) then
        same = .True.
     end if
  end if

end subroutine comparePolygons



subroutine computePolygons(npoly, nvert, nedge, Lx, Ly, &
     verts, edges, poly_vert, poly_edge)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) npoly, nvert, nedge, Lx, Ly, verts, edges
  !f2py intent(out) poly_vert, poly_edge
  !f2py depend(nvert) verts
  !f2py depend(nedge) edges
  !f2py depend(npoly) poly_vert, poly_edge

  !Input
  integer, intent(in) ::  npoly, nvert, nedge
  double precision, intent(in) ::  Lx, Ly
  double precision, intent(in) ::  verts(nvert,2), edges(nedge,5)

  !Output
  integer, intent(out) ::  poly_vert(npoly,5), poly_edge(npoly,5)

  !Working
  integer v0, v1, v2, e0, e1, e2
  integer count, p

  poly_vert(:,:) = 0
  poly_edge(:,:) = 0

  p = 1
  do v0=1,nvert
     do e0=1,nedge
        e1 = e0
        call getOtherV(v0, edges(e1,1:2), v1)
        if ((v1 .ne. 0) .and. (p .le. npoly)) then
           poly_vert(p,1) = v1
           poly_edge(p,1) = e1
           turns: do count=1,4
              call turnRight(nvert, nedge, v1, e1, Lx, Ly, verts, edges, v2, e2)
              poly_vert(p,count+1) = v2
              poly_edge(p,count+1) = e2
              if (v2 .eq. 0) then
                 exit turns
              else if (v2 .eq. v0) then
                 p = p + 1
                 exit turns
              else if (count .eq. 4) then
                 print *, 'Error: more than 5 sides'
                 exit turns
              else
                 v1 = v2
                 e1 = e2
              end if
           end do turns
        end if
     end do
  end do

end subroutine computePolygons



subroutine countPolygons(nvert, nedge, Lx, Ly, verts, edges, npent, nquad, ntri)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvert, nedge, Lx, Ly, verts, edges
  !f2py intent(out) npent, nquad, ntri
  !f2py depend(nvert) verts
  !f2py depend(nedge) edges

  !Input
  integer, intent(in) ::  nvert, nedge
  double precision, intent(in) ::  Lx, Ly
  double precision, intent(in) ::  verts(nvert,2), edges(nedge,5)

  !Output
  integer, intent(out) ::  npent, nquad, ntri

  !Working
  integer v0, v1, v2, e0, e1, e2
  integer count

  npent = 0
  nquad = 0
  ntri = 0

  do v0=1,nvert
     do e0=1,nedge
        e1 = e0
        call getOtherV(v0, edges(e1,1:2), v1)
        if (v1 .ne. 0) then
           turns: do count=1,4
              call turnRight(nvert, nedge, v1, e1, Lx, Ly, verts, edges, v2, e2)
              if (v2 .eq. 0) then
                 exit turns
              else if (v2 .eq. v0) then
                 if (count .eq. 2) then
                    ntri = ntri + 1
                 else if (count .eq. 3) then
                    nquad = nquad + 1
                 else if (count .eq. 4) then
                    npent = npent + 1
                 end if
                 exit turns
              else if (count .eq. 4) then
                 print *, 'Error: more than 5 sides'
                 exit turns
              else
                 v1 = v2
                 e1 = e2
              end if
           end do turns
        end if
     end do
  end do

  npent = int(npent/5)
  nquad = int(nquad/4)
  ntri = int(ntri/3)

end subroutine countPolygons



subroutine turnRight(nvert, nedge, v1, e1, Lx, Ly, verts, edges, v2, e2)

  implicit none

  !Input
  integer, intent(in) ::  nvert, nedge
  integer, intent(in) ::  v1, e1
  double precision, intent(in) ::  Lx, Ly
  double precision, intent(in) ::  verts(nvert,2), edges(nedge,5)

  !Output
  integer, intent(out) ::  v2, e2

  !Working
  integer v0, v, e, mine
  double precision t1, t, mint, pi

  pi = 2*acos(0.0)

  call getOtherV(v1, edges(e1,1:2), v0)
  call arc_tan(verts(v0,:)-verts(v1,:), Lx, Ly, t1)
  t1 = t1/pi
  
  mint = 4.0
  mine = 0
  do e=1,nedge
     call getOtherV(v1, edges(e,1:2), v)
     if ((v .ne. 0) .and. (e .ne. e1)) then
        call arc_tan(verts(v,:)-verts(v1,:), Lx, Ly, t)
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
