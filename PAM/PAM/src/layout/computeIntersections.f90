subroutine splitEdges(nedge, nvert, nedge0, verts, edges0, edges)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nedge, nvert, nedge0, verts, edges0
  !f2py intent(out) edges
  !f2py depend(nvert) verts
  !f2py depend(nedge0) edges0
  !f2py depend(nedge) edges

  !Input
  integer, intent(in) ::  nedge, nvert, nedge0
  double precision, intent(in) ::  verts(nvert,2), edges0(nedge0,5)

  !Output
  double precision, intent(out) ::  edges(nedge,5)

  !Working
  integer v, e, e0, c, index
  integer s, v1, v2, va, vb
  double precision u, u1, u2, ua, ub
  integer count
  integer, allocatable, dimension(:) ::  vs
  double precision, allocatable, dimension(:) ::  us

  e = 1
  do e0=1,nedge0
     count = 0
     do v=1,nvert
        call split(v, e0, nvert, nedge0, verts, edges0, u)
        if (u .ge. 0) then
           count = count + 1
        end if
     end do
     if (count .eq. 0) then
        edges(e,:) = edges0(e0,:)
        e = e + 1
     else
        allocate(us(count))
        allocate(vs(count))
        count = 0
        do v=1,nvert
           call split(v, e0, nvert, nedge0, verts, edges0, u)
           if (u .ge. 0) then           
              count = count + 1
              us(count) = u
              vs(count) = v
           end if
        end do
        v1 = int(edges0(e0,1))
        v2 = int(edges0(e0,2))
        s = int(edges0(e0,3))
        u1 = edges0(e0,4)
        u2 = edges0(e0,5)
        vb = v2
        ub = 1.0
        do c=1,count+1
           if (c.eq.count+1) then
              va = v1
              ua = 0.0
           else
              index = maxloc(us,1)
              va = vs(index)
              ua = us(index)
              us(index) = 0.0
           end if
           edges(e,1) = va
           edges(e,2) = vb
           edges(e,3) = s
           edges(e,4) = u1*(1-ua) + u2*ua
           edges(e,5) = u1*(1-ub) + u2*ub
           e = e + 1
           vb = va
           ub = ua
        end do
        deallocate(us)  
        deallocate(vs)        
     end if
  end do

end subroutine splitEdges



subroutine numSplits(nvert, nedge, verts, edges, nsplit)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvert, nedge, verts, edges
  !f2py intent(out) nsplit
  !f2py depend(nvert) verts
  !f2py depend(nedge) edges

  !Input
  integer, intent(in) ::  nvert, nedge
  double precision, intent(in) ::  verts(nvert,2), edges(nedge,5)

  !Output
  integer, intent(out) ::  nsplit

  !Working
  integer v, e
  double precision u

  nsplit = 0
  do e=1,nedge
     do v=1,nvert
        call split(v, e, nvert, nedge, verts, edges, u)
        if (u .ge. 0) then
           nsplit = nsplit + 1
        end if
     end do
  end do

end subroutine numSplits



subroutine split(v, e, nvert, nedge, verts, edges, u)

  implicit none

  !Input
  integer, intent(in) ::  v, e, nvert, nedge
  double precision, intent(in) ::  verts(nvert,2), edges(nedge,5)

  !Output
  double precision, intent(out) ::  u

  !Working
  integer v1, v2
  double precision r1(2), r2(2), L1, L2

  u = -1
  v1 = int(edges(e,1))
  v2 = int(edges(e,2))
  if ((v .ne. v1) .and. (v .ne. v2)) then
     r1 = (verts(v,:) - verts(v1,:))
     r2 = (verts(v,:) - verts(v2,:))
     L1 = dot_product(r1,r1)**0.5
     L2 = dot_product(r2,r2)**0.5
     r1 = r1/L1
     r2 = r2/L2
     if (abs(dot_product(r1,r2)+1) .lt. 1e-14) then
        u = L1/(L1+L2)
     end if
  end if

end subroutine split



subroutine addIntersections(nvert, nvert0, nedge, verts0, edges, verts)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvert, nvert0, nedge, verts0, edges
  !f2py intent(out) verts
  !f2py depend(nvert0) verts0
  !f2py depend(nedge) edges
  !f2py depend(nvert) verts
  
  !Input
  integer, intent(in) ::  nvert, nvert0, nedge
  double precision, intent(in) ::  verts0(nvert0,2), edges(nedge,5)

  !Output
  double precision, intent(out) ::  verts(nvert,2)

  !Working
  integer v, e1, e2
  double precision vert(2)

  v = nvert0 + 1
  verts(1:nvert0,:) = verts0(:,:)
  do e1=1,nedge
     do e2=e1+1,nedge
        call intersect(e1, e2, nvert0, nedge, verts0, edges, vert)
        if ((vert(1).ge.0) .and. (vert(2).ge.0)) then
           verts(v,:) = vert(:)
           v = v + 1
        end if
     end do
  end do

end subroutine addIntersections



subroutine numIntersections(nvert, nedge, verts, edges, nint)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvert, nedge, verts, edges
  !f2py intent(out) nint
  !f2py depend(nvert) verts
  !f2py depend(nedge) edges
  
  !Input
  integer, intent(in) ::  nvert, nedge
  double precision, intent(in) ::  verts(nvert,2), edges(nedge,5)

  !Output
  integer, intent(out) ::  nint

  !Working
  integer e1, e2
  double precision vert(2)

  nint = 0
  do e1=1,nedge
     do e2=e1+1,nedge
        call intersect(e1, e2, nvert, nedge, verts, edges, vert)
        if ((vert(1).ge.0) .and. (vert(2).ge.0)) then
           nint = nint + 1
        end if
     end do
  end do

end subroutine numIntersections



subroutine intersect(e1, e2, nvert, nedge, verts, edges, vert)

  implicit none

  !Input
  integer, intent(in) ::  e1, e2, nvert, nedge
  double precision, intent(in) ::  verts(nvert,2), edges(nedge,5)

  !Output
  double precision, intent(out) :: vert(2)

  !Working
  double precision a(2), b(2), u(2), w(2), bma(2), s, t, det

  vert(:) = -1
  a = verts(int(edges(e1,1)),:)
  u = verts(int(edges(e1,2)),:) - verts(int(edges(e1,1)),:)
  b = verts(int(edges(e2,1)),:)
  w = verts(int(edges(e2,2)),:) - verts(int(edges(e2,1)),:)
  det = u(2)*w(1) - u(1)*w(2)
  bma = b - a
  if (det .ne. 0) then
     s = 1/det*(-w(2)*bma(1) + w(1)*bma(2))
     t = 1/det*(-u(2)*bma(1) + u(1)*bma(2))
     if ((0.lt.s).and.(s.lt.1).and.(0.lt.t).and.(t.lt.1)) then
        vert = a + u*s
     end if
  end if

end subroutine intersect
