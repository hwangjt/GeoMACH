subroutine addEdgePts(nvert0, nedge, nvert, maxL, &
     edgeLengths, verts0, edges, verts)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvert0, nedge, nvert, maxL, edgeLengths, verts0, edges
  !f2py intent(out) verts
  !f2py depend(nvert0) verts0
  !f2py depend(nedge) edges
  !f2py depend(nvert) verts

  !Input
  integer, intent(in) ::  nvert0, nedge, nvert
  double precision, intent(in) ::  maxL, edgeLengths(2,2), verts0(nvert0,2)
  integer, intent(in) ::  edges(nedge,2)

  !Output
  double precision, intent(out) ::  verts(nvert,2)

  !Working
  logical validSplit, found, collinear(nvert0)
  integer iedge, ivert, ivert0, i, n
  double precision v(2), v1(2), v2(2), wtdDist, t, den
  double precision Lx, Ly, Lx0, Lx1, Ly0, Ly1

  verts(:,:) = 0.
  verts(1:nvert0,:) = verts0(:,:)

  Lx0 = edgeLengths(1,1)
  Lx1 = edgeLengths(1,2)
  Ly0 = edgeLengths(2,1)
  Ly1 = edgeLengths(2,2)

  ivert = nvert0
  do iedge=1,nedge
     v1 = verts0(edges(iedge,1),:)
     v2 = verts0(edges(iedge,2),:)
     v = 0.5 * (v1 + v2)
     Lx = Lx0*(1-v(2)) + Lx1*v(2)
     Ly = Ly0*(1-v(1)) + Ly1*v(1)
     n = floor(wtdDist(Lx,Ly,v1,v2)/maxL)
     if (n .ne. 0) then
        do ivert0=1,nvert0
           collinear(ivert0) = validSplit(v1, v2, verts0(ivert0,:))
        end do
        den = 1.0/(n+1)
        do i=1,n
           t = i * den
           v = v1 + t*(v2-v1)
           Lx = Lx0*(1-v(2)) + Lx1*v(2)
           Ly = Ly0*(1-v(1)) + Ly1*v(1)
           found = .False.
           do ivert0=1,nvert0
              if (collinear(ivert0)) then
                 if (wtdDist(Lx,Ly,verts0(ivert0,:),v) .lt. maxL/2) then
                    found = .True.
                 end if
              end if
           end do
           if (.not. found) then
              ivert = ivert + 1
              verts(ivert,:) = v
           end if
        end do
     end if
  end do

end subroutine addEdgePts



subroutine countEdgePts(nvert, nedge, maxL, &
     edgeLengths, verts, edges, npts)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvert, nedge, maxL, edgeLengths, verts, edges
  !f2py intent(out) npts
  !f2py depend(nvert) verts
  !f2py depend(nedge) edges

  !Input
  integer, intent(in) ::  nvert, nedge
  double precision, intent(in) ::  maxL, edgeLengths(2,2), verts(nvert,2)
  integer, intent(in) ::  edges(nedge,2)

  !Output
  integer, intent(out) ::  npts

  !Working
  integer iedge
  double precision v(2), v1(2), v2(2), wtdDist
  double precision Lx, Ly, Lx0, Lx1, Ly0, Ly1

  Lx0 = edgeLengths(1,1)
  Lx1 = edgeLengths(1,2)
  Ly0 = edgeLengths(2,1)
  Ly1 = edgeLengths(2,2)

  npts = 0
  do iedge=1,nedge
     v1 = verts(edges(iedge,1),:)
     v2 = verts(edges(iedge,2),:)
     v = 0.5 * (v1 + v2)
     Lx = Lx0*(1-v(2)) + Lx1*v(2)
     Ly = Ly0*(1-v(1)) + Ly1*v(1)
     npts = npts + floor(wtdDist(Lx,Ly,v1,v2)/maxL)
  end do

end subroutine countEdgePts




function wtdDist(Lx, Ly, A, B)

  implicit none

  !Input
  double precision, intent(in) ::  Lx, Ly, A(2), B(2)

  !Output
  double precision wtdDist

  !Working
  double precision V(2)

  V(1) = Lx*(B(1) - A(1))
  V(2) = Ly*(B(2) - A(2))  
  wtdDist = sqrt(dot_product(V,V))

end function wtdDist
