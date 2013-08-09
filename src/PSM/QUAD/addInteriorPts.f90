subroutine addInteriorPts(nvert0, nvert, maxL, edgeLengths, verts0, verts)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvert0, nvert, maxL, edgeLengths, verts0
  !f2py intent(out) verts
  !f2py depend(nvert0) verts0
  !f2py depend(nvert) verts

  !Input
  integer, intent(in) ::  nvert0, nvert
  double precision, intent(in) ::  maxL, edgeLengths(2,2), verts0(nvert0,2)

  !Output
  double precision ,intent(out) ::  verts(nvert,2)

  !Working
  logical found
  double precision Lx0, Lx1, Ly0, Ly1, Lx, Ly
  double precision denx, deny, v(2), wtdDist
  double precision nx0, nx1, ny0, ny1
  integer i, j, nx, ny, ivert, ivert0

  Lx0 = edgeLengths(1,1)
  Lx1 = edgeLengths(1,2)
  Ly0 = edgeLengths(2,1)
  Ly1 = edgeLengths(2,2)

  verts(:,:) = 0.
  verts(1:nvert0,:) = verts0(:,:)

  nx0 = Lx0/maxL
  nx1 = Lx1/maxL
  ny0 = Ly0/maxL
  ny1 = Ly1/maxL
  ny = floor(0.5*ny0 + 0.5*ny1)

  ivert = nvert0
  deny = 1.0/(ny+1)
  do j=1,ny
     nx = floor(nx0 + (nx1-nx0)*j*deny)
     denx = 1.0/(nx+1)
     do i=1,nx
        v(1) = i*denx
        v(2) = j*deny
        Lx = Lx0*(1-v(1)) + Lx1*v(1)
        Ly = Ly0*(1-v(2)) + Ly1*v(2)
        found = .False.
        do ivert0=1,nvert0
           if (wtdDist(Lx,Ly,verts0(ivert0,:),v) .lt. maxL/2) then
              found = .True.
           end if
        end do
        if (.not. found) then
           ivert = ivert + 1
           verts(ivert,:) = v
        end if
     end do
  end do

end subroutine addInteriorPts




subroutine countInteriorPts(maxL, edgeLengths, npts)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) maxL, edgeLengths
  !f2py intent(out) npts

  !Input
  double precision, intent(in) ::  maxL, edgeLengths(2,2)

  !Output
  integer, intent(out) ::  npts

  !Working
  double precision Lx0, Lx1, Ly0, Ly1, deny
  double precision nx0, nx1, ny0, ny1
  integer i, j, nx, ny

  Lx0 = edgeLengths(1,1)
  Lx1 = edgeLengths(1,2)
  Ly0 = edgeLengths(2,1)
  Ly1 = edgeLengths(2,2)

  nx0 = Lx0/maxL
  nx1 = Lx1/maxL
  ny0 = Ly0/maxL
  ny1 = Ly1/maxL
  ny = floor(0.5*ny0 + 0.5*ny1)

  deny = 1.0/(ny+1)

  npts = 0
  do j=1,ny
     nx = floor(nx0 + (nx1-nx0)*j*deny)
     do i=1,nx
        npts = npts + 1
     end do
  end do

end subroutine countInteriorPts




function dist(a,b)

  implicit none

  !Input
  double precision, intent(in) ::  a(2), b(2)

  !Output
  double precision dist

  dist = sqrt((b(1)-a(1))**2 + (b(2)-a(2))**2)

end function dist
