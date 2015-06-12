subroutine computeCDT(nvert, nedge, ntri, verts, edges, ntri0, triangles)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvert, nedge, ntri, verts, edges
  !f2py intent(out) ntri0, triangles
  !f2py depend(nvert) verts
  !f2py depend(nedge) edges
  !f2py depend(ntri) triangles

  !Input
  integer, intent(in) ::  nvert, nedge, ntri
  double precision, intent(in) ::  verts(nvert,2)
  integer, intent(in) ::  edges(nedge,2)

  !Output
  integer, intent(out) ::  ntri0, triangles(ntri,3)

  !Working
  integer ier

  !Working - trmesh
  double precision x(nvert), y(nvert), dist(nvert)
  integer list(6*nvert-12), lptr(6*nvert-12), lend(nvert), lnew
  integer near(nvert), next(nvert)

  !Working - edge
  integer lwk, iwk(2,nvert-3)

  !Working - trlist
  integer ncc, lcc(1), nrow
  integer nt, ltri(6,ntri), lct(1)
  integer iedge, itri

  x(:) = verts(:,1)
  y(:) = verts(:,2)

  call trmesh(nvert, x, y, list, lptr, lend, lnew, near, next, dist, ier)
  if (ier .ne. 0) then
     print *, 'CDT trmesh error:', ier
     call exit(1)
  end if

  do iedge=1,nedge
     lwk = nvert-3
     call edge(edges(iedge,1), edges(iedge,2), x, y, lwk, iwk, list, lptr, lend, ier)
     if (ier .ne. 0) then
        print *, 'CDT edge error:', ier
        call exit(1)
     end if
  end do

  ncc = 0
  lcc(1) = 0
  nrow = 6
  call trlist(ncc, lcc, nvert, list, lptr, lend, nrow, nt, ltri, lct, ier)
  if (ier .ne. 0) then
     print *, 'CDT trlist error:', ier
     call exit(1)
  end if

  ntri0=nt
  triangles(:,:) = 0
  do itri=1,ntri0
     triangles(itri,:) = ltri(1:3,itri)
  end do

end subroutine computeCDT
