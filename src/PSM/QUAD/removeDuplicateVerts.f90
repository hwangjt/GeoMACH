subroutine removeDuplicateVerts(nvert0, nvert, nedge, &
     ids, verts0, edges0, verts, edges)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvert0, nvert, nedge, ids, verts0, edges0
  !f2py intent(out) verts, edges
  !f2py depend(nvert0) ids, verts0
  !f2py depend(nedge) edges0
  !f2py depend(nvert) verts
  !f2py depend(nedge) edges

  !Input
  integer, intent(in) ::  nvert0, nvert, nedge, ids(nvert0)
  double precision, intent(in) ::  verts0(nvert0,2)
  integer, intent(in) ::  edges0(nedge,2)

  !Output
  double precision, intent(out) ::  verts(nvert,2)
  integer, intent(out) ::  edges(nedge,2)

  !Working
  integer i, j

  do i=1,nvert0
     verts(ids(i),:) = verts0(i,:)
  end do

  do i=1,nedge
     do j=1,2
        edges(i,j) = ids(edges0(i,j))
     end do
  end do

end subroutine removeDuplicateVerts




subroutine computeUniqueVerts(nvert, verts, nid, ids)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvert, verts
  !f2py intent(out) nid, ids
  !f2py depend(nvert) verts, ids

  !Input
  integer, intent(in) ::  nvert
  double precision, intent(in) ::  verts(nvert,2)

  !Output
  integer, intent(out) ::  nid, ids(nvert)

  !Working
  integer i1, i2
  double precision d(2)

  nid = 0
  ids(:) = 0
  do i1=1,nvert
     if (ids(i1) .eq. 0) then
        nid = nid + 1
        ids(i1) = nid
        do i2=i1+1,nvert
           if (ids(i2) .eq. 0) then
              d = verts(i1,:) - verts(i2,:)
              if (dot_product(d,d) .lt. 1e-12) then
                 ids(i2) = ids(i1)
              end if
           end if
        end do
     end if
  end do

end subroutine computeUniqueVerts
