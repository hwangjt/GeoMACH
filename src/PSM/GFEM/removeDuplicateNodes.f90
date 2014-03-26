subroutine removeDuplicateNodes(nnode0, nnode, nquad, &
     ids, nodes0, quads0, nodes, quads)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nnode0, nnode, nquad, ids, nodes0, quads0
  !f2py intent(out) nodes, quads
  !f2py depend(nnode0) ids, nodes0
  !f2py depend(nquad) quads0
  !f2py depend(nnode) nodes
  !f2py depend(nquad) quads

  !Input
  integer, intent(in) ::  nnode0, nnode, nquad, ids(nnode0)
  double precision, intent(in) ::  nodes0(nnode0,3)
  integer, intent(in) ::  quads0(nquad,4)

  !Output
  double precision, intent(out) ::  nodes(nnode,3)
  integer, intent(out) ::  quads(nquad,4)

  !Working
  integer i, j

  do i=1,nnode0
     nodes(ids(i),:) = nodes0(i,:)
  end do

  do i=1,nquad
     do j=1,4
        quads(i,j) = ids(quads0(i,j))
     end do
  end do

end subroutine removeDuplicateNodes




subroutine computeUniqueNodes(nnode, nodes, nid, ids)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nnode, nodes
  !f2py intent(out) nid, ids
  !f2py depend(nnode) nodes, ids

  !Input
  integer, intent(in) ::  nnode
  double precision, intent(in) ::  nodes(nnode,3)

  !Output
  integer, intent(out) ::  nid, ids(nnode)

  !Working
  integer i1, i2
  double precision d(3)

  nid = 0
  ids(:) = 0
  do i1=1,nnode
     if (ids(i1) .eq. 0) then
        nid = nid + 1
        ids(i1) = nid
        do i2=i1+1,nnode
           if (ids(i2) .eq. 0) then
              d = nodes(i1,:) - nodes(i2,:)
              if (dot_product(d,d) .lt. 1e-5) then
                 ids(i2) = ids(i1)
              end if
           end if
        end do
     end if
  end do

end subroutine computeUniqueNodes
