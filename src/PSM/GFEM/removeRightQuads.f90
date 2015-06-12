subroutine removeRightQuads(nquad0, nquad, &
     quads0, right, quads)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nquad0, nquad, quads0, right
  !f2py intent(out) quads
  !f2py depend(nquad0) quads0, right
  !f2py depend(nquad) quads

  !Input
  integer, intent(in) ::  nquad0, nquad
  integer, intent(in) ::  quads0(nquad0,4)
  logical, intent(in) ::  right(nquad0)

  !Output
  integer, intent(out) ::  quads(nquad,4)

  !Working
  integer iquad0, iquad

  iquad = 0
  do iquad0=1,nquad0
     if (.not. right(iquad0)) then
        iquad = iquad + 1
        quads(iquad,:) = quads0(iquad0,:)
     end if
  end do

end subroutine removeRightQuads




subroutine countRightQuads(nnode, nquad, nodes, quads, right)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nnode, nquad, nodes, quads
  !f2py intent(out) right
  !f2py depend(nnode) nodes
  !f2py depend(nquad) quads, right

  !Input
  integer, intent(in) ::  nnode, nquad
  double precision, intent(in) ::  nodes(nnode,3)
  integer, intent(in) ::  quads(nquad,4)

  !Output
  logical, intent(out) ::  right(nquad)

  !Working
  integer iquad, j
  double precision z

  do iquad=1,nquad
     z = 0
     do j=1,4
        z = z + nodes(quads(iquad,j),3)
     end do
     if (z .lt. -1e-5) then
        right(iquad) = .True.
     else
        right(iquad) = .False.
     end if
  end do

end subroutine countRightQuads
