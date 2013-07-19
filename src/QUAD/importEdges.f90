subroutine importEdges(nvert, nedge, edges0, verts, edges)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvert, nedge, edges0
  !f2py intent(out) verts, edges
  !f2py depend(nvert) verts
  !f2py depend(nedge) edges0, edges

  !Input
  integer, intent(in) ::  nvert, nedge
  double precision, intent(in) ::  edges0(nedge,2,2)

  !Output
  double precision, intent(out) ::  verts(nvert,2)
  integer, intent(out) ::  edges(nedge,2)

  !Working
  integer i, j

  do i=1,nedge
     do j=1,2
        verts(2*(i-1)+j,:) = edges0(i,j,:)
        edges(i,j) = 2*(i-1) + j
     end do
  end do

end subroutine importEdges
