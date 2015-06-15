subroutine computeConstrainedVerts(nvert, nedge, edges, edgeCon, vertCon)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvert, nedge, edges, edgeCon
  !f2py intent(out) vertCon
  !f2py depend(nvert) vertCon
  !f2py depend(nedge) edges, edgeCon

  !Input
  integer, intent(in) ::  nvert, nedge, edges(nedge,2)
  logical, intent(in) ::  edgeCon(nedge)

  !Output
  logical, intent(out) ::  vertCon(nvert)

  !Working
  integer iedge

  vertCon(:) = .False.
  do iedge=1,nedge
     if (edgeCon(iedge)) then
        vertCon(edges(iedge,1)) = .True.
        vertCon(edges(iedge,2)) = .True.
     end if
  end do

end subroutine computeConstrainedVerts




function sameVert(A, B)

  implicit none
  double precision, intent(in) ::  A(2), B(2)
  logical sameVert

  sameVert = sqrt(dot_product(B-A,B-A)) .lt. 1e-10

end function sameVert




subroutine computeConstrainedEdges(nedge0, nedge, edges0, edges, edgeCon)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nedge0, nedge, edges0, edges
  !f2py intent(out) edgeCon
  !f2py depend(nedge0) edges0
  !f2py depend(nedge) edges, edgeCon

  !Input
  integer, intent(in) ::  nedge0, nedge, edges0(nedge0,2), edges(nedge,2)

  !Output
  logical, intent(out) ::  edgeCon(nedge)

  !Working
  integer iedge0, iedge, sameEdge

  edgeCon(:) = .False.
  do iedge=1,nedge
     do iedge0=1,nedge0
        if (sameEdge(edges(iedge,:), edges0(iedge0,:)) .ne. 0) then
           edgeCon(iedge) = .True.
        end if
     end do
  end do

end subroutine computeConstrainedEdges



function sameEdge(A,B)

  implicit none
  integer, intent(in) ::  A(2), B(2)
  integer sameEdge

  if ((A(1) .eq. B(1)) .and. (A(2) .eq. B(2))) then
     sameEdge = 1
  else if ((A(1) .eq. B(2)) .and. (A(2) .eq. B(1))) then
     sameEdge = -1
  else
     sameEdge = 0
  end if

end function sameEdge
