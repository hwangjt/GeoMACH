subroutine removeDegenerateEdges(nedge0, nedge, edges0, edges)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nedge0, nedge, edges0
  !f2py intent(out) edges
  !f2py depend(nedge0) edges0
  !f2py depend(nedge) edges

  !Input
  integer, intent(in) ::  nedge0, nedge, edges0(nedge0,2)

  !Output
  integer, intent(out) ::  edges(nedge,2)

  !Working
  integer iedge0, iedge

  iedge = 0
  do iedge0=1,nedge0
     if (edges0(iedge0,1) .ne. edges0(iedge0,2)) then
        iedge = iedge + 1
        edges(iedge,:) = edges0(iedge0,:)
     end if
  end do

end subroutine removeDegenerateEdges



subroutine countDegenerateEdges(nedge, edges, ndeg)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nedge, edges
  !f2py intent(out) ndeg
  !f2py depend(nedge) edges

  !Input
  integer, intent(in) ::  nedge, edges(nedge,2)

  !Output
  integer, intent(out) ::  ndeg

  !Working
  integer iedge

  ndeg = 0
  do iedge=1,nedge
     if (edges(iedge,1) .eq. edges(iedge,2)) then
        ndeg = ndeg + 1
     end if
  end do

end subroutine countDegenerateEdges
