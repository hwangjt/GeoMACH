subroutine getEdgeIndices(nedge, ngroup, edge_group, group_n, edge_index)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nedge, ngroup, edge_group, group_n
  !f2py intent(out) edge_index
  !f2py depend(nedge) edge_group
  !f2py depend(ngroup) group_n
  !f2py depend(nedge) edge_index

  !Input
  integer, intent(in) ::  nedge, ngroup
  integer, intent(in) ::  edge_group(nedge), group_n(ngroup)

  !Output
  integer, intent(out) ::  edge_index(nedge,2)

  !Working
  integer edge1, edge2
  integer i1, i2

  do edge1=1,nedge
     i2 = 0
     do edge2=1,edge1
        i2 = i2 + group_n(edge_group(edge2)) - 2
     end do
     i1 = i2 - group_n(edge_group(edge1)) + 2
     edge_index(edge1,1) = i1
     edge_index(edge1,2) = i2
  end do

end subroutine getEdgeIndices



subroutine getSurfIndices(nsurf, nedge, ngroup, surf_edge, edge_group, group_n, surf_index)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nsurf, nedge, ngroup, surf_edge, edge_group, group_n
  !f2py intent(out) surf_index
  !f2py depend(nsurf) surf_edge
  !f2py depend(nedge) edge_group
  !f2py depend(ngroup) group_n
  !f2py depend(nsurf) surf_index

  !Input
  integer, intent(in) ::  nsurf, nedge, ngroup
  integer, intent(in) ::  surf_edge(nsurf,2,2), edge_group(nedge), group_n(ngroup)

  !Output
  integer, intent(out) ::  surf_index(nsurf,2)

  !Working
  integer surf1, surf2
  integer i1, i2

  do surf1=1,nsurf
     i2 = 0
     do surf2=1,surf1
        i2 = i2 + (group_n(edge_group(abs(surf_edge(surf2,1,1)))) - 2)*(group_n(edge_group(abs(surf_edge(surf2,2,1)))) - 2)
     end do
     i1 = i2 - (group_n(edge_group(abs(surf_edge(surf1,1,1)))) - 2)*(group_n(edge_group(abs(surf_edge(surf1,2,1)))) - 2)
     surf_index(surf1,1) = i1
     surf_index(surf1,2) = i2
  end do

end subroutine getSurfIndices
