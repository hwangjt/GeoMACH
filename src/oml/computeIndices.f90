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
  integer edge
  integer i1, i2

  i2 = 0
  do edge=1,nedge
     i1 = i2
     i2 = i2 + group_n(edge_group(edge)) - 2
     edge_index(edge,1) = i1
     edge_index(edge,2) = i2
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
  integer surf
  integer i1, i2

  i2 = 0
  do surf=1,nsurf
     i1 = i2
     i2 = i2 + (group_n(edge_group(abs(surf_edge(surf,1,1)))) - 2)*(group_n(edge_group(abs(surf_edge(surf,2,1)))) - 2)
     surf_index(surf,1) = i1
     surf_index(surf,2) = i2
  end do

end subroutine getSurfIndices



subroutine getEdgeIndicesQ(nsurf, nedge, ngroup, surf_edge, edge_group, group_n, surf_c1, edge_index)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nsurf, nedge, ngroup, surf_edge, edge_group, group_n, surf_c1
  !f2py intent(out) edge_index
  !f2py depend(nsurf) surf_edge
  !f2py depend(nedge) edge_group
  !f2py depend(ngroup) group_n
  !f2py depend(nsurf) surf_c1
  !f2py depend(nedge) edge_index

  !Input
  integer, intent(in) ::  nsurf, nedge, ngroup
  integer, intent(in) ::  surf_edge(nsurf,2,2), edge_group(nedge), group_n(ngroup)
  logical, intent(in) ::  surf_c1(nsurf,3,3)

  !Output
  integer, intent(out) ::  edge_index(nedge,2)

  !Working
  integer surf,edge
  integer i1, i2
  logical dof

  edge_index(:,:) = 0
  i2 = 0
  do edge=1,nedge
     dof = .true.
     do surf=1,nsurf
        if ((surf_edge(surf,1,1) .eq. edge) .and. surf_c1(surf,2,1)) then
           dof = .false.
        end if
        if ((surf_edge(surf,1,2) .eq. edge) .and. surf_c1(surf,2,3)) then
           dof = .false.
        end if
        if ((surf_edge(surf,2,1) .eq. edge) .and. surf_c1(surf,1,2)) then
           dof = .false.
        end if
        if ((surf_edge(surf,2,2) .eq. edge) .and. surf_c1(surf,3,2)) then
           dof = .false.
        end if
     end do
     if (dof) then
        i1 = i2
        i2 = i2 + group_n(edge_group(edge)) - 2
        edge_index(edge,1) = i1
        edge_index(edge,2) = i2
     end if
  end do

end subroutine getEdgeIndicesQ



subroutine getVertIndicesQ(nsurf, nedge, nvert, surf_vert, surf_edge, surf_c1, edge_c1, vert_index)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nsurf, nedge, nvert, surf_vert, surf_edge, surf_c1, edge_c1
  !f2py intent(out) vert_index
  !f2py depend(nsurf) surf_vert, surf_edge
  !f2py depend(nsurf) surf_c1
  !f2py depend(nedge) edge_c1
  !f2py depend(nvert) vert_index

  !Input  
  integer, intent(in) ::  nsurf, nedge, nvert
  integer, intent(in) ::  surf_vert(nsurf,2,2), surf_edge(nsurf,2,2)
  logical, intent(in) ::  surf_c1(nsurf,3,3), edge_c1(nedge,2)

  !Output
  integer, intent(out) ::  vert_index(nvert)

  !Working
  integer vert, surf
  integer i,j
  logical dof(nvert)

  dof(:) = .true.
  do surf=1,nsurf
     do i=1,2
        do j=1,2
           if (surf_c1(surf,2*i-1,2*j-1)) then
              dof(surf_vert(surf,i,j)) = .false.
           end if
        end do
     end do
     do i=1,2
        do j=1,2
           if ((edge_c1(abs(surf_edge(surf,1,i)),j)) .and. (surf_edge(surf,1,i) .gt. 0)) then
              dof(surf_vert(surf,j,i)) = .false.
           end if
           if ((edge_c1(abs(surf_edge(surf,1,i)),j)) .and. (surf_edge(surf,1,i) .lt. 0)) then
              dof(surf_vert(surf,3-j,i)) = .false.
           end if
           if ((edge_c1(abs(surf_edge(surf,2,i)),j)) .and. (surf_edge(surf,2,i) .gt. 0)) then
              dof(surf_vert(surf,i,j)) = .false.
           end if
           if ((edge_c1(abs(surf_edge(surf,2,i)),j)) .and. (surf_edge(surf,2,i) .lt. 0)) then
              dof(surf_vert(surf,i,3-j)) = .false.
           end if
        end do
     end do
  end do
  
  vert_index(:) = 0
  i = 1
  do vert=1,nvert
     if (dof(vert)) then
        vert_index(vert) = i
        i = i + 1
     end if
  end do

end subroutine getVertIndicesQ
