subroutine computeFilletedC(axis, nC, nsurf, nedge, ngroup, nvert, surf_vert, surf_edge, edge_group, group_m, vert_count, edge_count, vert_symm, edge_symm, surf_index_C, edge_index_C, C)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) axis, nC, nsurf, nedge, ngroup, nvert, surf_vert, surf_edge, edge_group, group_m, vert_count, edge_count, vert_symm, edge_symm, surf_index_C, edge_index_C
  !f2py intent(inout) C
  !f2py depend(nsurf) surf_vert
  !f2py depend(nsurf) surf_edge
  !f2py depend(nedge) edge_group
  !f2py depend(ngroup) group_m
  !f2py depend(nvert) vert_count
  !f2py depend(nedge) edge_count
  !f2py depend(nvert) vert_symm
  !f2py depend(nedge) edge_symm
  !f2py depend(nsurf) surf_index_C
  !f2py depend(nedge) edge_index_C
  !f2py depend(nC) C

  !Input
  integer, intent(in) ::  axis, nC, nsurf, nedge, ngroup, nvert
  integer, intent(in) ::  surf_vert(nsurf,2,2), surf_edge(nsurf,2,2), edge_group(nedge), group_m(ngroup), vert_count(nvert), edge_count(nedge)
  logical, intent(in) ::  vert_symm(nvert), edge_symm(nedge)
  logical, intent(in) ::  surf_index_C(nsurf,2), edge_index_C(nedge,2)

  !Output
  double precision, intent(inout) ::  C(nC,3)

  !Working 
  integer surf, edge, vert
  integer k
  integer u,v
  integer iC, iC1, iC2
  integer, allocatable, dimension(:,:) ::  mapping
  integer ugroup, vgroup, m, mu, mv

  iC2 = edge_index_C(nedge,2)
  C(1:iC2+nvert,:) = 0.0

  do k=1,3
     do surf=1,nsurf
        ugroup = edge_group(abs(surf_edge(surf,1,1)))
        vgroup = edge_group(abs(surf_edge(surf,2,1)))
        mu = group_m(ugroup)
        mv = group_m(vgroup)
        allocate(mapping(mu,mv))
        call getMapping(surf, mu, mv, nsurf, nedge, ngroup, nvert, surf_vert, surf_edge, edge_group, group_m, surf_index_C, edge_index_C, mapping)
        do u=2,mu-1
           do v=2,2
              C(mapping(u,v-1),k) = C(mapping(u,v-1),k) + C(mapping(u,v),k)
           end do
           do v=mv-1,mv-1
              C(mapping(u,v+1),k) = C(mapping(u,v+1),k) + C(mapping(u,v),k)
           end do
        end do
        do v=2,mv-1
           do u=2,2
              C(mapping(u-1,v),k) = C(mapping(u-1,v),k) + C(mapping(u,v),k)
           end do
           do u=mu-1,mu-1
              C(mapping(u+1,v),k) = C(mapping(u+1,v),k) + C(mapping(u,v),k)
           end do
        end do
        C(mapping(1,1),k) = C(mapping(1,1),k) + C(mapping(2,2),k)
        C(mapping(1,mv),k) = C(mapping(1,mv),k) + C(mapping(2,mv-1),k)
        C(mapping(mu,1),k) = C(mapping(mu,1),k) + C(mapping(mu-1,2),k)
        C(mapping(mu,mv),k) = C(mapping(mu,mv),k) + C(mapping(mu-1,mv-1),k)
        deallocate(mapping)
     end do
     do vert=1,nvert
        C(vert,k) = (1.0*C(vert,k))/vert_count(vert)
        if (vert_symm(vert)) then
           C(vert,axis) = 0.0
        end if
     end do
     do edge=1,nedge
        m = group_m(edge_group(edge))
        iC1 = edge_index_C(edge,1)
        iC = nvert + iC1 + 1
        do u=1,m-2
           C(iC,k) = (1.0*C(iC,k))/edge_count(edge)
           if (edge_symm(edge)) then
              C(iC,axis) = 0.0
           end if
           iC = iC + 1
        end do
     end do
  end do

end subroutine computeFilletedC
