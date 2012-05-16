subroutine computeTopology(nsurf, vtol, etol, P, nvert, nedge, surf_vert, surf_edge)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nsurf, vtol, etol, P
  !f2py intent(out) nvert, nedge, surf_vert, surf_edge
  !f2py depend(nsurf) P
  !f2py depend(nsurf) surf_vert
  !f2py depend(nsurf) surf_edge

  !Input
  integer, intent(in) ::  nsurf
  double precision, intent(in) ::  vtol, etol
  double precision, intent(in) ::  P(nsurf,3,3,3)

  !Output
  integer, intent(out) ::  nvert, nedge
  integer, intent(out) ::  surf_vert(nsurf,2,2)
  integer, intent(out) ::  surf_edge(nsurf,2,2)

  !Working
  integer surf_ptrs(nsurf,3,3)
  integer surf,surf1,surf2
  integer u1,v1,u2,v2
  integer du1,dv1,du2,dv2
  integer vert, edge
  double precision dP(3), norm

  surf_ptrs(:,:,:) = 0

  vert = 0
  do surf1=1,nsurf
     do v1=1,3,2
        do u1=1,3,2
           if (surf_ptrs(surf1,u1,v1) .eq. 0) then
              vert = vert + 1
              surf_ptrs(surf1,u1,v1) = vert
              do surf2=surf1+1,nsurf
                 do v2=1,3,2
                    do u2=1,3,2
                       if (surf_ptrs(surf2,u2,v2) .eq. 0) then
                          dP(:) = P(surf1,u1,v1,:) - P(surf2,u2,v2,:)
                          norm = (dP(1)**2+dP(2)**2+dP(3)**2)**0.5
                          if (norm .lt. vtol) then
                             surf_ptrs(surf2,u2,v2) = vert
                          end if
                       end if
                    end do
                 end do
              end do
           end if
        end do
     end do
  end do

  edge = 0
  do surf1=1,nsurf
     do v1=2,2
        do u1=1,3,2
           du1 = 0
           dv1 = 1
           if (surf_ptrs(surf1,u1,v1) .eq. 0) then
              edge = edge + 1
              surf_ptrs(surf1,u1,v1) = edge
              do surf2=surf1+1,nsurf
                 do v2=2,2
                    do u2=1,3,2
                       du2 = 0
                       dv2 = 1
                       call checkEdge(nsurf,edge,surf1,surf2,u1,v1,u2,v2,du1,dv1,du2,dv2,etol,P,surf_ptrs)
                    end do
                 end do
                 do v2=1,3,2
                    do u2=2,2
                       du2 = 1
                       dv2 = 0
                       call checkEdge(nsurf,edge,surf1,surf2,u1,v1,u2,v2,du1,dv1,du2,dv2,etol,P,surf_ptrs)
                    end do
                 end do
              end do
           end if
        end do
     end do
     do v1=1,3,2
        do u1=2,2
           du1 = 1
           dv1 = 0
           if (surf_ptrs(surf1,u1,v1) .eq. 0) then
              edge = edge + 1
              surf_ptrs(surf1,u1,v1) = edge
              do surf2=surf1+1,nsurf
                 do v2=2,2
                    do u2=1,3,2
                       du2 = 0
                       dv2 = 1
                       call checkEdge(nsurf,edge,surf1,surf2,u1,v1,u2,v2,du1,dv1,du2,dv2,etol,P,surf_ptrs)
                    end do
                 end do
                 do v2=1,3,2
                    do u2=2,2
                       du2 = 1
                       dv2 = 0
                       call checkEdge(nsurf,edge,surf1,surf2,u1,v1,u2,v2,du1,dv1,du2,dv2,etol,P,surf_ptrs)
                    end do
                 end do
              end do
           end if
        end do
     end do
  end do

  do surf=1,nsurf
     surf_vert(surf,1,1) = surf_ptrs(surf,1,1)
     surf_vert(surf,1,2) = surf_ptrs(surf,1,3)
     surf_vert(surf,2,1) = surf_ptrs(surf,3,1)
     surf_vert(surf,2,2) = surf_ptrs(surf,3,3)
     surf_edge(surf,1,1) = surf_ptrs(surf,2,1)
     surf_edge(surf,1,2) = surf_ptrs(surf,2,3)
     surf_edge(surf,2,1) = surf_ptrs(surf,1,2)
     surf_edge(surf,2,2) = surf_ptrs(surf,3,2)
  end do

  nvert = vert
  nedge = edge

end subroutine computeTopology



subroutine checkEdge(nsurf, edge, surf1, surf2, u1, v1, u2, v2, du1, dv1, du2, dv2, etol, P, surf_ptrs)

  implicit none

  !Input
  integer, intent(in) ::  nsurf, edge, surf1, surf2, u1, v1, u2, v2, du1, dv1, du2, dv2
  double precision, intent (in) ::  etol, P(nsurf,3,3,3)

  !Output
  integer, intent (inout) ::  surf_ptrs(nsurf,3,3)

  !Working
  double precision dP(3), norm

  if (surf_ptrs(surf2,u2,v2) .eq. 0) then
     dP(:) = P(surf1,u1,v1,:) - P(surf2,u2,v2,:)
     norm = (dP(1)**2+dP(2)**2+dP(3)**2)**0.5
     if ((norm .lt. etol) .and. (surf_ptrs(surf1,u1-du1,v1-dv1) .eq. surf_ptrs(surf2,u2-du2,v2-dv2)) .and. (surf_ptrs(surf1,u1+du1,v1+dv1) .eq. surf_ptrs(surf2,u2+du2,v2+dv2))) then
        surf_ptrs(surf2,u2,v2) = edge
     else if ((norm .lt. etol) .and. (surf_ptrs(surf1,u1-du1,v1-dv1) .eq. surf_ptrs(surf2,u2+du2,v2+dv2)) .and. (surf_ptrs(surf1,u1+du1,v1+dv1) .eq. surf_ptrs(surf2,u2-du2,v2-dv2))) then
        surf_ptrs(surf2,u2,v2) = -edge
     end if
  end if

end subroutine checkEdge



subroutine countVEptrs(nsurf, nvert, nedge, surf_vert, surf_edge, vert_count, edge_count)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nsurf, nvert, nedge, surf_vert, surf_edge
  !f2py intent(out) vert_count, edge_count
  !f2py depend(nsurf) surf_vert
  !f2py depend(nsurf) surf_edge
  !f2py depend(nvert) vert_count
  !f2py depend(nedge) edge_count

  !Input
  integer, intent(in) ::  nsurf, nvert, nedge
  integer, intent(in) ::  surf_vert(nsurf,2,2)
  integer, intent(in) ::  surf_edge(nsurf,2,2)

  !Output
  integer, intent(out) ::  vert_count(nvert), edge_count(nedge)

  !Working
  integer surf,u,v

  vert_count(:) = 0
  edge_count(:) = 0

  do surf=1,nsurf
     do v=1,2
        do u=1,2
           vert_count(abs(surf_vert(surf,u,v))) = vert_count(abs(surf_vert(surf,u,v))) + 1
        end do
     end do
     do v=1,2
        do u=1,2       
           edge_count(abs(surf_edge(surf,u,v))) = edge_count(abs(surf_edge(surf,u,v))) + 1  
        end do
     end do
  end do

end subroutine countVEptrs



subroutine computeGroups(nsurf, nedge, surf_edge, ngroup, edge_group)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nsurf, nedge, surf_edge
  !f2py intent(out) ngroup, edge_group
  !f2py depend(nsurf) surf_edge
  !f2py depend(nedge) edge_group

  !Input
  integer, intent(in) ::  nsurf, nedge
  integer, intent(in) ::  surf_edge(nsurf,2,2)

  !Output
  integer, intent(out) ::  ngroup
  integer, intent(out) ::  edge_group(nedge)

  !Working
  integer edge_group0(nedge)
  integer leftGroup, rightGroup, bottomGroup, topGroup
  integer edge, edge1, edge2, surf
  logical found

  do edge=1,nedge
     edge_group0(edge) = edge
  end do

  do surf=1,nsurf
     leftGroup = edge_group0(abs(surf_edge(surf,2,1)))
     rightGroup = edge_group0(abs(surf_edge(surf,2,2)))
     bottomGroup = edge_group0(abs(surf_edge(surf,1,1)))
     topGroup = edge_group0(abs(surf_edge(surf,1,2)))
     if (leftGroup .ne. rightGroup) then
        call setAll(nedge, leftGroup, rightGroup, edge_group0)
     end if
     if (bottomGroup .ne. topGroup) then
        call setAll(nedge, bottomGroup, topGroup, edge_group0)
     end if
  end do

  edge_group(:) = 0
  ngroup = 0
  do edge1=1,nedge
     found = .False.
     do edge2=1,nedge
        if (edge_group0(edge2) .eq. edge1) then
           if (.not. found) then
              ngroup = ngroup + 1
              found = .True.
           end if
           edge_group(edge2) = ngroup
        end if
     end do
  end do

end subroutine computeGroups



subroutine setAll(nedge, oldGroup, newGroup, edge_group0)

  implicit none

  !Input
  integer, intent(in) ::  nedge, oldGroup, newGroup

  !Output
  integer, intent(inout) ::  edge_group0(nedge)

  !Working
  integer edge

  do edge=1,nedge
     if (edge_group0(edge) .eq. oldGroup) then
        edge_group0(edge) = newGroup
     end if
  end do

end subroutine setAll
