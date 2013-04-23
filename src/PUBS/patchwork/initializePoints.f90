subroutine populateP(nP, n1, n2, surf, nsurf, nedge, ngroup, nvert, surf_vert, &
           surf_edge, edge_group, group_n, & 
           surf_index_P, edge_index_P, P0, P)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nP, n1, n2, surf, nsurf, nedge, ngroup, nvert, surf_vert, surf_edge, edge_group, group_n, surf_index_P, edge_index_P, P0
  !f2py intent(inout) P
  !f2py depend(nsurf) surf_vert
  !f2py depend(nsurf) surf_edge
  !f2py depend(nedge) edge_group
  !f2py depend(ngroup) group_n
  !f2py depend(nsurf) surf_index_P
  !f2py depend(nedge) edge_index_P
  !f2py depend(n1,n2) P0
  !f2py depend(nP) P

  !Input
  integer, intent(in) ::  nP, n1, n2, surf, nsurf, nedge, ngroup, nvert
  integer, intent(in) ::  surf_vert(nsurf,2,2), surf_edge(nsurf,2,2), & 
                          edge_group(nedge)
  integer, intent(in) ::  group_n(ngroup)
  integer, intent(in) ::  surf_index_P(nsurf,2), edge_index_P(nedge,2)
  double precision, intent(in) ::  P0(n1,n2,3)

  !Output
  double precision, intent(inout) ::  P(nP,3)

  !Working
  integer k
  integer u,v,w
  integer vert,edge,n,iP,iP1,iP2,dir

  do k=1,3
     do u=1,2
        do v=1,2
           vert = surf_vert(surf,u,v)
           P(vert,k) = P(vert,k) + P0(1+(n1-1)*(u-1),1+(n2-1)*(v-1),k)
        end do
     end do
     do u=1,2
        do v=1,2
           edge = abs(surf_edge(surf,u,v))
           n = group_n(edge_group(edge))
           iP1 = edge_index_P(edge,1)
           iP2 = edge_index_P(edge,2)
           if (surf_edge(surf,u,v) .gt. 0) then
              iP = nvert + iP1 + 1
              dir = 1
           else
              iP = nvert + iP2
              dir = -1
           end if
           do w=1,n-2
              if (u.eq.1) then
                 if (v.eq.1) then
                    P(iP,k) = P(iP,k) + P0(w+1,1,k)
                 else
                    P(iP,k) = P(iP,k) + P0(w+1,n2,k)
                 end if
              else
                 if (v.eq.1) then
                    P(iP,k) = P(iP,k) + P0(1,w+1,k)
                 else
                    P(iP,k) = P(iP,k) + P0(n1,w+1,k)
                 end if
              end if
              iP = iP + dir
           end do
        end do
     end do   
     iP2 = edge_index_P(nedge,2)
     iP = nvert + iP2
     iP1 = surf_index_P(surf,1)
     iP = iP + iP1 + 1
     do v=2,n2-1
        do u=2,n1-1
           P(iP,k) = P0(u,v,k)
           iP = iP + 1
        end do
     end do
  end do

end subroutine populateP



subroutine avgBoundaries(nP, nedge, ngroup, nvert, edge_group, group_n, & 
           vert_count, edge_count, edge_index_P, P)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nP, nedge, ngroup, nvert, edge_group, group_n, vert_count, edge_count, edge_index_P
  !f2py intent(inout) P
  !f2py depend(nedge) edge_group
  !f2py depend(ngroup) group_n
  !f2py depend(nvert) vert_count
  !f2py depend(nedge) edge_count
  !f2py depend(nedge) edge_index_P
  !f2py depend(nP) P

  !Input
  integer, intent(in) ::  nP, nedge, ngroup, nvert
  integer, intent(in) ::  edge_group(nedge), group_n(ngroup), & 
           vert_count(nvert), edge_count(nedge), edge_index_P(nedge,2)

  !Output
  double precision, intent(inout) ::  P(nP,3)

  !Working
  integer k,u
  integer vert, edge
  integer n,iP,iP1

  do k=1,3
     do vert=1,nvert
        P(vert,k) = P(vert,k)/vert_count(vert)
     end do
     do edge=1,nedge
        n = group_n(edge_group(edge))
        iP1 = edge_index_P(edge,1)
        iP = nvert + iP1 + 1
        do u=1,n-2
           P(iP,k) = P(iP,k)/edge_count(edge)
           iP = iP + 1
        end do
     end do
  end do

end subroutine avgBoundaries


subroutine determineSymm(axis, vtol, etol, nP, nedge, ngroup, nvert, & 
           edge_group, group_n, P, vert_symm, edge_symm)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) axis, vtol, etol, nP, nedge, ngroup, nvert, edge_group, group_n, P
  !f2py intent(out) vert_symm, edge_symm
  !f2py depend(nedge) edge_group
  !f2py depend(ngroup) group_n
  !f2py depend(nP) P
  !f2py depend(nvert) vert_symm
  !f2py depend(nedge) edge_symm

  !Input
  integer, intent(in) ::  axis
  double precision, intent(in) ::  vtol, etol
  integer, intent(in) ::  nP, nedge, ngroup, nvert
  integer, intent(in) ::  edge_group(nedge), group_n(ngroup)
  double precision, intent(in) ::  P(nP,3)

  !Output
  logical, intent(out) ::  vert_symm(nvert), edge_symm(nedge)

  !Working
  integer vert,edge,u
  integer iP
  double precision norm

  vert_symm(:) = .False.
  edge_symm(:) = .False.

  do vert=1,nvert
     if (abs(P(vert,axis)) .lt. vtol) then
        vert_symm(vert) = .True.
     end if
  end do

  iP = nvert + 1
  do edge=1,nedge
     norm = 0.0
     do u=1,group_n(edge_group(edge))-2
        norm = norm + P(iP,axis)**2
        iP = iP + 1
     end do
     norm = norm**0.5
     if (norm .lt. etol) then
        edge_symm(edge) = .True.
     end if
  end do

end subroutine determineSymm

  
