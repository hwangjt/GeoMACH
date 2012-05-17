subroutine getSurfaceSizes(nsurf, nedge, ngroup, surf_edge, edge_group, & 
           group_n, n)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nsurf, nedge, ngroup, surf_edge, edge_group, group_n
  !f2py intent(out) n
  !f2py depend(nsurf) surf_edge
  !f2py depend(nedge) edge_group
  !f2py depend(ngroup) group_n
  !f2py depend(nsurf) n

  !Input
  integer, intent(in) ::  nsurf, nedge, ngroup
  integer, intent(in) ::  surf_edge(nsurf,2,2), edge_group(nedge), &
                          group_n(ngroup)

  !Output
  integer, intent(out) ::  n(nsurf,2)

  !Working
  integer surf

  do surf=1,nsurf
     n(surf,1) = group_n(edge_group(abs(surf_edge(surf,1,1))))
     n(surf,2) = group_n(edge_group(abs(surf_edge(surf,2,1))))
  end do
  
end subroutine getSurfaceSizes



subroutine getSurfaceP(surf, nP, nu, nv, nsurf, nedge, ngroup, nvert, &
           surf_vert, surf_edge, edge_group, group_n, surf_index_P, &
           edge_index_P, P, surfP)

  implicit none
  
  !Fortran-python interface directives
  !f2py intent(in) surf, nP, nu, nv, nsurf, nedge, ngroup, nvert, surf_vert, surf_edge, edge_group, group_n, surf_index_P, edge_index_P, P
  !f2py intent(out) surfP
  !f2py depend(nsurf) surf_vert, surf_edge
  !f2py depend(nedge) edge_group
  !f2py depend(ngroup) group_n
  !f2py depend(nsurf) surf_index_P
  !f2py depend(nedge) edge_index_P
  !f2py depend(nP) P
  !f2py depend(nu,nv) surfP

  !Input
  integer, intent(in) ::  surf, nP, nu, nv, nsurf, nedge, ngroup, nvert
  integer, intent(in) ::  surf_vert(nsurf,2,2), surf_edge(nsurf,2,2), &
                          edge_group(nedge), group_n(ngroup), &
                          surf_index_P(nsurf,2), edge_index_P(nedge,2)
  double precision, intent(in) ::  P(nP,3)

  !Output
  double precision, intent(out) ::  surfP(nu,nv,3)

  !Working
  integer k,u,v
  integer iP1, iP2

  surfP(:,:,:) = 0.0
  do k=1,3
     do u=1,2
        do v=1,2
           surfP(1+(u-1)*(nu-1),1+(v-1)*(nv-1),k) = P(surf_vert(surf,u,v),k)
        end do
     end do     
  end do
  do k=1,3
     do v=1,2
        iP1 = edge_index_P(abs(surf_edge(surf,1,v)),1)
        iP2 = edge_index_P(abs(surf_edge(surf,1,v)),2)
        iP1 = iP1 + nvert + 1
        iP2 = iP2 + nvert
        if (surf_edge(surf,1,v) .gt. 0) then
           surfP(2:nu-1,1+(v-1)*(nv-1),k) = P(iP1:iP2,k)
        else
           surfP(2:nu-1,1+(v-1)*(nv-1),k) = P(iP2:iP1:-1,k)
        end if
     end do
     do u=1,2
        iP1 = edge_index_P(abs(surf_edge(surf,2,u)),1)
        iP2 = edge_index_P(abs(surf_edge(surf,2,u)),2)
        iP1 = iP1 + nvert + 1
        iP2 = iP2 + nvert
        if (surf_edge(surf,2,u) .gt. 0) then
           surfP(1+(u-1)*(nu-1),2:nv-1,k) = P(iP1:iP2,k)
        else
           surfP(1+(u-1)*(nu-1),2:nv-1,k) = P(iP2:iP1:-1,k)
        end if
     end do
  end do
  v = edge_index_P(nedge,2)
  iP1 = surf_index_P(surf,1)
  iP1 = iP1 + nvert + v + 1  
  do v=2,nv-1
     do u=2,nu-1
        do k=1,3
           surfP(u,v,k) = P(iP1,k)
        end do
        iP1 = iP1 + 1
     end do
  end do
  
end subroutine getSurfaceP
