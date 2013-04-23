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



subroutine getSurfaceT(surf, nu, nv, nT, nsurf, nedge, surf_edge, & 
           surf_index_P, edge_index_P, T, bufferT)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) surf, nu, nv, nT, nsurf, nedge, surf_edge, surf_index_P, edge_index_P, T
  !f2py intent(out) bufferT
  !f2py depend(nsurf) surf_edge, surf_index_P
  !f2py depend(nedge) edge_index_P
  !f2py depend(nT) T
  !f2py depend(nu,nv) bufferT

  !Input
  integer, intent(in) ::  surf, nu, nv, nT, nsurf, nedge
  integer, intent(in) ::  surf_edge(nsurf,2,2), surf_index_P(nsurf,2), & 
           edge_index_P(nedge,2)
  double precision, intent(in) ::  T(nT)

  !Output
  double precision, intent(out) ::  bufferT(nu,nv,2)

  !Working
  integer u,v
  integer iT1, iT2

  do u=1,2
     do v=1,nv
        bufferT(1+(u-1)*(nu-1),v,1) = u-1
     end do
  end do
  do u=1,nu
     do v=1,2
        bufferT(u,1+(v-1)*(nv-1),2) = v-1
     end do
  end do
  do v=1,2
     iT1 = edge_index_P(abs(surf_edge(surf,1,v)),1) + 1
     iT2 = edge_index_P(abs(surf_edge(surf,1,v)),2)
     if (surf_edge(surf,1,v) .gt. 0) then
        bufferT(2:nu-1,1+(v-1)*(nv-1),1) = T(iT1:iT2)
     else
        bufferT(2:nu-1,1+(v-1)*(nv-1),1) = 1 - T(iT2:iT1:-1)
     end if
  end do
  do u=1,2
     iT1 = edge_index_P(abs(surf_edge(surf,2,u)),1) + 1
     iT2 = edge_index_P(abs(surf_edge(surf,2,u)),2)
     if (surf_edge(surf,2,u) .gt. 0) then
        bufferT(1+(u-1)*(nu-1),2:nv-1,2) = T(iT1:iT2)
     else
        bufferT(1+(u-1)*(nu-1),2:nv-1,2) = 1 - T(iT2:iT1:-1)
     end if
  end do
  v = edge_index_P(nedge,2)
  iT1 = surf_index_P(surf,1)
  iT1 = 2*iT1
  iT1 = iT1 + v + 1
  do v=2,nv-1
     do u=2,nu-1
        bufferT(u,v,1) = T(iT1)
        iT1 = iT1 + 1
     end do
  end do
  do v=2,nv-1
     do u=2,nu-1
        bufferT(u,v,2) = T(iT1)
        iT1 = iT1 + 1
     end do
  end do

end subroutine getSurfaceT



subroutine getSurfaceP(surf, nP, nu, nv, nvar, nsurf, nedge, nvert, &
           surf_vert, surf_edge, surf_index_P, &
           edge_index_P, P, surfP)

  implicit none
  
  !Fortran-python interface directives
  !f2py intent(in) surf, nP, nu, nv, nvar, nsurf, nedge, nvert, surf_vert, surf_edge, surf_index_P, edge_index_P, P
  !f2py intent(out) surfP
  !f2py depend(nsurf) surf_vert, surf_edge
  !f2py depend(nsurf) surf_index_P
  !f2py depend(nedge) edge_index_P
  !f2py depend(nP,nvar) P
  !f2py depend(nu,nv,nvar) surfP

  !Input
  integer, intent(in) ::  surf, nP, nu, nv, nvar, nsurf, nedge, nvert
  integer, intent(in) ::  surf_vert(nsurf,2,2), surf_edge(nsurf,2,2), &
                          surf_index_P(nsurf,2), edge_index_P(nedge,2)
  double precision, intent(in) ::  P(nP,nvar)

  !Output
  double precision, intent(out) ::  surfP(nu,nv,nvar)

  !Working
  integer k,u,v
  integer iP1, iP2

  surfP(:,:,:) = 0.0
  do k=1,nvar
     do u=1,2
        do v=1,2
           surfP(1+(u-1)*(nu-1),1+(v-1)*(nv-1),k) = P(surf_vert(surf,u,v),k)
        end do
     end do     
  end do
  do k=1,nvar
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
        do k=1,nvar
           surfP(u,v,k) = P(iP1,k)
        end do
        iP1 = iP1 + 1
     end do
  end do
  
end subroutine getSurfaceP



subroutine getMapping(surf, mu, mv, nsurf, nedge, nvert, surf_vert, &
           surf_edge, surf_index_C, edge_index_C, mapping)

  implicit none

  !Input
  integer, intent(in) ::  surf, mu, mv, nsurf, nedge, nvert
  integer, intent(in) ::  surf_vert(nsurf,2,2), surf_edge(nsurf,2,2), & 
           surf_index_C(nsurf,2), edge_index_C(nedge,2)

  !Output
  integer, intent(out) ::  mapping(mu,mv)

  !Working
  integer u,v
  integer i, i1, i2
  
  mapping(1,1) = surf_vert(surf,1,1)
  mapping(1,mv) = surf_vert(surf,1,2)
  mapping(mu,1) = surf_vert(surf,2,1)
  mapping(mu,mv) = surf_vert(surf,2,2)
  do v=1,2
     i1 = edge_index_C(abs(surf_edge(surf,1,v)),1)
     i2 = edge_index_C(abs(surf_edge(surf,1,v)),2)
     if (surf_edge(surf,1,v) .gt. 0) then
        i = i1 + nvert + 1
        do u=2,mu-1
           mapping(u,1+(v-1)*(mv-1)) = i
           i = i + 1
        end do
     else
        i = i2 + nvert
        do u=2,mu-1
           mapping(u,1+(v-1)*(mv-1)) = i
           i = i - 1
        end do
     end if
  end do
  do u=1,2
     i1 = edge_index_C(abs(surf_edge(surf,2,u)),1)
     i2 = edge_index_C(abs(surf_edge(surf,2,u)),2)
     if (surf_edge(surf,2,u) .gt. 0) then
        i = i1 + nvert + 1
        do v=2,mv-1
           mapping(1+(u-1)*(mu-1),v) = i
           i = i + 1
        end do
     else
        i = i2 + nvert
        do v=2,mv-1
           mapping(1+(u-1)*(mu-1),v) = i
           i = i - 1
        end do
     end if
  end do
  v = edge_index_C(nedge,2)
  i1 = surf_index_C(surf,1)
  i = i1 + nvert + v + 1
  do v=2,mv-1
     do u=2,mu-1
        mapping(u,v) = i
        i = i + 1
     end do
  end do

end subroutine getMapping



subroutine inflateVector(ni, nj, nk, nV, V, Q)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) ni, nj, nk, nV, V
  !f2py intent(out) Q
  !f2py depend(nV,nk) V
  !f2py depend(ni,nj,nk) Q

  !Input
  integer, intent(in) ::  ni, nj, nk, nV
  double precision, intent(in) ::  V(nV,nk)

  !Output
  double precision, intent(out) ::  Q(ni,nj,nk)

  !Working
  integer i, j, k

  do k=1,nk
     do j=1,nj
        do i=1,ni
           Q(i,j,k) = V(ni*(j-1) + i, k)
        end do
     end do
  end do

end subroutine inflateVector
