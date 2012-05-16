subroutine getJacobian(nP, nJ, nT, nD, nsurf, nedge, ngroup, nvert, surf_vert, surf_edge, edge_group, group_k, group_m, group_n, group_d, surf_index_P, edge_index_P, surf_index_C, edge_index_C, edge_count, T, Ja, Ji, Jj)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nP, nJ, nT, nD, nsurf, nedge, ngroup, nvert, surf_vert, surf_edge, edge_group, group_k, group_m, group_n, group_d, surf_index_P, edge_index_P, surf_index_C, edge_index_C, edge_count, T
  !f2py intent(out) Ja, Ji, Jj
  !f2py depend(nsurf) surf_vert
  !f2py depend(nsurf) surf_edge
  !f2py depend(nedge) edge_group
  !f2py depend(ngroup) group_k
  !f2py depend(ngroup) group_m
  !f2py depend(ngroup) group_n
  !f2py depend(nD) group_d
  !f2py depend(nsurf) surf_index_P
  !f2py depend(nedge) edge_index_P
  !f2py depend(nsurf) surf_index_C
  !f2py depend(nedge) edge_index_C
  !f2py depend(nedge) edge_count
  !f2py depend(nT) T
  !f2py depend(nJ) Ja
  !f2py depend(nJ) Ji
  !f2py depend(nJ) Jj

  !Input
  integer, intent(in) ::  nP, nJ, nT, nD, nsurf, nedge, ngroup, nvert
  integer, intent(in) ::  surf_vert(nsurf,2,2), surf_edge(nsurf,2,2), edge_group(nedge), group_k(ngroup), group_m(ngroup), group_n(ngroup)
  double precision, intent(in) ::  group_d(nD)
  integer, intent(in) ::  surf_index_P(nsurf,2), edge_index_P(nedge,2), surf_index_C(nsurf,2), edge_index_C(nedge,2), edge_count(nedge)
  double precision, intent(in) ::  T(nT)

  !Output
  double precision, intent(out) ::  Ja(nJ)
  integer, intent(out) ::  Ji(nJ), Jj(nJ)

  !Working
  integer vert, surf
  integer iJ,u,v,k,k1,k2
  integer ugroup, vgroup
  integer ku,kv,mu,mv,nu,nv
  double precision, allocatable, dimension(:,:,:) ::  bufferT
  double precision, allocatable, dimension(:) ::  bufferD1, bufferD2
  double precision, allocatable, dimension(:) ::  Bu, Bv
  integer, allocatable, dimension(:,:) ::  mappingC, mappingP
  integer i0, j0

  Ja(:) = 0.0
  Ji(:) = 0
  Jj(:) = 0
  do vert=1,nvert
     Ja(vert) = 1.0
     Ji(vert) = vert
     Jj(vert) = vert
  end do

  iJ = nvert + 1
  do surf=1,nsurf
     ugroup = edge_group(abs(surf_edge(surf,1,1)))
     vgroup = edge_group(abs(surf_edge(surf,2,1)))
     ku = group_k(ugroup)
     kv = group_k(vgroup)
     mu = group_m(ugroup)
     mv = group_m(vgroup)
     nu = group_n(ugroup)
     nv = group_n(vgroup)
     allocate(bufferT(nu,nv,2))
     allocate(bufferD1(ku+mu))
     allocate(bufferD2(kv+mv))
     allocate(mappingC(mu,mv))
     allocate(mappingP(nu,nv))
     allocate(Bu(ku))
     allocate(Bv(kv))
     call extractSurfaceT(surf, nu, nv, nT, nsurf, nedge, surf_edge, surf_index_P, edge_index_P, T, bufferT)
     call extractD(ugroup, ngroup, nD, ku+mu, group_k, group_m, group_d, bufferD1)
     call extractD(vgroup, ngroup, nD, kv+mv, group_k, group_m, group_d, bufferD2)
     call getMapping(surf, mu, mv, nsurf, nedge, ngroup, nvert, surf_vert, surf_edge, edge_group, group_m, surf_index_C, edge_index_C, mappingC)
     call getMapping(surf, nu, nv, nsurf, nedge, ngroup, nvert, surf_vert, surf_edge, edge_group, group_n, surf_index_P, edge_index_P, mappingP)
     do v=1,2
        do u=2,nu-1
           call basis(ku,ku+mu,bufferT(u,1+(v-1)*(nv-1),1),bufferD1,Bu,i0)
           j0 = 1+(v-1)*(mv-1)
           do k=1,ku
              Ja(iJ) = Bu(k)/edge_count(abs(surf_edge(surf,1,v)))
              Ji(iJ) = mappingP(u,1+(v-1)*(nv-1))
              Jj(iJ) = mappingC(i0+k,j0)
              iJ = iJ + 1
           end do
        end do
     end do
     do u=1,2
        do v=2,nv-1
           call basis(kv,kv+mv,bufferT(1+(u-1)*(nu-1),v,2),bufferD2,Bv,j0)
           i0 = 1+(u-1)*(mu-1)
           do k=1,kv
              Ja(iJ) = Bv(k)/edge_count(abs(surf_edge(surf,2,u)))
              Ji(iJ) = mappingP(1+(u-1)*(nu-1),v)
              Jj(iJ) = mappingC(i0,j0+k)
              iJ = iJ + 1
           end do
        end do
     end do
     do u=2,nu-1
        do v=2,nv-1
           call basis(ku,ku+mu,bufferT(u,v,1),bufferD1,Bu,i0)
           call basis(kv,kv+mv,bufferT(u,v,2),bufferD2,Bv,j0)
           do k1=1,ku
              do k2=1,kv
                 Ja(iJ) = Bu(k1)*Bv(k2)
                 Ji(iJ) = mappingP(u,v)
                 Jj(iJ) = mappingC(i0+k1,j0+k2)
                 iJ = iJ + 1
              end do
           end do
        end do
     end do
     deallocate(mappingC)
     deallocate(mappingP)
     deallocate(Bu)
     deallocate(Bv)
     deallocate(bufferT)
     deallocate(bufferD1)
     deallocate(bufferD2)
  end do
  do iJ=1,nJ
     Ji(iJ) = Ji(iJ) - 1
     Jj(iJ) = Jj(iJ) - 1
  end do

end subroutine getJacobian



subroutine extractSurfaceT(surf, nu, nv, nT, nsurf, nedge, surf_edge, surf_index_P, edge_index_P, T, bufferT)

  implicit none

  !Input
  integer, intent(in) ::  surf, nu, nv, nT, nsurf, nedge
  integer, intent(in) ::  surf_edge(nsurf,2,2), surf_index_P(nsurf,2), edge_index_P(nedge,2)
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
     iT1 = edge_index_P(abs(surf_edge(surf,1,v)),1)
     iT2 = edge_index_P(abs(surf_edge(surf,1,v)),2)
     iT1 = iT1 + 1
     iT2 = iT2
     if (surf_edge(surf,1,v) .gt. 0) then
        bufferT(2:nu-1,1+(v-1)*(nv-1),1) = T(iT1:iT2)
     else
        bufferT(2:nu-1,1+(v-1)*(nv-1),1) = 1 - T(iT2:iT1:-1)
     end if
  end do
  do u=1,2
     iT1 = edge_index_P(abs(surf_edge(surf,2,u)),1)
     iT2 = edge_index_P(abs(surf_edge(surf,2,u)),2)
     iT1 = iT1 + 1
     iT2 = iT2
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

end subroutine extractSurfaceT



subroutine getMapping(surf, mu, mv, nsurf, nedge, ngroup, nvert, surf_vert, surf_edge, edge_group, group_m, surf_index_C, edge_index_C, mapping)

  implicit none

  !Input
  integer, intent(in) ::  surf, mu, mv, nsurf, nedge, ngroup, nvert
  integer, intent(in) ::  surf_vert(nsurf,2,2), surf_edge(nsurf,2,2), edge_group(nedge), group_m(ngroup), surf_index_C(nsurf,2), edge_index_C(nedge,2)

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



subroutine getJnnz(nsurf, nedge, ngroup, nvert, surf_edge, edge_group, group_k, group_n, vert_count, edge_count, nJ)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nsurf, nedge, ngroup, nvert, surf_edge, edge_group, group_k, group_n, vert_count, edge_count
  !f2py intent(out) nJ
  !f2py depend(nsurf) surf_edge
  !f2py depend(nedge) edge_group
  !f2py depend(ngroup) group_k
  !f2py depend(ngroup) group_n
  !f2py depend(nvert) vert_count
  !f2py depend(nedge) edge_count

  !Input
  integer, intent(in) ::  nsurf, nedge, ngroup, nvert
  integer, intent(in) ::  surf_edge(nsurf,2,2), edge_group(nedge), group_k(ngroup), group_n(ngroup), vert_count(nvert), edge_count(nedge)

  !Output
  integer, intent(out) ::  nJ

  !Working
  integer surf, edge
  integer ugroup, vgroup

  nJ = nvert
  do edge=1,nedge
     nJ = nJ + group_k(edge_group(edge))*(group_n(edge_group(edge)) - 2)*edge_count(edge)
  end do
  do surf=1,nsurf
     ugroup = edge_group(abs(surf_edge(surf,1,1)))
     vgroup = edge_group(abs(surf_edge(surf,2,1)))
     nJ = nJ + group_k(ugroup)*group_k(vgroup)*(group_n(ugroup)-2)*(group_n(vgroup)-2)
  end do

end subroutine getJnnz
