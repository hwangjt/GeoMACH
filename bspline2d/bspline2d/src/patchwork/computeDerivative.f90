subroutine computeDerivativeSurface(surf, uder, vder, nB, nT, nD, nsurf, nedge, ngroup, nvert, surf_vert, surf_edge, edge_group, group_k, group_m, group_n, group_d, surf_index_P, edge_index_P, surf_index_C, edge_index_C, T, Ba, Bi, Bj)

  implicit none
  
  !Fortran-python interface directives
  !f2py intent(in) surf, uder, vder, nB, nT, nD, nsurf, nedge, ngroup, nvert, surf_vert, surf_edge, edge_group, group_m, group_n, group_d, surf_index_C, edge_index_C
  !f2py intent(out) Ba, Bi, Bj
  !f2py depend(nsurf) surf_vert
  !f2py depend(nsurf) surf_edge
  !f2py depend(nedge) edge_group
  !f2py depend(ngroup) group_k, group_m, group_n
  !f2py depend(nD) group_d
  !f2py depend(nsurf) surf_index_P
  !f2py depend(nedge) edge_index_P
  !f2py depend(nsurf) surf_index_C
  !f2py depend(nedge) edge_index_C
  !f2py depend(nT) T
  !f2py depend(nB) Ba, Bi, Bj

  !Input
  integer, intent(in) ::  surf, uder, vder, nB, nT, nD, nsurf, nedge, ngroup, nvert
  integer, intent(in) ::  surf_vert(nsurf,2,2), surf_edge(nsurf,2,2), edge_group(nedge), group_k(ngroup), group_m(ngroup), group_n(ngroup)
  double precision, intent(in) ::  group_d(nD)
  integer, intent(in) ::  surf_index_P(nsurf,2), edge_index_P(nedge,2), surf_index_C(nsurf,2), edge_index_C(nedge,2)
  double precision, intent(in) ::  T(nT)

  !Output
  double precision, intent(out) ::  Ba(nB)
  integer, intent(out) ::  Bi(nB), Bj(nB)

  !Working
  integer iB
  integer u, v
  integer k1, k2
  integer ugroup, vgroup, ku, kv, mu, mv, nu, nv
  double precision, allocatable, dimension(:,:,:) ::  bufferT
  double precision, allocatable, dimension(:) ::  bufferD1, bufferD2
  double precision, allocatable, dimension(:) ::  Bu, Bv
  integer, allocatable, dimension(:,:) ::  mappingC, mappingP
  integer i0, j0

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
  iB = 1
  do u=1,nu
     do v=1,nv
        if (uder.eq.0) then
           call basis(ku,ku+mu,bufferT(u,v,1),bufferD1,Bu,i0)
        else if (uder.eq.1) then
           call basis1(ku,ku+mu,bufferT(u,v,1),bufferD1,Bu,i0)
        else if (uder.eq.2) then
           call basis2(ku,ku+mu,bufferT(u,v,1),bufferD1,Bu,i0)
        end if
        if (vder.eq.0) then
           call basis(kv,kv+mv,bufferT(u,v,2),bufferD2,Bv,j0)
        else if (vder.eq.1) then
           call basis1(kv,kv+mv,bufferT(u,v,2),bufferD2,Bv,j0)
        else if (vder.eq.2) then
           call basis2(kv,kv+mv,bufferT(u,v,2),bufferD2,Bv,j0)
        end if
        do k1=1,ku
           do k2=1,kv
              Ba(iB) = Bu(k1)*Bv(k2)
              Bi(iB) = mappingP(u,v)
              Bj(iB) = mappingC(i0+k1,j0+k2)
              iB = iB + 1
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
  do iB=1,nB
     Bi(iB) = Bi(iB) - 1
     Bj(iB) = Bj(iB) - 1
  end do

end subroutine computeDerivativeSurface



subroutine computeDerivativeSpecified(surf, uder, vder, nB, n1, n2, nD, nsurf, nedge, ngroup, nvert, surf_vert, surf_edge, edge_group, group_k, group_m, group_n, group_d, surf_index_P, edge_index_P, surf_index_C, edge_index_C, bufferT, Ba, Bi, Bj)

  implicit none
  
  !Fortran-python interface directives
  !f2py intent(in) surf, uder, vder, nB, n1, n2, nD, nsurf, nedge, ngroup, nvert, surf_vert, surf_edge, edge_group, group_m, group_n, group_d, surf_index_C, edge_index_C, bufferT
  !f2py intent(out) Ba, Bi, Bj
  !f2py depend(nsurf) surf_vert
  !f2py depend(nsurf) surf_edge
  !f2py depend(nedge) edge_group
  !f2py depend(ngroup) group_k, group_m, group_n
  !f2py depend(nD) group_d
  !f2py depend(nsurf) surf_index_P
  !f2py depend(nedge) edge_index_P
  !f2py depend(nsurf) surf_index_C
  !f2py depend(nedge) edge_index_C
  !f2py depend(n1,n2) bufferT
  !f2py depend(nB) Ba, Bi, Bj

  !Input
  integer, intent(in) ::  surf, uder, vder, nB, n1, n2, nD, nsurf, nedge, ngroup, nvert
  integer, intent(in) ::  surf_vert(nsurf,2,2), surf_edge(nsurf,2,2), edge_group(nedge), group_k(ngroup), group_m(ngroup), group_n(ngroup)
  double precision, intent(in) ::  group_d(nD)
  integer, intent(in) ::  surf_index_P(nsurf,2), edge_index_P(nedge,2), surf_index_C(nsurf,2), edge_index_C(nedge,2)
  double precision, intent(in) ::  bufferT(n1,n2,2)

  !Output
  double precision, intent(out) ::  Ba(nB)
  integer, intent(out) ::  Bi(nB), Bj(nB)

  !Working
  integer iB
  integer u, v
  integer k1, k2
  integer ugroup, vgroup, ku, kv, mu, mv, nu, nv
  double precision, allocatable, dimension(:) ::  bufferD1, bufferD2
  double precision, allocatable, dimension(:) ::  Bu, Bv
  integer, allocatable, dimension(:,:) ::  mappingC, mappingP
  integer i0, j0

  ugroup = edge_group(abs(surf_edge(surf,1,1)))
  vgroup = edge_group(abs(surf_edge(surf,2,1)))
  ku = group_k(ugroup)
  kv = group_k(vgroup)
  mu = group_m(ugroup)
  mv = group_m(vgroup)
  nu = group_n(ugroup)
  nv = group_n(vgroup)
  allocate(bufferD1(ku+mu))
  allocate(bufferD2(kv+mv))
  allocate(mappingC(mu,mv))
  allocate(mappingP(nu,nv))
  allocate(Bu(ku))
  allocate(Bv(kv))
  call extractD(ugroup, ngroup, nD, ku+mu, group_k, group_m, group_d, bufferD1)
  call extractD(vgroup, ngroup, nD, kv+mv, group_k, group_m, group_d, bufferD2)
  call getMapping(surf, mu, mv, nsurf, nedge, ngroup, nvert, surf_vert, surf_edge, edge_group, group_m, surf_index_C, edge_index_C, mappingC)
  call getMapping(surf, nu, nv, nsurf, nedge, ngroup, nvert, surf_vert, surf_edge, edge_group, group_n, surf_index_P, edge_index_P, mappingP)
  iB = 1
  do u=1,n1
     do v=1,n2
        if (uder.eq.0) then
           call basis(ku,ku+mu,bufferT(u,v,1),bufferD1,Bu,i0)
        else if (uder.eq.1) then
           call basis1(ku,ku+mu,bufferT(u,v,1),bufferD1,Bu,i0)
        else if (uder.eq.2) then
           call basis2(ku,ku+mu,bufferT(u,v,1),bufferD1,Bu,i0)
        end if
        if (vder.eq.0) then
           call basis(kv,kv+mv,bufferT(u,v,2),bufferD2,Bv,j0)
        else if (vder.eq.1) then
           call basis1(kv,kv+mv,bufferT(u,v,2),bufferD2,Bv,j0)
        else if (vder.eq.2) then
           call basis2(kv,kv+mv,bufferT(u,v,2),bufferD2,Bv,j0)
        end if
        do k1=1,ku
           do k2=1,kv
              Ba(iB) = Bu(k1)*Bv(k2)
              Bi(iB) = u + (v-1)*n1 !mappingP(u,v)
              Bj(iB) = mappingC(i0+k1,j0+k2)
              iB = iB + 1
           end do
        end do
     end do
  end do
  deallocate(mappingC)
  deallocate(mappingP)
  deallocate(Bu)
  deallocate(Bv)
  deallocate(bufferD1)
  deallocate(bufferD2)
  do iB=1,nB
     Bi(iB) = Bi(iB) - 1
     Bj(iB) = Bj(iB) - 1
  end do

end subroutine computeDerivativeSpecified



subroutine computeDerivativeSingle(surf, uder, vder, nB, nD, nsurf, nedge, ngroup, nvert, t, surf_vert, surf_edge, edge_group, group_k, group_m, group_n, group_d, surf_index_P, edge_index_P, surf_index_C, edge_index_C, Ba, Bi, Bj)

  implicit none
  
  !Fortran-python interface directives
  !f2py intent(in) surf, uder, vder, nB, nD, nsurf, nedge, ngroup, nvert, t, surf_vert, surf_edge, edge_group, group_k, group_m, group_n, group_d, surf_index_P, edge_index_P, surf_index_C, edge_index_C
  !f2py intent(out) Ba, Bi, Bj
  !f2py depend(nsurf) surf_vert
  !f2py depend(nsurf) surf_edge
  !f2py depend(nedge) edge_group
  !f2py depend(ngroup) group_k, group_m, group_n
  !f2py depend(nD) group_d
  !f2py depend(nsurf) surf_index_P
  !f2py depend(nedge) edge_index_P
  !f2py depend(nsurf) surf_index_C
  !f2py depend(nedge) edge_index_C
  !f2py depend(nB) Ba, Bi, Bj

  !Input
  integer, intent(in) ::  surf, uder, vder, nB, nD, nsurf, nedge, ngroup, nvert
  double precision, intent(in) ::  t(2)
  integer, intent(in) ::  surf_vert(nsurf,2,2), surf_edge(nsurf,2,2), edge_group(nedge), group_k(ngroup), group_m(ngroup), group_n(ngroup)
  double precision, intent(in) ::  group_d(nD)
  integer, intent(in) ::  surf_index_P(nsurf,2), edge_index_P(nedge,2), surf_index_C(nsurf,2), edge_index_C(nedge,2)

  !Output
  double precision, intent(out) ::  Ba(nB)
  integer, intent(out) ::  Bi(nB), Bj(nB)

  !Working
  integer iB
  integer u, v
  integer k1, k2
  integer ugroup, vgroup, ku, kv, mu, mv, nu, nv
  double precision, allocatable, dimension(:) ::  bufferD1, bufferD2
  double precision, allocatable, dimension(:) ::  Bu, Bv
  integer, allocatable, dimension(:,:) ::  mappingC, mappingP
  integer i0, j0

  ugroup = edge_group(abs(surf_edge(surf,1,1)))
  vgroup = edge_group(abs(surf_edge(surf,2,1)))
  ku = group_k(ugroup)
  kv = group_k(vgroup)
  mu = group_m(ugroup)
  mv = group_m(vgroup)
  nu = group_n(ugroup)
  nv = group_n(vgroup)
  allocate(bufferD1(ku+mu))
  allocate(bufferD2(kv+mv))
  allocate(mappingC(mu,mv))
  allocate(mappingP(nu,nv))
  allocate(Bu(ku))
  allocate(Bv(kv))
  call extractD(ugroup, ngroup, nD, ku+mu, group_k, group_m, group_d, bufferD1)
  call extractD(vgroup, ngroup, nD, kv+mv, group_k, group_m, group_d, bufferD2)
  call getMapping(surf, mu, mv, nsurf, nedge, ngroup, nvert, surf_vert, surf_edge, edge_group, group_m, surf_index_C, edge_index_C, mappingC)
  call getMapping(surf, nu, nv, nsurf, nedge, ngroup, nvert, surf_vert, surf_edge, edge_group, group_n, surf_index_P, edge_index_P, mappingP)
  iB = 1
  if (uder.eq.0) then
     call basis(ku,ku+mu,t(1),bufferD1,Bu,i0)
  else if (uder.eq.1) then
     call basis1(ku,ku+mu,t(1),bufferD1,Bu,i0)
  else if (uder.eq.2) then
     call basis2(ku,ku+mu,t(1),bufferD1,Bu,i0)
  end if
  if (vder.eq.0) then
     call basis(kv,kv+mv,t(2),bufferD2,Bv,j0)
  else if (vder.eq.1) then
     call basis1(kv,kv+mv,t(2),bufferD2,Bv,j0)
  else if (vder.eq.2) then
     call basis2(kv,kv+mv,t(2),bufferD2,Bv,j0)
  end if
  do k1=1,ku
     do k2=1,kv
        Ba(iB) = Bu(k1)*Bv(k2)
        Bi(iB) = 1
        Bj(iB) = mappingC(i0+k1,j0+k2)
        iB = iB + 1
     end do
  end do
  deallocate(mappingC)
  deallocate(mappingP)
  deallocate(Bu)
  deallocate(Bv)
  deallocate(bufferD1)
  deallocate(bufferD2)
  do iB=1,nB
     Bi(iB) = Bi(iB) - 1
     Bj(iB) = Bj(iB) - 1
  end do

end subroutine computeDerivativeSingle
