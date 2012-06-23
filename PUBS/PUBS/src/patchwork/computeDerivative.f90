subroutine computeDerivativeSurface(surf, uder, vder, nB, nT, nD, nsurf, nedge,& 
           ngroup, nvert, surf_vert, surf_edge, edge_group, group_k, group_m, &
           group_n, group_d, surf_index_P, edge_index_P, surf_index_C, &
           edge_index_C, knot_index, T, Ba, Bi, Bj)

  implicit none
  
  !Fortran-python interface directives
  !f2py intent(in) surf, uder, vder, nB, nT, nD, nsurf, nedge, ngroup, nvert, surf_vert, surf_edge, edge_group, group_m, group_n, group_d, surf_index_P, edge_index_P, surf_index_C, edge_index_C, knot_index, T
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
  !f2py depend(ngroup) knot_index
  !f2py depend(nT) T
  !f2py depend(nB) Ba, Bi, Bj

  !Input
  integer, intent(in) ::  surf, uder, vder, nB, nT, nD, nsurf, nedge, ngroup, & 
                          nvert
  integer, intent(in) ::  surf_vert(nsurf,2,2), surf_edge(nsurf,2,2), &
           edge_group(nedge), group_k(ngroup), group_m(ngroup), group_n(ngroup)
  double precision, intent(in) ::  group_d(nD)
  integer, intent(in) ::  surf_index_P(nsurf,2), edge_index_P(nedge,2), & 
           surf_index_C(nsurf,2), edge_index_C(nedge,2)
  integer, intent(in) ::  knot_index(ngroup,2)
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
  call extractSurfaceT(surf, nu, nv, nT, nsurf, nedge, surf_edge, surf_index_P,&
       edge_index_P, T, bufferT)
  bufferD1 = group_d(knot_index(ugroup,1)+1:knot_index(ugroup,2))
  bufferD2 = group_d(knot_index(vgroup,1)+1:knot_index(vgroup,2))
  call getMapping(surf, mu, mv, nsurf, nedge, nvert, surf_vert, &
       surf_edge, surf_index_C, edge_index_C, mappingC)
  call getMapping(surf, nu, nv, nsurf, nedge, nvert, surf_vert, & 
       surf_edge, surf_index_P, edge_index_P, mappingP)
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



subroutine computeDerivativeSpecified(surf, uder, vder, nB, n1, n2, nD, nsurf, &
           nedge, ngroup, nvert, surf_vert, surf_edge, edge_group, group_k, & 
           group_m, group_n, group_d, surf_index_P, edge_index_P, surf_index_C,&
           edge_index_C, knot_index, bufferT, Ba, Bi, Bj)

  implicit none
  
  !Fortran-python interface directives
  !f2py intent(in) surf, uder, vder, nB, n1, n2, nD, nsurf, nedge, ngroup, nvert, surf_vert, surf_edge, edge_group, group_m, group_n, group_d, surf_index_C, edge_index_C, knot_index, bufferT
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
  !f2py depend(ngroup) knot_index
  !f2py depend(n1,n2) bufferT
  !f2py depend(nB) Ba, Bi, Bj

  !Input
  integer, intent(in) ::  surf, uder, vder, nB, n1, n2, nD, nsurf, nedge, & 
           ngroup, nvert
  integer, intent(in) ::  surf_vert(nsurf,2,2), surf_edge(nsurf,2,2), & 
           edge_group(nedge), group_k(ngroup), group_m(ngroup), group_n(ngroup)
  double precision, intent(in) ::  group_d(nD)
  integer, intent(in) ::  surf_index_P(nsurf,2), edge_index_P(nedge,2), & 
           surf_index_C(nsurf,2), edge_index_C(nedge,2)
  integer, intent(in) ::  knot_index(ngroup,2)
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
  bufferD1 = group_d(knot_index(ugroup,1)+1:knot_index(ugroup,2))
  bufferD2 = group_d(knot_index(vgroup,1)+1:knot_index(vgroup,2))
  call getMapping(surf, mu, mv, nsurf, nedge, nvert, surf_vert, & 
       surf_edge, surf_index_C, edge_index_C, mappingC)
  call getMapping(surf, nu, nv, nsurf, nedge, nvert, surf_vert, & 
       surf_edge, surf_index_P, edge_index_P, mappingP)
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



subroutine computeDerivativeSingle(surf, uder, vder, ku, kv, mu, mv, kpmu, &
           kpmv, nB, nD, nsurf, nedge, ngroup, nvert, t, surf_vert, &
           surf_edge, edge_group, group_d, surf_index_C, &
           edge_index_C, knot_index, Ba, Bi, Bj)

  implicit none
  
  !Fortran-python interface directives
  !f2py intent(in) surf, uder, vder, ku, kv, mu, mv, kpmu, kpmv, nB, nD, nsurf, nedge, ngroup, nvert, t, surf_vert, surf_edge, edge_group, group_d, surf_index_C, edge_index_C, knot_index
  !f2py intent(out) Ba, Bi, Bj
  !f2py depend(nsurf) surf_vert
  !f2py depend(nsurf) surf_edge
  !f2py depend(nedge) edge_group
  !f2py depend(nD) group_d
  !f2py depend(nsurf) surf_index_C
  !f2py depend(nedge) edge_index_C
  !f2py depend(ngroup) knot_index
  !f2py depend(nB) Ba, Bi, Bj

  !Input
  integer, intent(in) ::  surf, uder, vder, ku, kv, mu, mv, kpmu, kpmv, nB, &
                          nD, nsurf, nedge, ngroup, nvert
  double precision, intent(in) ::  t(2)
  integer, intent(in) ::  surf_vert(nsurf,2,2), surf_edge(nsurf,2,2), & 
                          edge_group(nedge)
  double precision, intent(in) ::  group_d(nD)
  integer, intent(in) ::  surf_index_C(nsurf,2), edge_index_C(nedge,2), &
                          knot_index(ngroup,2)

  !Output
  double precision, intent(out) ::  Ba(nB)
  integer, intent(out) ::  Bi(nB), Bj(nB)

  !Working
  integer iB
  integer k1, k2
  integer ugroup, vgroup
  double precision du(kpmu), dv(kpmv)
  double precision Bu(ku), Bv(kv)
  integer i0, j0
  integer index

  ugroup = edge_group(abs(surf_edge(surf,1,1)))
  vgroup = edge_group(abs(surf_edge(surf,2,1)))

  du = group_d(knot_index(ugroup,1)+1:knot_index(ugroup,2))
  dv = group_d(knot_index(vgroup,1)+1:knot_index(vgroup,2))

  if (uder.eq.0) then
     call basis(ku,kpmu,t(1),du,Bu,i0)
  else if (uder.eq.1) then
     call basis1(ku,kpmu,t(1),du,Bu,i0)
  else if (uder.eq.2) then
     call basis2(ku,kpmu,t(1),du,Bu,i0)
  end if
  if (vder.eq.0) then
     call basis(kv,kpmv,t(2),dv,Bv,j0)
  else if (vder.eq.1) then
     call basis1(kv,kpmv,t(2),dv,Bv,j0)
  else if (vder.eq.2) then
     call basis2(kv,kpmv,t(2),dv,Bv,j0)
  end if
  iB = 1
  do k1=1,ku
     do k2=1,kv
        call computeIndex(surf, i0+k1, j0+k2, mu, mv, nsurf, nedge, nvert, &
             surf_vert, surf_edge, surf_index_C, edge_index_C, index)
        Ba(iB) = Bu(k1)*Bv(k2)
        Bi(iB) = 1
        Bj(iB) = index
        iB = iB + 1
     end do
  end do
  do iB=1,nB
     Bi(iB) = Bi(iB) - 1
     Bj(iB) = Bj(iB) - 1
  end do

end subroutine computeDerivativeSingle



subroutine computeBnnz(nP, nsurf, nedge, ngroup, surf_edge, edge_group, &
     group_k, s, nB)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nP, nsurf, nedge, ngroup, surf_edge, edge_group, group_k, s
  !f2py intent(out) nB
  !f2py depend(nsurf) surf_edge
  !f2py depend(nedge) edge_group
  !f2py depend(ngroup) group_k
  !f2py depend(nP) s

  !Input
  integer, intent(in) ::  nP, nsurf, nedge, ngroup
  integer, intent(in) ::  surf_edge(nsurf,2,2), edge_group(nedge), &
       group_k(ngroup), s(nP)
  
  !Output
  integer, intent(out) ::  nB

  !Working
  integer iP, surf, ugroup, vgroup, ku, kv

  nB = 0
  do iP=1,nP
     surf = s(iP)
     ugroup = edge_group(abs(surf_edge(surf,1,1)))
     vgroup = edge_group(abs(surf_edge(surf,2,1)))
     ku = group_k(ugroup)
     kv = group_k(vgroup)
     nB = nB + ku*kv
  end do

end subroutine computeBnnz



subroutine computeBases(uder, vder, nP, nB, nD, nsurf, nedge, ngroup, nvert,&
     surf_vert, surf_edge, edge_group, group_k, group_m, group_d, &
     surf_index_C, edge_index_C, knot_index, s, u, v, Ba, Bi, Bj)

  implicit none
  
  !Fortran-python interface directives
  !f2py intent(in) uder, vder, nP, nB, nD, nsurf, nedge, ngroup, nvert, surf_vert, surf_edge, edge_group, group_k, group_m, group_d, surf_index_C, edge_index_C, knot_index, s, u, v
  !f2py intent(out) Ba, Bi, Bj
  !f2py depend(nsurf) surf_vert
  !f2py depend(nsurf) surf_edge
  !f2py depend(nedge) edge_group
  !f2py depend(ngroup) group_k, group_m
  !f2py depend(nD) group_d
  !f2py depend(nsurf) surf_index_C
  !f2py depend(nedge) edge_index_C
  !f2py depend(ngroup) knot_index
  !f2py depend(nP) s, u, v
  !f2py depend(nB) Ba, Bi, Bj

  !Input
  integer, intent(in) ::  uder, vder, nP, nB, nD, nsurf, nedge, ngroup, nvert
  integer, intent(in) ::  surf_vert(nsurf,2,2), surf_edge(nsurf,2,2), & 
                          edge_group(nedge), group_k(ngroup), group_m(ngroup)
  double precision, intent(in) ::  group_d(nD)
  integer, intent(in) ::  surf_index_C(nsurf,2), edge_index_C(nedge,2), &
                          knot_index(ngroup,2)
  integer, intent(in) ::  s(nP)
  double precision, intent(in) ::  u(nP), v(nP)

  !Output
  double precision, intent(out) ::  Ba(nB)
  integer, intent(out) ::  Bi(nB), Bj(nB)

  !Working
  integer iB, iP, k1, k2
  integer surf, ugroup, vgroup, ku, kv, mu, mv
  integer i0, j0
  integer index
  double precision, allocatable, dimension(:) ::  du, dv, Bu, Bv

  iB = 1
  do iP=1,nP
     surf = s(iP)
     ugroup = edge_group(abs(surf_edge(surf,1,1)))
     vgroup = edge_group(abs(surf_edge(surf,2,1)))
     ku = group_k(ugroup)
     kv = group_k(vgroup)
     mu = group_m(ugroup)
     mv = group_m(vgroup)

     allocate(du(ku+mu))
     allocate(dv(kv+mv))
     allocate(Bu(ku))
     allocate(Bv(kv))

     du = group_d(knot_index(ugroup,1)+1:knot_index(ugroup,2))
     dv = group_d(knot_index(vgroup,1)+1:knot_index(vgroup,2))

     if (uder.eq.0) then
        call basis(ku,ku+mu,u(iP),du,Bu,i0)
     else if (uder.eq.1) then
        call basis1(ku,ku+mu,u(iP),du,Bu,i0)
     else if (uder.eq.2) then
        call basis2(ku,ku+mu,u(iP),du,Bu,i0)
     end if
     if (vder.eq.0) then
        call basis(kv,kv+mv,v(iP),dv,Bv,j0)
     else if (vder.eq.1) then
        call basis1(kv,kv+mv,v(iP),dv,Bv,j0)
     else if (vder.eq.2) then
        call basis2(kv,kv+mv,v(iP),dv,Bv,j0)
     end if

     do k1=1,ku
        do k2=1,kv
           call computeIndex(surf, i0+k1, j0+k2, mu, mv, nsurf, nedge, nvert, &
                surf_vert, surf_edge, surf_index_C, edge_index_C, index)
           Ba(iB) = Bu(k1)*Bv(k2)
           Bi(iB) = iP
           Bj(iB) = index
           iB = iB + 1
        end do
     end do
     deallocate(du)
     deallocate(dv)
     deallocate(Bu)
     deallocate(Bv)
  end do        

  do iB=1,nB
     Bi(iB) = Bi(iB) - 1
     Bj(iB) = Bj(iB) - 1
  end do 

end subroutine computeBases
