subroutine computeJacobian(nJ, nT, nD, nsurf, nedge, ngroup, nvert, surf_vert, & 
           surf_edge, edge_group, group_k, group_m, group_n, group_d, & 
           surf_index_P, edge_index_P, surf_index_C, edge_index_C, & 
           knot_index, edge_count, T, Ja, Ji, Jj)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nJ, nT, nD, nsurf, nedge, ngroup, nvert, surf_vert, surf_edge, edge_group, group_k, group_m, group_n, group_d, surf_index_P, edge_index_P, surf_index_C, edge_index_C, knot_index, edge_count, T
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
  !f2py depend(ngroup) knot_index
  !f2py depend(nedge) edge_count
  !f2py depend(nT) T
  !f2py depend(nJ) Ja
  !f2py depend(nJ) Ji
  !f2py depend(nJ) Jj

  !Input
  integer, intent(in) ::  nJ, nT, nD, nsurf, nedge, ngroup, nvert
  integer, intent(in) ::  surf_vert(nsurf,2,2), surf_edge(nsurf,2,2), & 
           edge_group(nedge), group_k(ngroup), group_m(ngroup), group_n(ngroup)
  double precision, intent(in) ::  group_d(nD)
  integer, intent(in) ::  surf_index_P(nsurf,2), edge_index_P(nedge,2), & 
           surf_index_C(nsurf,2), edge_index_C(nedge,2), &
           knot_index(ngroup,2), edge_count(nedge)
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
     call getSurfaceT(surf, nu, nv, nT, nsurf, nedge, & 
          surf_edge, surf_index_P, edge_index_P, T, bufferT)
     bufferD1 = group_d(knot_index(ugroup,1)+1:knot_index(ugroup,2))
     bufferD2 = group_d(knot_index(vgroup,1)+1:knot_index(vgroup,2))
     call getMapping(surf, mu, mv, nsurf, nedge, nvert, surf_vert, & 
          surf_edge, surf_index_C, edge_index_C, mappingC)
     call getMapping(surf, nu, nv, nsurf, nedge, nvert, surf_vert, & 
          surf_edge, surf_index_P, edge_index_P, mappingP)
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

end subroutine computeJacobian



subroutine computeJnnz(nsurf, nedge, ngroup, nvert, surf_edge, edge_group, & 
           group_k, group_n, edge_count, nJ)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nsurf, nedge, ngroup, nvert, surf_edge, edge_group, group_k, group_n, edge_count
  !f2py intent(out) nJ
  !f2py depend(nsurf) surf_edge
  !f2py depend(nedge) edge_group
  !f2py depend(ngroup) group_k
  !f2py depend(ngroup) group_n
  !f2py depend(nedge) edge_count

  !Input
  integer, intent(in) ::  nsurf, nedge, ngroup, nvert
  integer, intent(in) ::  surf_edge(nsurf,2,2), edge_group(nedge), & 
                          group_k(ngroup), group_n(ngroup), &
                          edge_count(nedge)

  !Output
  integer, intent(out) ::  nJ

  !Working
  integer surf, edge
  integer ugroup, vgroup

  nJ = nvert
  do edge=1,nedge
     nJ = nJ + group_k(edge_group(edge)) * (group_n(edge_group(edge)) - 2) * & 
          edge_count(edge)
  end do
  do surf=1,nsurf
     ugroup = edge_group(abs(surf_edge(surf,1,1)))
     vgroup = edge_group(abs(surf_edge(surf,2,1)))
     nJ = nJ + group_k(ugroup) * group_k(vgroup) * (group_n(ugroup)-2) * & 
         (group_n(vgroup)-2)
  end do

end subroutine computeJnnz
