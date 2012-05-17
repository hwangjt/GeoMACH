subroutine extractD(group, ngroup, nD, kpm, group_k, group_m, group_d, d)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) group, ngroup, nD, kpm, group_k, group_m, group_d
  !f2py intent(out) d
  !f2py depend(ngroup) group_k, group_m
  !f2py depend(nD) group_d
  !f2py depend(kpm) d

  !Input
  integer, intent(in) ::  group, ngroup, nD, kpm
  integer, intent(in) ::  group_k(ngroup), group_m(ngroup)
  double precision, intent(in) ::  group_d(nD)

  !Output
  double precision, intent(out) ::  d(kpm)

  !Working
  integer group0, d1, d2

  d2 = 0
  do group0=1,group
     d2 = d2 + group_k(group0) + group_m(group0)
  end do
  d1 = d2 - kpm + 1
  d(:) = group_d(d1:d2)
  
end subroutine extractD



subroutine initializeT(nT, nD, nsurf, nedge, ngroup, surf_edge, edge_group, & 
           group_k, group_m, group_n, group_d, T)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nT, nD, nsurf, nedge, ngroup, surf_edge
  !f2py intent(in) edge_group, group_k, group_m, group_n, group_d
  !f2py intent(out) T
  !f2py depend(nsurf) surf_edge
  !f2py depend(nedge) edge_group
  !f2py depend(ngroup) group_k,group_m,group_n
  !f2py depend(nD) group_d
  !f2py depend(nT) T

  !Input
  integer, intent(in) ::  nT, nD, nsurf, nedge, ngroup
  integer, intent(in) ::  surf_edge(nsurf,2,2), edge_group(nedge), & 
                          group_k(ngroup), group_m(ngroup), group_n(ngroup)
  double precision, intent(in) ::  group_d(nD)

  !Output
  double precision, intent(out) ::  T(nT)

  !Working
  double precision, allocatable, dimension(:) ::  bufferD, bufferT
  integer edge, surf
  integer u, v
  integer group, ugroup, vgroup
  integer k, m, n, ku, kv, mu, mv, nu, nv 
  integer iT1, iT2, d1, d2

  iT1 = 1
  iT2 = 1
  do edge=1,nedge
     group = edge_group(edge)
     k = group_k(group)
     m = group_m(group)
     n = group_n(group)
     allocate(bufferD(k+m))
     allocate(bufferT(n))
     call extractD(group, ngroup, nD, k+m, group_k, group_m, group_d, bufferD)
     call paramuni(k+m,m,n,bufferD,bufferT)
     iT2 = iT1 + n - 2 - 1
     T(iT1:iT2) = bufferT(2:n-1)
     iT1 = iT1 + n - 2
     deallocate(bufferD)
     deallocate(bufferT)
  end do
  do surf=1,nsurf
     ugroup = edge_group(abs(surf_edge(surf,1,1)))
     vgroup = edge_group(abs(surf_edge(surf,2,1)))
     ku = group_k(ugroup)
     kv = group_k(vgroup)
     mu = group_m(ugroup)
     mv = group_m(vgroup)
     nu = group_n(ugroup)
     nv = group_n(vgroup)
     allocate(bufferD(ku+mu))
     allocate(bufferT(nu))
     call extractD(ugroup, ngroup, nD, ku+mu, group_k, group_m, & 
          group_d, bufferD)
     call paramuni(ku+mu,mu,nu,bufferD,bufferT)
     do v=2,nv-1
        do u=2,nu-1
           T(iT1) = bufferT(u)
           iT1 = iT1 + 1
        end do
     end do
     deallocate(bufferD)
     deallocate(bufferT)
     allocate(bufferD(kv+mv))
     allocate(bufferT(nv))
     call extractD(vgroup, ngroup, nD, kv+mv, group_k, group_m, & 
          group_d, bufferD)
     call paramuni(kv+mv,mv,nv,bufferD,bufferT)
     do v=2,nv-1
        do u=2,nu-1
           T(iT1) = bufferT(v)
           iT1 = iT1 + 1
        end do
     end do
     deallocate(bufferD)
     deallocate(bufferT)
  end do

end subroutine initializeT
