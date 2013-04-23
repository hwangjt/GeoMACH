subroutine initializeKMN(k, nsurf, nedge, ngroup, ratio, ns, surf_edge, edge_group, &
           group_k, group_m, group_n)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) k, nsurf, nedge, ngroup, ratio, ns, surf_edge, edge_group
  !f2py intent(out) group_k, group_m, group_n
  !f2py depend(nsurf) ns
  !f2py depend(nsurf) surf_edge
  !f2py depend(nedge) edge_group
  !f2py depend(ngroup) group_k
  !f2py depend(ngroup) group_m
  !f2py depend(ngroup) group_n

  !Input
  integer, intent(in) ::  k, nsurf, nedge, ngroup
  double precision, intent(in) ::  ratio
  integer, intent(in) ::  ns(nsurf,2), surf_edge(nsurf,2,2)
  integer, intent(in) ::  edge_group(nedge)

  !Output
  integer, intent(out) ::  group_k(ngroup)
  integer, intent(out) ::  group_m(ngroup)
  integer, intent(out) ::  group_n(ngroup)

  !Working
  integer surf
  integer m, ugroup, vgroup

  group_k(:) = k
  group_m(:) = k
  group_n(:) = k
  do surf=1,nsurf
     ugroup = edge_group(abs(surf_edge(surf,1,1)))
     vgroup = edge_group(abs(surf_edge(surf,2,1)))
     group_n(ugroup) = ns(surf,1)
     group_n(vgroup) = ns(surf,2)
     m = nint(ns(surf,1)/ratio)
     if (m .lt. k) then
        m = k
     end if
     group_m(ugroup) = m
     m = nint(ns(surf,2)/ratio)
     if (m .lt. k) then
        m = k
     end if
     group_m(vgroup) = m
  end do

end subroutine initializeKMN
