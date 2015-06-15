subroutine computebnnz(nsurf, ngroup, &
     surf_group, order, num_pt, nB)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nsurf, ngroup, surf_group, order, num_pt
  !f2py intent(out) nB
  !f2py depend(nsurf) surf_group
  !f2py depend(ngroup) order, num_pt

  !Input
  integer, intent(in) ::  nsurf, ngroup
  integer, intent(in) ::  surf_group(nsurf, 2)
  integer, intent(in) ::  order(ngroup), num_pt(ngroup)

  !Output
  integer, intent(out) ::  nB

  !Working
  integer isurf, igroup_u, igroup_v

  nB = 0
  do isurf = 1, nsurf
     igroup_u = surf_group(isurf, 1)
     igroup_v = surf_group(isurf, 2)
     nB = nB + &
          order(igroup_u) * order(igroup_v) * &
          num_pt(igroup_u) * num_pt(igroup_v)
  end do

end subroutine computebnnz



subroutine computebmtx(nB, nsurf, ngroup, uder, vder, &
     str_indices_cp, str_indices_pt, &
     surf_group, order, num_cp, num_pt, &
     Ba, Bi, Bj)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nB, nsurf, ngroup, uder, vder, str_indices_cp, str_indices_pt, surf_group, order, num_cp, num_pt
  !f2py intent(out) Ba, Bi, Bj
  !f2py depend(nsurf) str_indices_cp, str_indices_pt, surf_group
  !f2py depend(ngroup) order, num_cp, num_pt
  !f2py depend(nB) Ba, Bi, Bj

  !Input
  integer, intent(in) ::  nB, nsurf, ngroup, uder, vder
  integer, intent(in) ::  str_indices_cp(nsurf, 2)
  integer, intent(in) ::  str_indices_pt(nsurf, 2)
  integer, intent(in) ::  surf_group(nsurf, 2)
  integer, intent(in) ::  order(ngroup)
  integer, intent(in) ::  num_cp(ngroup), num_pt(ngroup)

  !Output
  double precision, intent(out) ::  Ba(nB)
  integer, intent(out) ::  Bi(nB), Bj(nB)

  !Working
  integer iB, isurf, ind_u, ind_v, u, v, uu, vv
  integer ku, kv, mu, mv, nu, nv
  double precision, allocatable, dimension(:) ::  knot_u, knot_v
  double precision, allocatable, dimension(:) ::  param_u, param_v
  double precision, allocatable, dimension(:, :) ::  basis_u, basis_v
  integer, allocatable, dimension(:) ::  offset_u, offset_v

  iB = 0
  Ba(:) = 0.0
  Bi(:) = 0
  Bj(:) = 0

  do isurf = 1, nsurf
     ku = order(surf_group(isurf, 1))
     kv = order(surf_group(isurf, 2))
     mu = num_cp(surf_group(isurf, 1))
     mv = num_cp(surf_group(isurf, 2))
     nu = num_pt(surf_group(isurf, 1))
     nv = num_pt(surf_group(isurf, 2))

     allocate(knot_u(ku + mu))
     call knotopen(ku, ku+mu, knot_u)

     allocate(knot_v(kv + mv))
     call knotopen(kv, kv+mv, knot_v)

     allocate(param_u(nu))
     call paramuni(ku+mu, mu, nu, knot_u, param_u)

     allocate(param_v(nv))
     call paramuni(kv+mv, mv, nv, knot_v, param_v)

     allocate(basis_u(ku, nu))
     allocate(basis_v(kv, nv))
     allocate(offset_u(nu))
     allocate(offset_v(nv))

     if (uder .eq. 0) then
        do u = 1, nu
           call basis0(ku, ku+mu, param_u(u), knot_u, &
                basis_u(:, u), offset_u(u))
        end do
     else if (uder .eq. 1) then
        do u = 1, nu
           call basis1(ku, ku+mu, param_u(u), knot_u, &
                basis_u(:, u), offset_u(u))
        end do
     else if (uder .eq. 2) then
        do u = 1, nu
           call basis2(ku, ku+mu, param_u(u), knot_u, &
                basis_u(:, u), offset_u(u))
        end do
     end if

     if (vder .eq. 0) then
        do v = 1, nv
           call basis0(kv, kv+mv, param_v(v), knot_v, &
                basis_v(:, v), offset_v(v))
        end do
     else if (vder .eq. 1) then
        do v = 1, nv
           call basis1(kv, kv+mv, param_v(v), knot_v, &
                basis_v(:, v), offset_v(v))
        end do
     else if (vder .eq. 2) then
        do v = 1, nv
           call basis2(kv, kv+mv, param_v(v), knot_v, &
                basis_v(:, v), offset_v(v))
        end do
     end if
     
     do u = 1, nu
        do v = 1, nv
           do ind_u = 1, ku
              do ind_v = 1, kv
                 uu = offset_u(u) + ind_u
                 vv = offset_v(v) + ind_v
                 iB = iB + 1
                 Ba(iB) = basis_u(ind_u, u) * basis_v(ind_v, v)
                 Bi(iB) = str_indices_pt(isurf, 1) + (v-1) * nu + u
                 Bj(iB) = str_indices_cp(isurf, 1) + (vv-1) * mu + uu
              end do
           end do
        end do
     end do

     deallocate(knot_u)
     deallocate(knot_v)
     deallocate(param_u)
     deallocate(param_v)
     deallocate(basis_u)
     deallocate(basis_v)
     deallocate(offset_u)
     deallocate(offset_v)
  end do

  if (iB .ne. nB) then
     print *, 'Error in computeBmtx', iB, nB
     call exit(1)
  end if

end subroutine computebmtx
     
