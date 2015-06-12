subroutine computeproj(npt_proj, nsurf_proj, &
     ncp_str, npt_str, nsurf, ngroup, &
     surf_proj, str_indices_cp, str_indices_pt, &
     surf_group, order, num_cp, num_pt, &
     cp_str, pt_str, pt_proj, &
     min_s, min_tu, min_tv)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) npt_proj, nsurf_proj, ncp_str, npt_str, nsurf, ngroup, surf_proj, str_indices_cp, str_indices_pt, surf_group, order, num_cp, num_pt, cp_str, pt_str, pt_proj
  !f2py intent(out) min_s, min_tu, min_tv
  !f2py depend(nsurf_proj) surf_proj
  !f2py depend(nsurf) str_indices_cp, str_indices_pt, surf_group
  !f2py depend(ngroup) order, num_cp, num_pt
  !f2py depend(ncp_str) cp_str
  !f2py depend(npt_str) pt_str
  !f2py depend(npt_proj) pt_proj
  !f2py depend(npt_proj) min_s, min_tu, min_tv

  !Input
  integer, intent(in) ::  npt_proj, nsurf_proj
  integer, intent(in) ::  ncp_str, npt_str, nsurf, ngroup
  integer, intent(in) ::  surf_proj(nsurf_proj)
  integer, intent(in) ::  str_indices_cp(nsurf, 2)
  integer, intent(in) ::  str_indices_pt(nsurf, 2)
  integer, intent(in) ::  surf_group(nsurf, 2), order(ngroup)
  integer, intent(in) ::  num_cp(ngroup), num_pt(ngroup)
  double precision, intent(in) ::  cp_str(ncp_str, 3)
  double precision, intent(in) ::  pt_str(npt_str, 3)
  double precision, intent(in) ::  pt_proj(npt_proj, 3)

  !Output
  integer, intent(out) ::  min_s(npt_proj)
  double precision, intent(out) ::  min_tu(npt_proj)
  double precision, intent(out) ::  min_tv(npt_proj)

  !Working
  double precision min_d(npt_proj)
  integer isurf_proj, isurf, ipt_proj
  integer offset_cp, offset_pt, offset_u, offset_v
  integer ku, kv, mu, mv, nu, nv
  integer u, v, min_u, min_v
  integer lu, lv, hu, hv, k, counter
  double precision pt0(3), pt(3), min_d_surf, d
  double precision P(3), Pu(3), Pv(3), Puu(3), Puv(3), Pvv(3)
  double precision f(3), x(2), g(2), H(2, 2), det
  double precision norm_g, norm_dx, W(2, 2), dx(2)
  double precision, allocatable, dimension(:) ::  knot_u, knot_v
  double precision, allocatable, dimension(:) ::  param_u, param_v
  double precision, allocatable, dimension(:) ::  bu0, bu1, bu2
  double precision, allocatable, dimension(:) ::  bv0, bv1, bv2
  double precision, allocatable, dimension(:, :, :) ::  cp

  min_d(:) = 1e10
  min_s(:) = 1
  min_tu(:) = 1.0
  min_tv(:) = 1.0

  do isurf_proj = 1, nsurf_proj
     isurf = surf_proj(isurf_proj)
     offset_cp = str_indices_cp(isurf, 1)
     offset_pt = str_indices_pt(isurf, 1)

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

     allocate(bu0(ku))
     allocate(bu1(ku))
     allocate(bu2(ku))
     allocate(bv0(kv))
     allocate(bv1(kv))
     allocate(bv2(kv))

     allocate(cp(mu, mv, 3))
     do u = 1, mu
        do v = 1, mv
           cp(u, v, :) = &
                cp_str(offset_cp + (v-1)*mu + u, :)
        end do
     end do

     do ipt_proj = 1, npt_proj
        pt0 = pt_proj(ipt_proj, :)

        min_d_surf = 1e10
        min_u = 1
        min_v = 1
        do u = 1, nu ! ceiling(nu/100.0)
           do v = 1, nv ! ceiling(nv/100.0)
              pt = pt_str(offset_pt + (v-1)*nu + u, :)
              d = abs(dot_product(pt - pt0, pt - pt0))
              if (d .lt. min_d_surf) then
                 min_d_surf = d
                 min_u = u
                 min_v = v
              end if
           end do
        end do
        if (min_d_surf .lt. min_d(ipt_proj)) then
           min_d(ipt_proj) = min_d_surf
           min_s(ipt_proj) = isurf
           min_tu(ipt_proj) = param_u(min_u)
           min_tv(ipt_proj) = param_v(min_v)
        end if

        x(1) = param_u(min_u)
        x(2) = param_v(min_v)
        do counter = 0, 40
           call basis0(ku, ku+mu, x(1), knot_u, bu0, offset_u)
           call basis1(ku, ku+mu, x(1), knot_u, bu1, offset_u)
           call basis2(ku, ku+mu, x(1), knot_u, bu2, offset_u)
           call basis0(kv, kv+mv, x(2), knot_v, bv0, offset_v)
           call basis1(kv, kv+mv, x(2), knot_v, bv1, offset_v)
           call basis2(kv, kv+mv, x(2), knot_v, bv2, offset_v)
           P(:) = 0.0
           Pu(:) = 0.0
           Pv(:) = 0.0
           Puu(:) = 0.0
           Puv(:) = 0.0
           Pvv(:) = 0.0
           do lu = 1, ku
              hu = lu + offset_u
              do lv = 1, kv
                 hv = lv + offset_v
                 P   = P   + bu0(lu) * bv0(lv) * cp(hu, hv, :)
                 Pu  = Pu  + bu1(lu) * bv0(lv) * cp(hu, hv, :)
                 Pv  = Pv  + bu0(lu) * bv1(lv) * cp(hu, hv, :)
                 Puu = Puu + bu2(lu) * bv0(lv) * cp(hu, hv, :)
                 Puv = Puv + bu1(lu) * bv1(lv) * cp(hu, hv, :)
                 Pvv = Pvv + bu0(lu) * bv2(lv) * cp(hu, hv, :)
              end do
           end do

           f = P - pt0
           g(1) = 2 * dot_product(f, Pu)
           g(2) = 2 * dot_product(f, Pv)
           H(1, 1) = 2 * dot_product(Pu, Pu) + 2 * dot_product(f, Puu)
           H(1, 2) = 2 * dot_product(Pu, Pv) + 2 * dot_product(f, Puv)
           H(2, 2) = 2 * dot_product(Pv, Pv) + 2 * dot_product(f, Pvv)
           H(2, 1) = H(1, 2)
           do k = 1, 2
              if (((x(k) .eq. 0) .and. (g(k) .gt. 0)) .or. &
                   ((x(k) .eq. 1) .and. (g(k) .lt. 0))) then
                 g(k) = 0.0
                 H(1, 2) = 0.0
                 H(2, 1) = 0.0
                 H(k, k) = 1.0
              end if
           end do
           det = H(1, 1) * H(2, 2) - H(1, 2) * H(2, 1)
           W(1, 1) = H(2, 2) / det
           W(1, 2) = -H(1, 2) / det
           W(2, 2) = H(1, 1) / det
           W(2, 1) = W(1, 2)
           norm_g = sqrt(dot_product(g, g))
           dx(1) = -dot_product(W(1, :), g)
           dx(2) = -dot_product(W(2, :), g)
           do k = 1, 2
              if (x(k) + dx(k) .lt. 0) then
                 dx(k) = -x(k)
              else if (x(k) + dx(k) .gt. 1) then
                 dx(k) = 1 - x(k)
              end if
           end do
           norm_dx = sqrt(dot_product(dx, dx))

           ! print *, counter, norm
           if ((norm_g .lt. 1e-13) .or. (norm_dx .lt. 1e-13)) then
              exit
           end if
           x = x + dx
        end do
        
        call basis0(ku, ku+mu, x(1), knot_u, bu0, offset_u)
        call basis0(kv, kv+mv, x(2), knot_v, bv0, offset_v)
        P(:) = 0.0
        do lu = 1, ku
           hu = lu + offset_u
           do lv = 1, kv
              hv = lv + offset_v
              P = P + bu0(lu) * bv0(lv) * cp(hu, hv, :)
           end do
        end do
        d = abs(dot_product(P - pt0, P - pt0))
        if (d .lt. min_d(ipt_proj)) then
           min_d(ipt_proj) = d
           min_s(ipt_proj) = isurf
           min_tu(ipt_proj) = x(1)
           min_tv(ipt_proj) = x(2)
        end if
     end do

     deallocate(knot_u)
     deallocate(knot_v)

     deallocate(param_u)
     deallocate(param_v)

     deallocate(bu0)
     deallocate(bu1)
     deallocate(bu2)
     deallocate(bv0)
     deallocate(bv1)
     deallocate(bv2)

     deallocate(cp)
  end do

end subroutine computeproj
