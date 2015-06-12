subroutine computesnnz(npt, nsurf, ngroup, &
     surf_group, order, surfs, nS)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) npt, nsurf, ngroup, surf_group, order, surfs
  !f2py intent(out) nS
  !f2py depend(nsurf) surf_group
  !f2py depend(ngroup) order
  !f2py depend(npt) surfs

  !Input
  integer, intent(in) ::  npt, nsurf, ngroup
  integer, intent(in) ::  surf_group(nsurf, 2)
  integer, intent(in) ::  order(ngroup)
  integer, intent(in) ::  surfs(npt)

  !Output
  integer, intent(out) ::  nS

  !Working
  integer ipt, isurf
  integer ku, kv
  logical found(nsurf)

  found(:) = .False.
  do ipt = 1, npt
     found(surfs(ipt)) = .True.
  end do

  nS = 0
  do isurf = 1, nsurf
     if (found(isurf)) then
        ku = order(surf_group(isurf, 1))
        kv = order(surf_group(isurf, 2))

        do ipt = 1, npt
           if (surfs(ipt) .eq. isurf) then
              nS = nS + ku * kv
           end if
        end do
     end if
  end do

end subroutine computesnnz



subroutine computesmtx( &
     uder, vder, nS, npt, nsurf, ngroup, &
     str_indices_cp, surf_group, &
     order, num_cp, surfs, ind_u, ind_v, &
     Sa, Si, Sj)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) uder, vder, nS, npt, nsurf, ngroup, str_indices_cp, surf_group, order, num_cp, surfs, ind_u, ind_v
  !f2py intent(out) Sa, Si, Sj
  !f2py depend(nsurf) str_indices_cp, surf_group
  !f2py depend(ngroup) order, num_cp
  !f2py depend(npt) surfs, ind_u, ind_v
  !f2py depend(nS) Sa, Si, Sj

  !Input
  integer, intent(in) ::  uder, vder, nS, npt, nsurf, ngroup
  integer, intent(in) ::  str_indices_cp(nsurf, 2)
  integer, intent(in) ::  surf_group(nsurf, 2)
  integer, intent(in) ::  order(ngroup), num_cp(ngroup)
  integer, intent(in) ::  surfs(npt)
  double precision, intent(in) ::  ind_u(npt), ind_v(npt)

  !Output
  double precision, intent(out) ::  Sa(nS)
  integer, intent(out) ::  Si(nS), Sj(nS)

  !Working
  integer iS, ipt, isurf
  integer ku, kv, mu, mv, u, v, lu, lv
  integer offset_cp, fu, fv
  logical found(nsurf)
  double precision tu, tv
  double precision, allocatable, dimension(:) ::  knot_u, knot_v
  double precision, allocatable, dimension(:) ::  bu, bv

  found(:) = .False.
  do ipt = 1, npt
     found(surfs(ipt)) = .True.
  end do

  iS = 0
  Sa(:) = 0.0
  Si(:) = 0
  Sj(:) = 0

  do isurf = 1, nsurf
     if (found(isurf)) then
        offset_cp = str_indices_cp(isurf, 1)

        ku = order(surf_group(isurf, 1))
        kv = order(surf_group(isurf, 2))
        mu = num_cp(surf_group(isurf, 1))
        mv = num_cp(surf_group(isurf, 2))

        allocate(knot_u(ku + mu))
        call knotopen(ku, ku+mu, knot_u)
        
        allocate(knot_v(kv + mv))
        call knotopen(kv, kv+mv, knot_v)

        allocate(bu(ku))
        allocate(bv(kv))

        do ipt = 1, npt
           if (surfs(ipt) .eq. isurf) then
              tu = ind_u(ipt)
              if (uder .eq. 0) then
                 call basis0(ku, ku+mu, tu, knot_u, bu, fu)
              else if (uder .eq. 1) then
                 call basis1(ku, ku+mu, tu, knot_u, bu, fu)
              else if (uder .eq. 2) then
                 call basis2(ku, ku+mu, tu, knot_u, bu, fu)
              end if
              
              tv = ind_v(ipt)
              if (vder .eq. 0) then
                 call basis0(kv, kv+mv, tv, knot_v, bv, fv)
              else if (vder .eq. 1) then
                 call basis1(kv, kv+mv, tv, knot_v, bv, fv)
              else if (vder .eq. 2) then
                 call basis2(kv, kv+mv, tv, knot_v, bv, fv)
              end if

              do lu = 1, ku
                 u = lu + fu
                 do lv = 1, kv
                    v = lv + fv

                    iS = iS + 1
                    Sa(iS) = bu(lu) * bv(lv)
                    Si(iS) = ipt
                    Sj(iS) = offset_cp + (v-1) * mu + u
                 end do
              end do
           end if
        end do

        deallocate(knot_u)
        deallocate(knot_v)

        deallocate(bu)
        deallocate(bv)
     end if
  end do

  if (iS .ne. nS) then
     print *, 'Error in computeSmtx', iS, nS
     call exit(1)
  end if

end subroutine computesmtx
