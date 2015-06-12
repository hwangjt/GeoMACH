subroutine computebspline(nB, &
     ku, kv, mu, mv, nu, nv, &
     cp_indices, pt_indices, &
     pos_u, pos_v, Ba, Bi, Bj)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nB, ku, kv, mu, mv, nu, nv, cp_indices, pt_indices, pos_u, pos_v
  !f2py intent(out) Ba, Bi, Bj
  !f2py depend(mu, mv) cp_indices
  !f2py depend(nu, nv) pt_indices
  !f2py depend(mu) pos_u
  !f2py depend(mv) pos_v
  !f2py depend(nB) Ba, Bi, Bj

  !Input
  integer, intent(in) ::  nB, ku, kv, mu, mv, nu, nv
  integer, intent(in) ::  cp_indices(mu, mv)
  integer, intent(in) ::  pt_indices(nu, nv)
  double precision, intent(in) ::  pos_u(mu), pos_v(mv)

  !Output
  double precision, intent(out) ::  Ba(nB)
  integer, intent(out) ::  Bi(nB), Bj(nB)

  !Working
  integer iB, cp_u, cp_v, pt_u, pt_v, ind_u, ind_v
  double precision knot_u(ku + mu), knot_v(kv + mv)
  double precision param_u(nu), param_v(nv)
  double precision basis_u(ku, nu), basis_v(kv, nv)
  integer offset_u(nu), offset_v(nv)

  iB = 0
  Ba(:) = 0.0
  Bi(:) = 0
  Bj(:) = 0

  if (mu .eq. 1) then
     basis_u(:, :) = 1.0
     offset_u(:) = 0
  else
     call knotopen(ku, ku + mu, knot_u)
     call param_custom(ku + mu, mu, nu, knot_u, pos_u, param_u)
     do pt_u = 1, nu
        if (param_u(pt_u) .ge. 0) then
           call basis(ku, ku + mu, param_u(pt_u), knot_u, &
                basis_u(:, pt_u), offset_u(pt_u))
        else
           basis_u(:, pt_u) = 0.0
           offset_u(pt_u) = 0
        end if
     end do
  end if

  if (mv .eq. 1) then
     basis_v(:, :) = 1.0
     offset_v(:) = 0
  else
     call knotopen(kv, kv + mv, knot_v)
     call param_custom(kv + mv, mv, nv, knot_v, pos_v, param_v)
     do pt_v = 1, nv
        if (param_v(pt_v) .ge. 0) then
           call basis(kv, kv + mv, param_v(pt_v), knot_v, &
                basis_v(:, pt_v), offset_v(pt_v))
        else
           basis_v(:, pt_v) = 0.0
           offset_v(pt_v) = 0
        end if
     end do
  end if

  do pt_v = 1, nv
     do pt_u = 1, nu
        do ind_v = 1, kv
           cp_v = offset_v(pt_v) + ind_v
           do ind_u = 1, ku
              cp_u = offset_u(pt_u) + ind_u

              iB = iB + 1
              Ba(iB) = &
                   basis_u(ind_u, pt_u) * &
                   basis_v(ind_v, pt_v)
              Bi(iB) = pt_indices(pt_u, pt_v)
              Bj(iB) = cp_indices(cp_u, cp_v)
           end do
        end do
     end do
  end do

  if (iB .ne. nB) then
     print *, 'Error in computeBspline', iB, nB
     call exit(1)
  end if

end subroutine computebspline
