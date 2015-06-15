subroutine computePmtx( &
     nP, nsurf, nedge, nvert, ngroup, &
     surf_ptrs, surf_group, surf_indices_pt, &
     edge_indices_pt, str_indices_pt, &
     vert_count, edge_count, num_pt, &
     Pa, Pi, Pj)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nP, nsurf, nedge, nvert, ngroup, surf_ptrs, surf_group, surf_indices_pt, edge_indices_pt, str_indices_pt, vert_count, edge_count, num_pt
  !f2py intent(out) Pa, Pi, Pj
  !f2py depend(nsurf) surf_ptrs, surf_group
  !f2py depend(nsurf) surf_indices_pt, str_indices_pt
  !f2py depend(nedge) edge_indices_pt, edge_count
  !f2py depend(nvert) vert_count
  !f2py depend(ngroup) num_pt

  !Input
  integer, intent(in) ::  nP, nsurf, nedge, nvert, ngroup
  integer, intent(in) ::  surf_ptrs(nsurf, 3, 3)
  integer, intent(in) ::  surf_group(nsurf, 2)
  integer, intent(in) ::  surf_indices_pt(nsurf, 2)
  integer, intent(in) ::  edge_indices_pt(nedge, 2)
  integer, intent(in) ::  str_indices_pt(nsurf, 2)
  integer, intent(in) ::  vert_count(nvert), edge_count(nedge)
  integer, intent(in) ::  num_pt(ngroup)

  !Output
  double precision, intent(out) ::  Pa(nP)
  integer, intent(out) ::  Pi(nP), Pj(nP)

  !Working
  integer iP, isurf, iedge, ivert
  integer i_list(4), j_list(4), u_list(4), v_list(4), k
  integer i, j, u, v
  integer nu, nv

  iP = 0
  Pa(:) = 0.0
  Pi(:) = 0
  Pj(:) = 0

  ! Loop through surfs and add verts
  i_list = (/ 1, 3, 1, 3 /)
  j_list = (/ 1, 1, 3, 3 /)
  do isurf = 1, nsurf
     nu = num_pt(surf_group(isurf, 1))
     nv = num_pt(surf_group(isurf, 2))
     u_list = (/ 1, nu, 1, nu /)
     v_list = (/ 1, 1, nv, nv /)
     do k = 1, 4
        i = i_list(k)
        j = j_list(k)
        u = u_list(k)
        v = v_list(k)
        ivert = surf_ptrs(isurf, i, j)
        iP = iP + 1
        Pa(iP) = 1.0 / vert_count(ivert)
        Pi(iP) = ivert
        Pj(iP) = str_indices_pt(isurf, 1) &
             + (v-1) * nu + u
     end do
  end do

  ! Loop through surfs and add edges
  do isurf = 1, nsurf
     nu = num_pt(surf_group(isurf, 1))
     nv = num_pt(surf_group(isurf, 2))

     iedge = surf_ptrs(isurf, 2, 1)
     v = 1
     do k = 1, nu-2
        if (iedge .gt. 0) then
           u = k + 1
        else
           u = nu - k
        end if
        iP = iP + 1
        Pa(iP) = 1.0 / edge_count(abs(iedge))
        Pi(iP) = edge_indices_pt(abs(iedge), 1) + k
        Pj(iP) = str_indices_pt(isurf, 1) &
             + (v-1) * nu + u
     end do

     iedge = surf_ptrs(isurf, 2, 3)
     v = nv
     do k = 1, nu-2
        if (iedge .gt. 0) then
           u = k + 1
        else
           u = nu - k
        end if
        iP = iP + 1
        Pa(iP) = 1.0 / edge_count(abs(iedge))
        Pi(iP) = edge_indices_pt(abs(iedge), 1) + k
        Pj(iP) = str_indices_pt(isurf, 1) &
             + (v-1) * nu + u
     end do

     iedge = surf_ptrs(isurf, 1, 2)
     u = 1
     do k = 1, nv-2
        if (iedge .gt. 0) then
           v = k + 1
        else
           v = nv - k
        end if
        iP = iP + 1
        Pa(iP) = 1.0 / edge_count(abs(iedge))
        Pi(iP) = edge_indices_pt(abs(iedge), 1) + k
        Pj(iP) = str_indices_pt(isurf, 1) &
             + (v-1) * nu + u
     end do

     iedge = surf_ptrs(isurf, 3, 2)
     u = nu
     do k = 1, nv-2
        if (iedge .gt. 0) then
           v = k + 1
        else
           v = nv - k
        end if
        iP = iP + 1
        Pa(iP) = 1.0 / edge_count(abs(iedge))
        Pi(iP) = edge_indices_pt(abs(iedge), 1) + k
        Pj(iP) = str_indices_pt(isurf, 1) &
             + (v-1) * nu + u
     end do
  end do

  ! Loop through surfaces and add interior pts
  do isurf = 1, nsurf
     nu = num_pt(surf_group(isurf, 1))
     nv = num_pt(surf_group(isurf, 2))
     do u = 1, nu-2
        do v = 1, nv-2
           iP = iP + 1
           Pa(iP) = 1.0
           Pi(iP) = surf_indices_pt(isurf, 1) &
                + (v-1) * (nu-2) + u
           Pj(iP) = str_indices_pt(isurf, 1) &
                + (v) * nu + u+1
        end do
     end do
  end do

  if (iP .ne. nP) then
     print *, 'Error in computePmtx', iP, nP
     call exit(1)
  end if

end subroutine computePmtx
