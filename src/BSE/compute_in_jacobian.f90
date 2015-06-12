subroutine computeEnnz( &
     nsurf, nedge, nvert, ngroup, &
     surf_ptrs, surf_group, &
     sm_vert_count, sm_edge_count, num_cp, nE)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nsurf, nedge, nvert, ngroup, surf_ptrs, surf_group, sm_vert_count, sm_edge_count, num_cp
  !f2py intent(out) nE
  !f2py depend(nsurf) surf_ptrs, surf_group
  !f2py depend(nvert) sm_vert_count
  !f2py depend(nedge) sm_edge_count
  !f2py depend(ngroup) num_cp

  !Input
  integer, intent(in) ::  nsurf, nedge, nvert, ngroup
  integer, intent(in) ::  surf_ptrs(nsurf, 3, 3)
  integer, intent(in) ::  surf_group(nsurf, 2)
  integer, intent(in) ::  sm_vert_count(nvert), sm_edge_count(nedge)
  integer, intent(in) ::  num_cp(ngroup)

  !Output
  integer, intent(out) ::  nE

  !Working
  integer isurf, iedge, ivert
  integer i, j
  integer mu, mv

  nE = 0

  ! Loop through surfs and add verts
  do isurf = 1, nsurf
     do i = 1, 3, 2
        do j = 1, 3, 2
           ivert = surf_ptrs(isurf, i, j)
           if (sm_vert_count(ivert) .eq. 0) then
              nE = nE + 1
           end if
        end do
     end do
  end do

  ! Loop through surfs and add edges
  do isurf = 1, nsurf
     mu = num_cp(surf_group(isurf, 1))
     mv = num_cp(surf_group(isurf, 2))

     iedge = surf_ptrs(isurf, 2, 1)
     if (sm_edge_count(abs(iedge)) .eq. 0) then
        nE = nE + mu-2
     end if

     iedge = surf_ptrs(isurf, 2, 3)
     if (sm_edge_count(abs(iedge)) .eq. 0) then
        nE = nE + mu-2
     end if

     iedge = surf_ptrs(isurf, 1, 2)
     if (sm_edge_count(abs(iedge)) .eq. 0) then
        nE = nE + mv-2
     end if

     iedge = surf_ptrs(isurf, 3, 2)
     if (sm_edge_count(abs(iedge)) .eq. 0) then
        nE = nE + mv-2
     end if
  end do

  ! Loop through surfaces and add interior pts
  do isurf = 1, nsurf
     mu = num_cp(surf_group(isurf, 1))
     mv = num_cp(surf_group(isurf, 2))
     nE = nE + (mu-2) * (mv-2)
  end do

end subroutine computeEnnz



subroutine computeEmtx( &
     nE, nsurf, nedge, nvert, ngroup, &
     surf_ptrs, surf_group, &
     surf_indices_df, edge_indices_df, &
     str_indices_df, vert_indices, &
     vert_count, edge_count, &
     sm_vert_count, sm_edge_count, num_cp, &
     Ea, Ei, Ej)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nE, nsurf, nedge, nvert, ngroup, surf_ptrs, surf_group, surf_indices_df, edge_indices_df, str_indices_df, vert_indices, vert_count, edge_count, sm_vert_count, sm_edge_count, num_cp
  !f2py intent(out) Ea, Ei, Ej
  !f2py depend(nsurf) surf_ptrs, surf_group, surf_indices_df
  !f2py depend(nedge) edge_indices_df
  !f2py depend(nsurf) str_indices_df
  !f2py depend(nvert) vert_count, sm_vert_count, vert_indices
  !f2py depend(nedge) edge_count, sm_edge_count
  !f2py depend(ngroup) num_cp

  !Input
  integer, intent(in) ::  nE, nsurf, nedge, nvert, ngroup
  integer, intent(in) ::  surf_ptrs(nsurf, 3, 3)
  integer, intent(in) ::  surf_group(nsurf, 2)
  integer, intent(in) ::  surf_indices_df(nsurf, 2)
  integer, intent(in) ::  edge_indices_df(nedge, 2)
  integer, intent(in) ::  str_indices_df(nsurf, 2)
  integer, intent(in) ::  vert_indices(nvert)
  integer, intent(in) ::  vert_count(nvert), edge_count(nedge)
  integer, intent(in) ::  sm_vert_count(nvert), sm_edge_count(nedge)
  integer, intent(in) ::  num_cp(ngroup)

  !Output
  double precision, intent(out) ::  Ea(nE)
  integer, intent(out) ::  Ei(nE), Ej(nE)

  !Working
  integer iE, isurf, iedge, ivert
  integer i_list(4), j_list(4), u_list(4), v_list(4), k
  integer i, j, u, v
  integer mu, mv

  iE = 0
  Ea(:) = 0.0
  Ei(:) = 0
  Ej(:) = 0

  ! Loop through surfs and add verts
  i_list = (/ 1, 3, 1, 3 /)
  j_list = (/ 1, 1, 3, 3 /)
  do isurf = 1, nsurf
     mu = num_cp(surf_group(isurf, 1))
     mv = num_cp(surf_group(isurf, 2))
     u_list = (/ 1, mu, 1, mu /)
     v_list = (/ 1, 1, mv, mv /)
     do k = 1, 4
        i = i_list(k)
        j = j_list(k)
        u = u_list(k)
        v = v_list(k)
        ivert = surf_ptrs(isurf, i, j)
        if (sm_vert_count(ivert) .eq. 0) then
           iE = iE + 1
           Ea(iE) = 1.0 / vert_count(ivert)
           Ei(iE) = vert_indices(ivert)
           Ej(iE) = str_indices_df(isurf, 1) &
                + (v-1) * mu + u
        end if
     end do
  end do

  ! Loop through surfs and add edges
  do isurf = 1, nsurf
     mu = num_cp(surf_group(isurf, 1))
     mv = num_cp(surf_group(isurf, 2))

     iedge = surf_ptrs(isurf, 2, 1)
     if (sm_edge_count(abs(iedge)) .eq. 0) then
        v = 1
        do k = 1, mu-2
           if (iedge .gt. 0) then
              u = k + 1
           else
              u = mu - k
           end if
           iE = iE + 1
           Ea(iE) = 1.0 / edge_count(abs(iedge))
           Ei(iE) = edge_indices_df(abs(iedge), 1) + k
           Ej(iE) = str_indices_df(isurf, 1) &
                + (v-1) * mu + u
        end do
     end if

     iedge = surf_ptrs(isurf, 2, 3)
     if (sm_edge_count(abs(iedge)) .eq. 0) then
        v = mv
        do k = 1, mu-2
           if (iedge .gt. 0) then
              u = k + 1
           else
              u = mu - k
           end if
           iE = iE + 1
           Ea(iE) = 1.0 / edge_count(abs(iedge))
           Ei(iE) = edge_indices_df(abs(iedge), 1) + k
           Ej(iE) = str_indices_df(isurf, 1) &
                + (v-1) * mu + u
        end do
     end if

     iedge = surf_ptrs(isurf, 1, 2)
     if (sm_edge_count(abs(iedge)) .eq. 0) then
        u = 1
        do k = 1, mv-2
           if (iedge .gt. 0) then
              v = k + 1
           else
              v = mv - k
           end if
           iE = iE + 1
           Ea(iE) = 1.0 / edge_count(abs(iedge))
           Ei(iE) = edge_indices_df(abs(iedge), 1) + k
           Ej(iE) = str_indices_df(isurf, 1) &
                + (v-1) * mu + u
        end do
     end if

     iedge = surf_ptrs(isurf, 3, 2)
     if (sm_edge_count(abs(iedge)) .eq. 0) then
        u = mu
        do k = 1, mv-2
           if (iedge .gt. 0) then
              v = k + 1
           else
              v = mv - k
           end if
           iE = iE + 1
           Ea(iE) = 1.0 / edge_count(abs(iedge))
           Ei(iE) = edge_indices_df(abs(iedge), 1) + k
           Ej(iE) = str_indices_df(isurf, 1) &
                + (v-1) * mu + u
        end do
     end if
  end do

  ! Loop through surfaces and add interior pts
  do isurf = 1, nsurf
     mu = num_cp(surf_group(isurf, 1))
     mv = num_cp(surf_group(isurf, 2))
     do u = 1, mu-2
        do v = 1, mv-2
           iE = iE + 1
           Ea(iE) = 1.0
           Ei(iE) = surf_indices_df(isurf, 1) &
                + (v-1) * (mu-2) + u
           Ej(iE) = str_indices_df(isurf, 1) &
                + (v) * mu + u+1
        end do
     end do
  end do

  if (iE .ne. nE) then
     print *, 'Error in computeEmtx', iE, nE
     call exit(1)
  end if

end subroutine computeEmtx
