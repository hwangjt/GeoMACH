subroutine computevertindices(nsurf, nedge, nvert, &
     surf_ptrs, edge_ptrs, sm_surf, sm_edge, vert_indices)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nsurf, nedge, nvert, surf_ptrs, edge_ptrs, sm_surf, sm_edge
  !f2py intent(out) vert_indices
  !f2py depend(nsurf) surf_ptrs, sm_surf
  !f2py depend(nedge) edge_ptrs, sm_edge
  !f2py depend(nvert) vert_indices

  !Input
  integer, intent(in) ::  nsurf, nedge, nvert
  integer, intent(in) ::  surf_ptrs(nsurf, 3, 3)
  integer, intent(in) ::  edge_ptrs(nedge, 2)
  logical, intent(in) ::  sm_surf(nsurf, 3, 3)
  logical, intent(in) ::  sm_edge(nedge, 2)

  !Output
  integer, intent(out) ::  vert_indices(nvert)

  !Working
  integer isurf, iedge, ivert, ivert_df, i, j, k
  logical dof(nvert)

  ! Assume it's a dof until at least one non-dof found
  dof(:) = .true.

  do isurf = 1, nsurf
     do i = 1, 3, 2
        do j = 1, 3, 2
           if (sm_surf(isurf, i, j)) then
              dof(surf_ptrs(isurf, i, j)) = .false.
           end if
        end do
     end do
  end do

  do iedge = 1, nedge
     do k = 1, 2
        if (sm_edge(iedge, k)) then
           dof(edge_ptrs(iedge, k)) = .false.
        end if
     end do
  end do

  ! For each dof, assigned unique vert ID
  vert_indices(:) = 0
  ivert_df = 0
  do ivert = 1, nvert
     if (dof(ivert)) then
        ivert_df = ivert_df + 1
        vert_indices(ivert) = ivert_df
     end if
  end do

end subroutine computevertindices



subroutine computeedgeindices(nsurf, nedge, ngroup, &
     surf_ptrs, sm_surf, edge_group, num_cp, num_pt, &
     edge_indices_df, edge_indices_cp, edge_indices_pt)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nsurf, nedge, ngroup, surf_ptrs, sm_surf, edge_group, num_cp, num_pt
  !f2py intent(out) edge_indices_df, edge_indices_cp, edge_indices_pt
  !f2py depend(nsurf) surf_ptrs, sm_surf
  !f2py depend(ngroup) num_cp, num_pt
  !f2py depend(nedge) edge_group, edge_indices_df, edge_indices_cp, edge_indices_pt

  !Input
  integer, intent(in) ::  nsurf, nedge, ngroup
  integer, intent(in) ::  surf_ptrs(nsurf, 3, 3)
  logical, intent(in) ::  sm_surf(nsurf, 3, 3)
  integer, intent(in) ::  edge_group(nedge)
  integer, intent(in) ::  num_cp(ngroup), num_pt(ngroup)

  !Output
  integer, intent(out) ::  edge_indices_df(nedge, 2)
  integer, intent(out) ::  edge_indices_cp(nedge, 2)
  integer, intent(out) ::  edge_indices_pt(nedge, 2)

  !Working
  integer isurf, iedge, ind_df, ind_cp, ind_pt
  integer i_list(4), j_list(4), i, j, k
  logical dof(nedge)

  ! Assume it's a dof until at least one non-dof found
  dof(:) = .true.

  i_list = (/ 2, 2, 1, 3 /)
  j_list = (/ 1, 3, 2, 2 /)
  do isurf = 1, nsurf
     do k = 1, 4
        i = i_list(k)
        j = j_list(k)
        if (sm_surf(isurf, i, j)) then
           dof(abs(surf_ptrs(isurf, i, j))) = .false.
        end if
     end do
  end do

  ! Loop through edges and assemble edge indices maps
  edge_indices_df(:, :) = 0
  edge_indices_cp(:, :) = 0
  edge_indices_pt(:, :) = 0

  ind_df = 0
  ind_cp = 0
  ind_pt = 0
  do iedge = 1, nedge
     if (dof(iedge)) then
        edge_indices_df(iedge, 1) = ind_df
        ind_df = ind_df + num_cp(edge_group(iedge)) - 2
        edge_indices_df(iedge, 2) = ind_df
     end if

     edge_indices_cp(iedge, 1) = ind_cp
     ind_cp = ind_cp + num_cp(edge_group(iedge)) - 2
     edge_indices_cp(iedge, 2) = ind_cp

     edge_indices_pt(iedge, 1) = ind_pt
     ind_pt = ind_pt + num_pt(edge_group(iedge)) - 2
     edge_indices_pt(iedge, 2) = ind_pt
  end do

end subroutine computeedgeindices



subroutine computesurfindices(nsurf, ngroup, &
     surf_group, num_cp, num_pt, &
     surf_indices_df, surf_indices_cp, surf_indices_pt)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nsurf, ngroup, surf_group, num_cp, num_pt
  !f2py intent(out) surf_indices_df, surf_indices_cp, surf_indices_pt
  !f2py depend(nsurf) surf_group
  !f2py depend(ngroup) num_cp, num_pt
  !f2py depend(nsurf) surf_indices_df, surf_indices_cp, surf_indices_pt

  !Input
  integer, intent(in) ::  nsurf, ngroup
  integer, intent(in) ::  surf_group(nsurf, 2)
  integer, intent(in) ::  num_cp(ngroup), num_pt(ngroup)

  !Output
  integer, intent(out) ::  surf_indices_df(nsurf, 2)
  integer, intent(out) ::  surf_indices_cp(nsurf, 2)
  integer, intent(out) ::  surf_indices_pt(nsurf, 2)

  !Working
  integer isurf, ind_df, ind_cp, ind_pt
  integer mu, mv, nu, nv

  ind_df = 0
  ind_cp = 0
  ind_pt = 0
  do isurf = 1, nsurf
     mu = num_cp(surf_group(isurf, 1))
     mv = num_cp(surf_group(isurf, 2))
     nu = num_pt(surf_group(isurf, 1))
     nv = num_pt(surf_group(isurf, 2))

     surf_indices_df(isurf, 1) = ind_df
     surf_indices_cp(isurf, 1) = ind_cp
     surf_indices_pt(isurf, 1) = ind_pt
     ind_df = ind_df + (mu - 2) * (mv - 2)
     ind_cp = ind_cp + (mu - 2) * (mv - 2)
     ind_pt = ind_pt + (nu - 2) * (nv - 2)
     surf_indices_df(isurf, 2) = ind_df
     surf_indices_cp(isurf, 2) = ind_cp
     surf_indices_pt(isurf, 2) = ind_pt
  end do

end subroutine computesurfindices



subroutine computestrindices(nsurf, ngroup, &
     surf_group, num_cp, num_pt, &
     str_indices_df, str_indices_cp, str_indices_pt)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nsurf, ngroup, surf_group, num_cp, num_pt
  !f2py intent(out) str_indices_df, str_indices_cp, str_indices_pt
  !f2py depend(nsurf) surf_group
  !f2py depend(ngroup) num_cp, num_pt
  !f2py depend(nsurf) str_indices_df, str_indices_cp, str_indices_pt

  !Input
  integer, intent(in) ::  nsurf, ngroup
  integer, intent(in) ::  surf_group(nsurf, 2)
  integer, intent(in) ::  num_cp(ngroup), num_pt(ngroup)

  !Output
  integer, intent(out) ::  str_indices_df(nsurf, 2)
  integer, intent(out) ::  str_indices_cp(nsurf, 2)
  integer, intent(out) ::  str_indices_pt(nsurf, 2)

  !Working
  integer isurf, ind_df, ind_cp, ind_pt
  integer mu, mv, nu, nv

  ind_df = 0
  ind_cp = 0
  ind_pt = 0
  do isurf = 1, nsurf
     mu = num_cp(surf_group(isurf, 1))
     mv = num_cp(surf_group(isurf, 2))
     nu = num_pt(surf_group(isurf, 1))
     nv = num_pt(surf_group(isurf, 2))

     str_indices_df(isurf, 1) = ind_df
     str_indices_cp(isurf, 1) = ind_cp
     str_indices_pt(isurf, 1) = ind_pt
     ind_df = ind_df + mu * mv
     ind_cp = ind_cp + mu * mv
     ind_pt = ind_pt + nu * nv
     str_indices_df(isurf, 2) = ind_df
     str_indices_cp(isurf, 2) = ind_cp
     str_indices_pt(isurf, 2) = ind_pt
  end do

end subroutine computestrindices




subroutine computemults( &
     nsurf, nedge, nvert, &
     surf_ptrs, edge_ptrs, &
     sm_surf, sm_edge, &
     vert_count, edge_count, &
     sm_vert_count, sm_edge_count)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nsurf, nedge, nvert, surf_ptrs, edge_ptrs, sm_surf, sm_edge
  !f2py intent(out) vert_count, edge_count, sm_vert_count, sm_edge_count
  !f2py depend(nsurf) surf_ptrs, sm_surf
  !f2py depend(nedge) edge_ptrs, sm_edge
  !f2py depend(nvert) vert_count, sm_vert_count
  !f2py depend(nedge) edge_count, sm_edge_count

  !Input
  integer, intent(in) ::  nsurf, nedge, nvert
  integer, intent(in) ::  surf_ptrs(nsurf, 3, 3)
  integer, intent(in) ::  edge_ptrs(nedge, 2)
  logical, intent(in) ::  sm_surf(nsurf, 3, 3)
  logical, intent(in) ::  sm_edge(nedge, 2)

  !Output
  integer, intent(out) ::  vert_count(nvert), edge_count(nedge)
  integer, intent(out) ::  sm_vert_count(nvert), sm_edge_count(nedge)

  !Working
  integer isurf, ivert, iedge, i_list(4), j_list(4), i, j, k

  vert_count(:) = 0
  edge_count(:) = 0
  sm_vert_count(:) = 0
  sm_edge_count(:) = 0

  ! Loop over surfaces and verts, sum appearances of each vertex
  i_list = (/ 1, 3, 1, 3 /)
  j_list = (/ 1, 1, 3, 3 /)
  do isurf = 1, nsurf
     do k = 1, 4
        i = i_list(k)
        j = j_list(k)
        ivert = surf_ptrs(isurf, i, j)
        vert_count(ivert) = vert_count(ivert) + 1
        if (sm_surf(isurf, i, j)) then
           sm_vert_count(ivert) = sm_vert_count(ivert) + 1
        end if
     end do
  end do

  ! Loop over surfaces and edges, sum appearances of each edge
  i_list = (/ 2, 2, 1, 3 /)
  j_list = (/ 1, 3, 2, 2 /)
  do isurf = 1, nsurf
     do k = 1, 4
        i = i_list(k)
        j = j_list(k)
        iedge = abs(surf_ptrs(isurf, i, j))
        edge_count(iedge) = edge_count(iedge) + 1
        if (sm_surf(isurf, i, j)) then
           sm_edge_count(iedge) = sm_edge_count(iedge) + 1
        end if
     end do
  end do

  ! Loop over edges, sum appearance of each vertex
  do iedge = 1, nedge
     do k = 1, 2
        if (sm_edge(iedge, k)) then
           ivert = edge_ptrs(iedge, k)
           sm_vert_count(ivert) = sm_vert_count(ivert) + 1
        end if
     end do
  end do

end subroutine computemults
