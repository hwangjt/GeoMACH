subroutine computesurfconnectivities( &
     nsurf, vtol, etol, surfaces, &
     nvert, nedge, surf_ptrs)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nsurf, vtol, etol, surfaces
  !f2py intent(out) nvert, nedge, surf_ptrs
  !f2py depend(nsurf) surfaces, surf_ptrs

  !Input
  integer, intent(in) ::  nsurf
  double precision, intent(in) ::  vtol, etol
  double precision, intent(in) ::  surfaces(nsurf, 3, 3, 3)

  !Output
  integer, intent(out) ::  nvert, nedge
  integer, intent(out) ::  surf_ptrs(nsurf, 3, 3)

  !Working
  integer isurf1, isurf2
  integer ivert, iedge
  integer i_list(4), j_list(4), di_list(4), dj_list(4)
  integer i1, i2, j1, j2, di1, di2, dj1, dj2, k1, k2
  integer i_m1, i_m2, i_p1, i_p2
  double precision norm_sq, compute_dist_sq

  surf_ptrs(:, :, :) = 0

  i_list = (/ 1, 3, 1, 3 /)
  j_list = (/ 1, 1, 3, 3 /)
  ivert = 0
  ! Loop over all surfaces and their vertices
  do isurf1 = 1, nsurf
     do k1 = 1, 4
        i1 = i_list(k1)
        j1 = j_list(k1)
        ! If vertex 1 is new, assign ID ++ivert to vertex 1
        if (surf_ptrs(isurf1, i1, j1) .eq. 0) then
           ivert = ivert + 1
           surf_ptrs(isurf1, i1, j1) = ivert
           ! Loop over all remaining surfaces and their vertices
           do isurf2 = isurf1 + 1, nsurf
              do k2 = 1, 4
                 i2 = i_list(k2)
                 j2 = j_list(k2)
                 ! If vertex 2 is new, compare vertices 1 and 2
                 if (surf_ptrs(isurf2, i2, j2) .eq. 0) then
                    norm_sq = compute_dist_sq( &
                         surfaces(isurf1, i1, j1, :), &
                         surfaces(isurf2, i2, j2, :))
                    ! If they match, assign ivert to vertex 2
                    if (norm_sq .lt. vtol) then
                       surf_ptrs(isurf2, i2, j2) = ivert
                    end if
                 end if
              end do
           end do
        end if
     end do
  end do
  nvert = ivert

  i_list  = (/ 1, 3, 2, 2 /)
  j_list  = (/ 2, 2, 1, 3 /)
  di_list = (/ 0, 0, 1, 1 /)
  dj_list = (/ 1, 1, 0, 0 /)
  iedge = 0
  ! Loop over all surfaces and their edges
  do isurf1 = 1, nsurf
     do k1 = 1, 4
        i1 = i_list(k1)
        j1 = j_list(k1)
        di1 = di_list(k1)
        dj1 = dj_list(k1)
        ! If edge 1 is new, assign ID ++iedge to edge 1
        if (surf_ptrs(isurf1, i1, j1) .eq. 0) then
           iedge = iedge + 1
           surf_ptrs(isurf1, i1, j1) = iedge
           ! Loop over all remaining surfaces and their edges
           do isurf2 = isurf1 + 1, nsurf
              do k2 = 1, 4
                 i2 = i_list(k2)
                 j2 = j_list(k2)
                 di2 = di_list(k2)
                 dj2 = dj_list(k2)
                 ! If edge 2 is new, compare edges 1 and 2
                 if (surf_ptrs(isurf2, i2, j2) .eq. 0) then
                    norm_sq = compute_dist_sq(&
                         surfaces(isurf1, i1, j1, :), &
                         surfaces(isurf2, i2, j2, :))
                    if (norm_sq .lt. etol) then
                       i_m1 = surf_ptrs(isurf1, i1 - di1, j1 - dj1)
                       i_m2 = surf_ptrs(isurf2, i2 - di2, j2 - dj2)
                       i_p1 = surf_ptrs(isurf1, i1 + di1, j1 + dj1)
                       i_p2 = surf_ptrs(isurf2, i2 + di2, j2 + dj2)

                       if ((i_m1 .eq. i_m2) .and. (i_p1 .eq. i_p2)) then
                          surf_ptrs(isurf2, i2, j2) = iedge
                       end if
                       if ((i_m1 .eq. i_p2) .and. (i_p1 .eq. i_m2)) then
                          surf_ptrs(isurf2, i2, j2) = -iedge
                       end if
                    end if
                 end if
              end do
           end do
        end if
     end do
  end do
  nedge = iedge

end subroutine computesurfconnectivities



subroutine computeedgeconnectivities( &
     nsurf, nedge, surf_ptrs, edge_ptrs)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nsurf, nedge, surf_ptrs
  !f2py intent(out) edge_ptrs
  !f2py depend(nsurf) surf_ptrs
  !f2py depend(nedge) edge_ptrs

  !Input
  integer, intent(in) ::  nsurf, nedge
  integer, intent(in) ::  surf_ptrs(nsurf, 3, 3)

  !Output
  integer, intent(out) ::  edge_ptrs(nedge, 2)

  !Working
  integer i_list(4), j_list(4), di_list(4), dj_list(4), k
  integer i, j, i_m, i_p, j_m, j_p
  integer isurf, iedge

  edge_ptrs(:, :) = 0

  i_list  = (/ 1, 3, 2, 2 /)
  j_list  = (/ 2, 2, 1, 3 /)
  di_list = (/ 0, 0, 1, 1 /)
  dj_list = (/ 1, 1, 0, 0 /)
  do isurf = 1, nsurf
     do k = 1, 4
        i = i_list(k)
        j = j_list(k)
        iedge = surf_ptrs(isurf, i, j)
        i_m = i - di_list(k)
        i_p = i + di_list(k)
        j_m = j - dj_list(k)
        j_p = j + dj_list(k)
        if (iedge .gt. 0) then
           edge_ptrs(abs(iedge), 1) = surf_ptrs(isurf, i_m, j_m)
           edge_ptrs(abs(iedge), 2) = surf_ptrs(isurf, i_p, j_p)
        else
           edge_ptrs(abs(iedge), 1) = surf_ptrs(isurf, i_p, j_p)
           edge_ptrs(abs(iedge), 2) = surf_ptrs(isurf, i_m, j_m)
        end if
     end do
  end do

end subroutine computeedgeconnectivities



function compute_dist_sq(vec1, vec2)

  implicit none

  !Input
  double precision, intent(in) ::  vec1(3), vec2(3)

  !Output
  double precision ::  compute_dist_sq

  ! Compute the square of the norm of vec2 - vec1
  compute_dist_sq = &
       (vec2(1) - vec1(1)) * (vec2(1) - vec1(1)) + &
       (vec2(2) - vec1(2)) * (vec2(2) - vec1(2)) + &
       (vec2(3) - vec1(3)) * (vec2(3) - vec1(3))

end function compute_dist_sq



subroutine computegroups( &
     nsurf, nedge, surf_ptrs, &
     ngroup, surf_group, edge_group)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nsurf, nedge, surf_ptrs
  !f2py intent(out) ngroup, surf_group, edge_group
  !f2py depend(nsurf) surf_ptrs, surf_group
  !f2py depend(nedge) edge_group

  !Input
  integer, intent(in) ::  nsurf, nedge
  integer, intent(in) ::  surf_ptrs(nsurf, 3, 3)

  !Output
  integer, intent(out) ::  ngroup
  integer, intent(out) ::  surf_group(nsurf, 2)
  integer, intent(out) ::  edge_group(nedge)

  !Working
  integer group_u0, group_u1, group_v0, group_v1
  integer iedge, isurf, igroup
  integer group_id(nedge)

  ! First, assume each edge has its own unique group
  do iedge = 1, nedge
     edge_group(iedge) = iedge
  end do

  ! Loop over all surfaces, making opposite edges have the same group
  do isurf = 1, nsurf
     group_v0 = edge_group(abs(surf_ptrs(isurf, 2, 1)))
     group_v1 = edge_group(abs(surf_ptrs(isurf, 2, 3)))
     group_u0 = edge_group(abs(surf_ptrs(isurf, 1, 2)))
     group_u1 = edge_group(abs(surf_ptrs(isurf, 3, 2)))

     if (group_v0 .ne. group_v1) then
        call set_all(nedge, group_v0, group_v1, edge_group)
     end if

     if (group_u0 .ne. group_u1) then
        call set_all(nedge, group_u0, group_u1, edge_group)
     end if
  end do

  ! Compute mapping (group_id) to remove gaps in group numbering
  group_id(:) = 0
  igroup = 0
  do iedge = 1, nedge
     if (group_id(edge_group(iedge)) .eq. 0) then
        igroup = igroup + 1
        group_id(edge_group(iedge)) = igroup
     end if
  end do
  ngroup = igroup

  ! Apply this mapping
  do iedge = 1, nedge
     edge_group(iedge) = group_id(edge_group(iedge))
  end do

  ! Map edge_group to surf_group
  surf_group(:, :) = 0
  do isurf = 1, nsurf
     surf_group(isurf, 1) = edge_group(abs(surf_ptrs(isurf, 2, 1)))
     surf_group(isurf, 2) = edge_group(abs(surf_ptrs(isurf, 1, 2)))
  end do

end subroutine computegroups



subroutine set_all(nedge, old_group, new_group, edge_group)

  implicit none

  !Input
  integer, intent(in) ::  nedge, old_group, new_group

  !Output
  integer, intent(inout) ::  edge_group(nedge)

  !Working
  integer iedge

  do iedge = 1, nedge
     if (edge_group(iedge) .eq. old_group) then
        edge_group(iedge) = new_group
     end if
  end do

end subroutine set_all
