subroutine computeDnnz(nsurf, nedge, nvert, ngroup, &
     surf_group, edge_group, sm_surf, sm_edge, &
     sm_vert_count, sm_edge_count, num_cp, &
     nD)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nsurf, nedge, nvert, ngroup, surf_group, edge_group, sm_surf, sm_edge, sm_vert_count, sm_edge_count, num_cp
  !f2py intent(out) nD
  !f2py depend(nsurf) surf_group, sm_surf
  !f2py depend(nedge) edge_group, sm_edge
  !f2py depend(nvert) sm_vert_count
  !f2py depend(nedge) sm_edge_count
  !f2py depend(ngroup) num_cp

  !Input
  integer, intent(in) ::  nsurf, nedge, nvert, ngroup
  integer, intent(in) ::  surf_group(nsurf, 2), edge_group(nedge)
  logical, intent(in) ::  sm_surf(nsurf, 3, 3), sm_edge(nedge, 2)
  integer, intent(in) ::  sm_vert_count(nvert), sm_edge_count(nedge)
  integer, intent(in) ::  num_cp(ngroup)

  !Output
  integer, intent(out) ::  nD

  !Working
  integer isurf, iedge, ivert
  integer m, mu, mv
  integer i, j

  nD = 0

  ! Verts: entries for free vertices
  do ivert = 1, nvert
     if (sm_vert_count(ivert) .eq. 0) then
        nD = nD + 1
     end if
  end do

  ! Verts: entries for vertices' dependence on edges
  do iedge = 1, nedge
     if (sm_edge(iedge, 1)) then
        nD = nD + 1
     end if
     if (sm_edge(iedge, 2)) then
        nD = nD + 1
     end if
  end do

  ! Verts: entries for vertices' dependence on surfaces
  do isurf = 1, nsurf
     do i = 1, 3, 2
        do j = 1, 3, 2
           if (sm_surf(isurf, i, j)) then
              nD = nD + 1
           end if
        end do
     end do
  end do

  ! Edges: entries for free edges
  do iedge = 1, nedge
     if (sm_edge_count(iedge) .eq. 0) then
        m = num_cp(edge_group(iedge))
        nD = nD + m-2
     end if
  end do

  ! Edges: entries for edges' dependence on surfaces
  do isurf = 1, nsurf
     mu = num_cp(surf_group(isurf, 1))
     mv = num_cp(surf_group(isurf, 2))
     if (sm_surf(isurf, 2, 1)) then
        nD = nD + mu-2
     end if
     if (sm_surf(isurf, 2, 3)) then
        nD = nD + mu-2
     end if
     if (sm_surf(isurf, 1, 2)) then
        nD = nD + mv-2
     end if
     if (sm_surf(isurf, 3, 2)) then
        nD = nD + mv-2
     end if
  end do

  ! Surfaces: entries for interior surface control points
  do isurf = 1, nsurf
     mu = num_cp(surf_group(isurf, 1))
     mv = num_cp(surf_group(isurf, 2))
     nD = nD + (mu-2) * (mv-2)
  end do

end subroutine computeDnnz



subroutine computeDmtx(nD, nsurf, nedge, nvert, ngroup, &
     surf_group, edge_group, surf_ptrs, edge_ptrs, &
     sm_surf, sm_edge, sm_vert_count, sm_edge_count, &
     surf_indices_df, surf_indices_cp, &
     edge_indices_df, edge_indices_cp, &
     vert_indices, num_cp, &
     Da, Di, Dj)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nD, nsurf, nedge, nvert, ngroup, surf_group, edge_group, surf_ptrs, edge_ptrs, sm_surf, sm_edge, sm_vert_count, sm_edge_count, surf_indices_df, surf_indices_cp, edge_indices_df, edge_indices_cp, vert_indices, num_cp
  !f2py intent(out) Da, Di, Dj
  !f2py depend(nsurf) surf_group, surf_ptrs, sm_surf
  !f2py depend(nedge) edge_group, edge_ptrs, sm_edge
  !f2py depend(nsurf) surf_indices_df, surf_indices_cp
  !f2py depend(nedge) edge_indices_df, edge_indices_cp
  !f2py depend(nedge) sm_edge_count
  !f2py depend(nvert) sm_vert_count, vert_indices
  !f2py depend(ngroup) num_cp
  !f2py depend(nD) Da, Di, Dj

  !Input
  integer, intent(in) ::  nD, nsurf, nedge, nvert, ngroup
  integer, intent(in) ::  surf_group(nsurf, 2), edge_group(nedge)
  integer, intent(in) ::  surf_ptrs(nsurf, 3, 3), edge_ptrs(nedge, 2)
  logical, intent(in) ::  sm_surf(nsurf, 3, 3), sm_edge(nedge, 2)
  integer, intent(in) ::  sm_vert_count(nvert), sm_edge_count(nedge)
  integer, intent(in) ::  surf_indices_df(nsurf, 2)
  integer, intent(in) ::  surf_indices_cp(nsurf, 2)
  integer, intent(in) ::  edge_indices_df(nedge, 2)
  integer, intent(in) ::  edge_indices_cp(nedge, 2)
  integer, intent(in) ::  vert_indices(nvert)
  integer, intent(in) ::  num_cp(ngroup)

  !Output
  double precision, intent(out) ::  Da(nD)
  integer, intent(out) ::  Di(nD), Dj(nD)

  !Working
  integer iD, isurf, iedge, ivert
  integer m, mu, mv
  integer i_list(4), j_list(4), u_list(4), v_list(4), k
  integer i, j, u, v

  iD = 0
  Da(:) = 0.0
  Di(:) = 0
  Dj(:) = 0

  ! Verts: entries for free vertices
  do ivert = 1, nvert
     if (sm_vert_count(ivert) .eq. 0) then
        iD = iD + 1
        Da(iD) = 1.0
        Di(iD) = ivert
        Dj(iD) = vert_indices(ivert)
     end if
  end do

  ! Verts: entries for vertices' dependence on edges
  do iedge = 1, nedge
     if (sm_edge(iedge, 1)) then
        ivert = edge_ptrs(iedge, 1)
        iD = iD + 1
        Da(iD) = 1.0 / sm_vert_count(ivert)
        Di(iD) = ivert
        Dj(iD) = edge_indices_df(iedge, 1) + 1
     end if
     if (sm_edge(iedge, 2)) then
        ivert = edge_ptrs(iedge, 2)
        iD = iD + 1
        Da(iD) = 1.0 / sm_vert_count(ivert)
        Di(iD) = ivert
        Dj(iD) = edge_indices_df(iedge, 2)
     end if
  end do

  ! Verts: entries for vertices' dependence on surfaces
  i_list = (/ 1, 3, 1, 3 /)
  j_list = (/ 1, 1, 3, 3 /)
  do isurf = 1, nsurf
     mu = num_cp(surf_group(isurf, 1))
     mv = num_cp(surf_group(isurf, 2))
     u_list = (/ 1, mu-2, 1, mu-2 /)
     v_list = (/ 1, 1, mv-2, mv-2 /)
     do k = 1, 4
        i = i_list(k)
        j = j_list(k)
        u = u_list(k)
        v = v_list(k)
        if (sm_surf(isurf, i, j)) then
           ivert = surf_ptrs(isurf, i, j)
           iD = iD + 1
           Da(iD) = 1.0 / sm_vert_count(ivert)
           Di(iD) = ivert
           Dj(iD) = surf_indices_df(isurf, 1) &
                + (v-1)*(mu-2) + u
        end if
     end do
  end do

  ! Edges: entries for free edges
  do iedge = 1, nedge
     if (sm_edge_count(iedge) .eq. 0) then
        m = num_cp(edge_group(iedge))
        do k = 1, m-2
           iD = iD + 1
           Da(iD) = 1.0
           Di(iD) = edge_indices_cp(iedge, 1) + k
           Dj(iD) = edge_indices_df(iedge, 1) + k
        end do
     end if
  end do

  ! Edges: entries for edges' dependence on surfaces
  do isurf = 1, nsurf
     mu = num_cp(surf_group(isurf, 1))
     mv = num_cp(surf_group(isurf, 2))
     if (sm_surf(isurf, 2, 1)) then
        iedge = surf_ptrs(isurf, 2, 1)
        v = 1
        do k = 1, mu-2
           if (iedge .gt. 0) then
              u = k
           else
              u = mu-1 - k
           end if
           iD = iD + 1
           Da(iD) = 1.0 / sm_edge_count(abs(iedge))
           Di(iD) = edge_indices_cp(abs(iedge), 1) + k
           Dj(iD) = surf_indices_df(isurf, 1) &
                + (v-1) * (mu-2) + u
        end do
     end if
     if (sm_surf(isurf, 2, 3)) then
        iedge = surf_ptrs(isurf, 2, 3)
        v = mv - 2
        do k = 1, mu-2
           if (iedge .gt. 0) then
              u = k
           else
              u = mu-1 - k
           end if
           iD = iD + 1
           Da(iD) = 1.0 / sm_edge_count(abs(iedge))
           Di(iD) = edge_indices_cp(abs(iedge), 1) + k
           Dj(iD) = surf_indices_df(isurf, 1) &
                + (v-1) * (mu-2) + u
        end do
     end if
     if (sm_surf(isurf, 1, 2)) then
        iedge = surf_ptrs(isurf, 1, 2)
        u = 1
        do k = 1, mv-2
           if (iedge .gt. 0) then
              v = k
           else
              v = mv-1 - k
           end if
           iD = iD + 1
           Da(iD) = 1.0 / sm_edge_count(abs(iedge))
           Di(iD) = edge_indices_cp(abs(iedge), 1) + k
           Dj(iD) = surf_indices_df(isurf, 1) &
                + (v-1) * (mu-2) + u
        end do
     end if
     if (sm_surf(isurf, 3, 2)) then
        iedge = surf_ptrs(isurf, 3, 2)
        u = mu - 2
        do k = 1, mv-2
           if (iedge .gt. 0) then
              v = k
           else
              v = mv-1 - k
           end if
           iD = iD + 1
           Da(iD) = 1.0 / sm_edge_count(abs(iedge))
           Di(iD) = edge_indices_cp(abs(iedge), 1) + k
           Dj(iD) = surf_indices_df(isurf, 1) &
                + (v-1) * (mu-2) + u
        end do
     end if
  end do

  ! Surfaces: entries for interior surface control points
  do isurf = 1, nsurf
     mu = num_cp(surf_group(isurf, 1))
     mv = num_cp(surf_group(isurf, 2))
     do u = 1, mu-2
        do v = 1, mv-2
           iD = iD + 1
           Da(iD) = 1.0
           Di(iD) = surf_indices_cp(isurf, 1) &
                + (v-1)*(mu-2) + u
           Dj(iD) = surf_indices_df(isurf, 1) &
                + (v-1)*(mu-2) + u
        end do
     end do
  end do

  if (iD .ne. nD) then
     print *, 'Error in computeFmtx', iD, nD
     call exit(1)
  end if

end subroutine computeDmtx
