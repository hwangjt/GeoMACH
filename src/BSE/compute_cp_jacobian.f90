subroutine computeCmtx( &
     nC, nsurf, nedge, ngroup, &
     surf_ptrs, surf_group, surf_indices_cp, &
     edge_indices_cp, str_indices_cp, num_cp, &
     Ca, Ci, Cj)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nC, nsurf, nedge, ngroup, surf_ptrs, surf_group, surf_indices_cp, edge_indices_cp, str_indices_cp, num_cp
  !f2py intent(out) Ca, Ci, Cj
  !f2py depend(nsurf) surf_ptrs, surf_group
  !f2py depend(nsurf) surf_indices_cp, str_indices_cp
  !f2py depend(nedge) edge_indices_cp
  !f2py depend(ngroup) num_cp

  !Input
  integer, intent(in) ::  nC, nsurf, nedge, ngroup
  integer, intent(in) ::  surf_ptrs(nsurf, 3, 3)
  integer, intent(in) ::  surf_group(nsurf, 2)
  integer, intent(in) ::  surf_indices_cp(nsurf, 2)
  integer, intent(in) ::  edge_indices_cp(nedge, 2)
  integer, intent(in) ::  str_indices_cp(nsurf, 2)
  integer, intent(in) ::  num_cp(ngroup)

  !Output
  double precision, intent(out) ::  Ca(nC)
  integer, intent(out) ::  Ci(nC), Cj(nC)

  !Working
  integer iC, isurf, iedge
  integer i_list(4), j_list(4), u_list(4), v_list(4), k
  integer i, j, u, v
  integer mu, mv

  iC = 0
  Ca(:) = 1.0
  Ci(:) = 0
  Cj(:) = 0

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
        iC = iC + 1
        Ci(iC) = str_indices_cp(isurf, 1) &
             + (v-1) * mu + u
        Cj(iC) = surf_ptrs(isurf, i, j)
     end do
  end do

  ! Loop through surfs and add edges
  do isurf = 1, nsurf
     mu = num_cp(surf_group(isurf, 1))
     mv = num_cp(surf_group(isurf, 2))

     iedge = surf_ptrs(isurf, 2, 1)
     v = 1
     do k = 1, mu-2
        if (iedge .gt. 0) then
           u = k + 1
        else
           u = mu - k
        end if
        iC = iC + 1
        Ci(iC) = str_indices_cp(isurf, 1) &
             + (v-1) * mu + u
        Cj(iC) = edge_indices_cp(abs(iedge), 1) + k
     end do

     iedge = surf_ptrs(isurf, 2, 3)
     v = mv
     do k = 1, mu-2
        if (iedge .gt. 0) then
           u = k + 1
        else
           u = mu - k
        end if
        iC = iC + 1
        Ci(iC) = str_indices_cp(isurf, 1) &
             + (v-1) * mu + u
        Cj(iC) = edge_indices_cp(abs(iedge), 1) + k
     end do

     iedge = surf_ptrs(isurf, 1, 2)
     u = 1
     do k = 1, mv-2
        if (iedge .gt. 0) then
           v = k + 1
        else
           v = mv - k
        end if
        iC = iC + 1
        Ci(iC) = str_indices_cp(isurf, 1) &
             + (v-1) * mu + u
        Cj(iC) = edge_indices_cp(abs(iedge), 1) + k
     end do

     iedge = surf_ptrs(isurf, 3, 2)
     u = mu
     do k = 1, mv-2
        if (iedge .gt. 0) then
           v = k + 1
        else
           v = mv - k
        end if
        iC = iC + 1
        Ci(iC) = str_indices_cp(isurf, 1) &
             + (v-1) * mu + u
        Cj(iC) = edge_indices_cp(abs(iedge), 1) + k
     end do
  end do

  ! Loop through surfaces and add interior cps
  do isurf = 1, nsurf
     mu = num_cp(surf_group(isurf, 1))
     mv = num_cp(surf_group(isurf, 2))
     do u = 1, mu-2
        do v = 1, mv-2
           iC = iC + 1
           Ci(iC) = str_indices_cp(isurf, 1) &
                + (v) * mu + u+1
           Cj(iC) = surf_indices_cp(isurf, 1) &
                + (v-1) * (mu-2) + u
        end do
     end do
  end do

  if (iC .ne. nC) then
     print *, 'Error in computeCmtx', iC, nC
     call exit(1)
  end if

end subroutine computeCmtx
