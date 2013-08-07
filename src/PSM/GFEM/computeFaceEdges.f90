subroutine computeFaceEdges(comp, face, ni, nj, nsurf, nmem, nadj, nedge, &
     idims, jdims, face_surf, surf_group, mem_group, adjoiningInt, adjoiningFlt, &
     surfEdgeLengths, memEdgeLengths, edge_group, edgeLengths, edges)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) comp, face, ni, nj, nsurf, nmem, nadj, nedge, idims, jdims, face_surf, surf_group, mem_group, adjoiningInt, adjoiningFlt, surfEdgeLengths, memEdgeLengths
  !f2py intent(out) edge_group, edgeLengths, edges
  !f2py depend(ni) idims
  !f2py depend(nj) jdims
  !f2py depend(ni,nj) face_surf
  !f2py depend(nsurf) surf_group
  !f2py depend(nmem) mem_group
  !f2py depend(nadj) adjoiningInt, adjoiningFlt
  !f2py depend(nsurf) surfEdgeLengths
  !f2py depend(nmem) memEdgeLengths
  !f2py depend(nedge) edge_group, edgeLengths, edges

  !Input
  integer, intent(in) ::  comp, face, ni, nj, nsurf, nmem, nadj, nedge
  double precision, intent(in) ::  idims(ni+1), jdims(nj+1)
  integer, intent(in) ::  face_surf(ni,nj), surf_group(nsurf,2,2)
  integer, intent(in) ::  mem_group(nmem,2,2), adjoiningInt(nadj,5)
  double precision, intent(in) ::  adjoiningFlt(nadj,4)
  double precision, intent(in) ::  surfEdgeLengths(nsurf,2,2)
  double precision, intent(in) ::  memEdgeLengths(nmem,2,2)

  !Output
  integer, intent(out) ::  edge_group(nedge)
  double precision, intent(out) ::  edgeLengths(nedge), edges(nedge,2,2)

  !Working
  double precision verts(ni+1,nj+1,2)
  integer iadj, iedge, isurf, i, j, imem

  do i=1,ni+1
     do j=1,nj+1
        verts(i,j,1) = idims(i)
        verts(i,j,2) = jdims(j)
     end do
  end do

  do i=1,ni
     do j=1,nj
        isurf = face_surf(i,j)

        iedge = (j-1)*ni + i
        edges(iedge,1,:) = verts(i,j,:)        
        edges(iedge,2,:) = verts(i+1,j,:)
        edge_group(iedge) = surf_group(isurf,1,1)
        edgeLengths(iedge) = surfEdgeLengths(isurf,1,1)

        iedge = j*ni + i
        edges(iedge,1,:) = verts(i,j+1,:)        
        edges(iedge,2,:) = verts(i+1,j+1,:)
        edge_group(iedge) = surf_group(isurf,1,2)
        edgeLengths(iedge) = surfEdgeLengths(isurf,1,2)

        iedge = (nj+1)*ni + (i-1)*nj + j
        edges(iedge,1,:) = verts(i,j,:)
        edges(iedge,2,:) = verts(i,j+1,:)
        edge_group(iedge) = surf_group(isurf,2,1)
        edgeLengths(iedge) = surfEdgeLengths(isurf,2,1)

        iedge = (nj+1)*ni + i*nj + j
        edges(iedge,1,:) = verts(i+1,j,:)
        edges(iedge,2,:) = verts(i+1,j+1,:)
        edge_group(iedge) = surf_group(isurf,2,2)
        edgeLengths(iedge) = surfEdgeLengths(isurf,2,2)
     end do
  end do

  iedge = (nj+1)*ni + (ni+1)*nj + 1

  do iadj=1,nadj
     if ((adjoiningInt(iadj,1) .eq. comp) .and. (adjoiningInt(iadj,2) .eq. face)) then
        imem = adjoiningInt(iadj,3)
        i = adjoiningInt(iadj,4)
        j = adjoiningInt(iadj,5)
        edges(iedge,1,:) = adjoiningFlt(iadj,1:2)
        edges(iedge,2,:) = adjoiningFlt(iadj,3:4)
        edge_group(iedge) = mem_group(imem,i,j)
        edgeLengths(iedge) = memEdgeLengths(imem,i,j)
        iedge = iedge + 1
     end if
  end do

end subroutine computeFaceEdges



subroutine countFaceEdges(comp, face, ni, nj, nadj, adjoiningInt, nedge)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) comp, face, ni, nj, nadj, adjoiningInt
  !f2py intent(out) nedge
  !f2py depend(nadj) adjoiningInt

  !Input
  integer, intent(in) ::  comp, face, ni, nj, nadj
  integer, intent(in) ::  adjoiningInt(nadj,5)

  !Output
  integer, intent(out) ::  nedge

  !Working
  integer iadj

  nedge = ni*(nj+1) + nj*(ni+1)

  do iadj=1,nadj
     if ((adjoiningInt(iadj,1) .eq. comp) .and. &
          (adjoiningInt(iadj,2) .eq. face)) then
        nedge = nedge + 1
     end if
  end do

end subroutine countFaceEdges
