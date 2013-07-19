subroutine computeMemberTopology(nmem, membersInt, membersFlt, nvert, ngroup, mem_vert, mem_group)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nmem, membersInt, membersFlt
  !f2py intent(out) nvert, ngroup, mem_vert, mem_group
  !f2py depend(nmem) membersInt, membersFlt, mem_vert, mem_group

  !Input
  integer, intent(in) ::  nmem, membersInt(nmem,4,2)
  double precision, intent(in) ::  membersFlt(nmem,4,2,2,6)

  !Output
  integer, intent(out) ::  nvert, ngroup
  integer, intent(out) ::  mem_vert(nmem,2,2), mem_group(nmem,2,2)

  !Working
  integer imem1, imem2, i1, i2, j1, j2
  integer ivert1, ivert2, jvert1, jvert2
  integer ivert, igroup
  logical sameVert

  mem_vert(:,:,:) = 0
  ivert = 0
  do imem1=1,nmem
     do i1=1,2
        do j1=1,2
           if (mem_vert(imem1,i1,j1) .eq. 0) then
              ivert = ivert + 1
              mem_vert(imem1,i1,j1) = ivert
              do imem2=imem1+1,nmem
                 do i2=1,2
                    do j2=1,2
                       if (mem_vert(imem2,i2,j2) .eq. 0) then
                          if (sameVert(imem1,imem2,i1,j1,i2,j2,nmem,membersInt,membersFlt)) then
                             mem_vert(imem2,i2,j2) = ivert
                          end if
                       end if
                    end do
                 end do
              end do
           end if
        end do
     end do
  end do
  nvert = ivert

  mem_group(:,:,:) = 0
  igroup = 0
  do imem1=1,nmem
     do i1=1,2
        do j1=1,2
           if (mem_group(imem1,i1,j1) .eq. 0) then
              igroup = igroup + 1
              mem_group(imem1,i1,j1) = igroup
              do imem2=imem1+1,nmem
                 do i2=1,2
                    do j2=1,2
                       if (mem_group(imem2,i2,j2) .eq. 0) then
                          call getGroupVerts(imem1, i1, j1, nmem, mem_vert, ivert1, jvert1)
                          call getGroupVerts(imem2, i2, j2, nmem, mem_vert, ivert2, jvert2)
                          if ((ivert1 .eq. ivert2) .and. (jvert1 .eq. jvert2)) then
                             mem_group(imem2,i2,j2) = igroup
                          else if ((ivert1 .eq. jvert2) .and. (jvert1 .eq. ivert2)) then
                             mem_group(imem2,i2,j2) = -igroup
                          end if
                       end if
                    end do
                 end do
              end do
           end if
        end do
     end do
  end do
  ngroup = igroup

end subroutine computeMemberTopology




subroutine getGroupVerts(imem, i, j, nmem, mem_vert, ivert, jvert)

  implicit none

  !Input
  integer, intent(in) ::  imem, i, j, nmem, mem_vert(nmem,2,2)

  !Output
  integer, intent(out) ::  ivert, jvert

  if (i .eq. 1) then
     ivert = mem_vert(imem,1,j)
     jvert = mem_vert(imem,2,j)
  else
     ivert = mem_vert(imem,j,1)
     jvert = mem_vert(imem,j,2)
  end if

end subroutine getGroupVerts




function sameVert(imem1,imem2,i1,j1,i2,j2,nmem,membersInt,membersFlt)

  implicit none

  !Input
  integer, intent(in) ::  imem1, imem2, i1, j1, i2, j2, nmem, membersInt(nmem,4,2)
  double precision, intent(in) ::  membersFlt(nmem,4,2,2,6)

  !Output
  logical sameVert

  !Working
  double precision verts1(4,5), verts2(4,5)
  integer i, j, k1, k2, k3, k4
  logical zeroVec

  do i=1,4
     do j=1,2
        verts1(i,j) = dble(membersInt(imem1,i,j))
        verts2(i,j) = dble(membersInt(imem2,i,j))
     end do
  end do
  verts1(:,3:5) = membersFlt(imem1,:,i1,j1,1:3)
  verts2(:,3:5) = membersFlt(imem2,:,i2,j2,1:3)

  sameVert = .False.
  do k1=1,4
     do k2=1,4
        do k3=1,4
           do k4=1,4
              if ((k1+k2+k3+k4) .eq. 10) then
                 if ((k1*k2*k3*k4) .eq. 24) then
                    sameVert = sameVert .or. (&
                         zeroVec(verts1(1,:)-verts2(k1,:)) .and. &
                         zeroVec(verts1(2,:)-verts2(k2,:)) .and. &
                         zeroVec(verts1(3,:)-verts2(k3,:)) .and. &
                         zeroVec(verts1(4,:)-verts2(k4,:)))
                 end if
              end if
           end do
        end do
     end do
  end do

end function sameVert




function zeroVec(V)

  implicit none

  double precision, intent(in) ::  V(5)
  logical zeroVec
  
  if (dot_product(V,V) .lt. 1e-10) then
     zeroVec = .True.
  else
     zeroVec = .False.
  end if

end function zeroVec
