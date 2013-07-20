subroutine importMembers(nmem, members, membersInt, membersFlt)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nmem, members
  !f2py intent(out) membersInt, membersFlt
  !f2py depend(nmem) members, membersInt, membersFlt

  !Input
  integer, intent(in) ::  nmem
  double precision, intent(in) ::  members(nmem,4,5,3)

  !Output
  integer, intent(out) ::  membersInt(nmem,4,2)
  double precision, intent(out) ::  membersFlt(nmem,4,2,2,6)

  !Working
  integer imem, isrc

  membersInt(:,:,:) = 0
  membersFlt(:,:,:,:,:) = 0.

  do imem=1,nmem
     do isrc=1,4
        if (members(imem,isrc,1,1) .gt. 0) then
           membersInt(imem,isrc,1) = int(members(imem,isrc,1,1)) + 1
           membersInt(imem,isrc,2) = int(members(imem,isrc,1,2)) + 1
           membersFlt(imem,isrc,1,1,1:3) = members(imem,isrc,2,:)
           membersFlt(imem,isrc,2,1,1:3) = members(imem,isrc,3,:)
           membersFlt(imem,isrc,1,2,1:3) = members(imem,isrc,4,:)
           membersFlt(imem,isrc,2,2,1:3) = members(imem,isrc,5,:)
        end if
     end do
  end do

end subroutine importMembers




subroutine computePreviewMemberCoords(comp, face, ni, nj, nmem, &
     surfs, idims, jdims, membersInt, membersFlt)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) comp, face, ni, nj, nmem, surfs, idims, jdims
  !f2py intent(inout) membersInt, membersFlt
  !f2py depend(ni,nj) surfs
  !f2py depend(ni) idims
  !f2py depend(nj) jdims
  !f2py depend(nmem) membersInt, membersFlt

  !Input
  integer, intent(in) ::  comp, face, ni, nj, nmem, surfs(ni,nj)
  double precision, intent(in) ::  idims(ni+1), jdims(nj+1)

  !Output
  integer, intent(inout) ::  membersInt(nmem,4,2)
  double precision, intent(inout) ::  membersFlt(nmem,4,2,2,6)

  !Working
  integer imem, isrc, ic, jc, i, j, ifound, jfound
  double precision a, b, u, v

  do imem=1,nmem
     do isrc=1,4
        if ((membersInt(imem,isrc,1) .eq. comp) .and. &
             (membersInt(imem,isrc,2) .eq. face)) then
           do ic=1,2
              do jc=1,2
                 a = membersFlt(imem,isrc,ic,jc,2)
                 b = membersFlt(imem,isrc,ic,jc,3)
                 iloop: do i=1,ni
                    if (a .le. idims(i+1)) then
                       ifound = i
                       u = (a-idims(i))/(idims(i+1)-idims(i))
                       exit iloop
                    end if
                 end do iloop
                 jloop: do j=1,nj
                    if (b .le. jdims(j+1)) then
                       jfound = j
                       v = (b-jdims(j))/(jdims(j+1)-jdims(j))
                       exit jloop
                    end if
                 end do jloop
                 membersFlt(imem,isrc,ic,jc,4) = surfs(ifound,jfound)
                 membersFlt(imem,isrc,ic,jc,5) = u
                 membersFlt(imem,isrc,ic,jc,6) = v
              end do
           end do
        end if
     end do
  end do

end subroutine computePreviewMemberCoords
