subroutine countPreviewMembers(isurf, isrc, nmem, membersFlt, npts)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) isurf, isrc, nmem, membersFlt
  !f2py intent(out) npts
  !f2py depend(nmem) membersFlt

  !Input
  integer, intent(in) ::  isurf, isrc, nmem
  double precision, intent(in) ::  membersFlt(nmem,4,2,2,6)

  !Output
  integer, intent(out) ::  npts

  !Working
  integer imem, ic, jc

  npts = 0
  do imem=1,nmem
     do ic=1,2
        do jc=1,2
           if (int(membersFlt(imem,isrc,ic,jc,4)) .eq. isurf) then
              npts = npts + 1
           end if
        end do
     end do
  end do

end subroutine countPreviewMembers




subroutine computePreviewMemberProj(isurf, isrc, nmem, npts, membersFlt, inds, P, Q)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) isurf, isrc, nmem, npts, membersFlt
  !f2py intent(out) inds, P, Q
  !f2py depend(nmem) membersFlt
  !f2py depend(npts) inds, P, Q

  !Input
  integer, intent(in) ::  isurf, isrc, nmem, npts
  double precision, intent(in) ::  membersFlt(nmem,4,2,2,6)

  !Output
  integer, intent(out) ::  inds(npts)
  double precision, intent(out) ::  P(npts,3), Q(npts,3)

  !Working
  integer imem, ic, jc, ipts

  P(:,:) = 0.
  Q(:,:) = 0.
  Q(:,3) = 1.

  ipts = 1
  do imem=1,nmem
     do ic=1,2
        do jc=1,2
           if (int(membersFlt(imem,isrc,ic,jc,4)) .eq. isurf) then
              inds(ipts) = 4*imem - 4 + 2*(jc-1) + ic
              P(ipts,1) = membersFlt(imem,isrc,ic,jc,5)
              P(ipts,2) = membersFlt(imem,isrc,ic,jc,6)
              ipts = ipts + 1
           end if
        end do
     end do
  end do

end subroutine computePreviewMemberProj




subroutine computePreviewMemberWeights(nmem, nnode, membersFlt, quads, W)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nmem, nnode, membersFlt
  !f2py intent(out) quads, W
  !f2py depend(nmem) membersFlt, quads
  !f2py depend(nnode) W

  !Input
  integer, intent(in) ::  nmem, nnode
  double precision, intent(in) ::  membersFlt(nmem,4,2,2,6)

  !Output
  integer, intent(out) ::  quads(nmem,4)
  double precision, intent(out) ::  W(nnode,4)
  
  !Working
  integer imem, isrc, index, ic, jc

  do imem=1,nmem
     index = 4*imem - 4
     quads(imem,:) = (/ index+1, index+2, index+4, index+3 /)
     do isrc=1,4
        do ic=1,2
           do jc=1,2
              index = 4*imem - 4 + 2*(jc-1) + ic
              W(index,isrc) = membersFlt(imem,isrc,ic,jc,1)
              if (W(index,isrc) .lt. 0) then
                 W(index,isrc) = 0.
              end if
           end do
        end do
     end do
  end do

end subroutine computePreviewMemberWeights
