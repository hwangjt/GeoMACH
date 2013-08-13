subroutine computeGroupSplits(nsurf, nmem, ngroup, nint, nsplit, maxL, &
     surf_group, mem_group, surfEdgeLengths, memEdgeLengths, &
     groupIntPtr, groupInts, groupSplitPtr, groupSplits)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nsurf, nmem, ngroup, nint, nsplit, maxL, surf_group, mem_group, surfEdgeLengths, memEdgeLengths, groupIntPtr, groupInts, groupSplitPtr
  !f2py intent(out) groupSplits
  !f2py depend(nsurf) surf_group, surfEdgeLengths
  !f2py depend(nmem) mem_group, memEdgeLengths
  !f2py depend(ngroup) groupIntPtr, groupSplitPtr
  !f2py depend(nint) groupInts
  !f2py depend(nsplit) groupSplits

  !Input
  integer, intent(in) ::  nsurf, nmem, ngroup, nint, nsplit
  double precision, intent(in) ::  maxL
  integer, intent(in) ::  surf_group(nsurf,2,2), mem_group(nmem,2,2)
  double precision, intent(in) ::  surfEdgeLengths(nsurf,2,2)
  double precision, intent(in) ::  memEdgeLengths(nmem,2,2)
  integer, intent(in) ::  groupIntPtr(ngroup,2)
  double precision, intent(in) ::  groupInts(nint)
  integer, intent(in) ::  groupSplitPtr(ngroup,2)

  !Output
  double precision, intent(out) ::  groupSplits(nsplit)

  !Working
  integer isurf, imem, igroup, i, j, n, k, i1, i2, m, l, counter(ngroup), index
  double precision len, den
  logical done(ngroup)
  double precision, allocatable, dimension(:) ::  t

  done = .False.
  counter(:) = 0
  do isurf=1,nsurf
     do i=1,2
        do j=1,2
           igroup = abs(surf_group(isurf,i,j))
           if (.not. done(igroup)) then
              done(igroup) = .True.
              len = surfEdgeLengths(isurf,i,j)
              i1 = groupIntPtr(igroup,1)
              i2 = groupIntPtr(igroup,2)
              n = i2 - i1 + 1
              allocate(t(n+2))
              t(1) = 0.
              if (n .gt. 0) then
                 t(2:n+1) = groupInts(i1:i2)
              end if
              t(n+2) = 1.
              do k=1,n+1
                 m = floor(len*(t(k+1)-t(k))/maxL)
                 den = 1.0/(m+1)
                 do l=1,m
                    index = groupSplitPtr(igroup,1) + counter(igroup)
                    groupSplits(index) = t(k) + (t(k+1)-t(k))*l*den
                    counter(igroup) = counter(igroup) + 1
                 end do
              end do
              deallocate(t)
           end if
        end do
     end do
  end do
  do imem=1,nmem
     do i=1,2
        do j=1,2
           igroup = abs(mem_group(imem,i,j))
           if (.not. done(igroup)) then
              done(igroup) = .True.
              len = memEdgeLengths(imem,i,j)
              i1 = groupIntPtr(igroup,1)
              i2 = groupIntPtr(igroup,2)
              n = i2 - i1 + 1
              allocate(t(n+2))
              t(1) = 0.
              if (n .gt. 0) then
                 t(2:n+1) = groupInts(i1:i2)
              end if
              t(n+2) = 1.
              do k=1,n+1
                 m = floor(len*(t(k+1)-t(k))/maxL)
                 den = 1.0/(m+1)
                 do l=1,m
                    index = groupSplitPtr(igroup,1) + counter(igroup)
                    groupSplits(index) = t(k) + (t(k+1)-t(k))*l*den
                    counter(igroup) = counter(igroup) + 1
                 end do
              end do
              deallocate(t)
           end if
        end do
     end do
  end do

end subroutine computeGroupSplits




subroutine countGroupSplits(nsurf, nmem, ngroup, nint, maxL, &
     surf_group, mem_group, surfEdgeLengths, memEdgeLengths, &
     groupIntPtr, groupInts, groupSplitCount)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nsurf, nmem, ngroup, nint, maxL, surf_group, mem_group, surfEdgeLengths, memEdgeLengths, groupIntPtr, groupInts
  !f2py intent(out) groupSplitCount
  !f2py depend(nsurf) surf_group, surfEdgeLengths
  !f2py depend(nmem) mem_group, memEdgeLengths
  !f2py depend(ngroup) groupIntPtr, groupSplitCount
  !f2py depend(nint) groupInts

  !Input
  integer, intent(in) ::  nsurf, nmem, ngroup, nint
  double precision, intent(in) ::  maxL
  integer, intent(in) ::  surf_group(nsurf,2,2), mem_group(nmem,2,2)
  double precision, intent(in) ::  surfEdgeLengths(nsurf,2,2)
  double precision, intent(in) ::  memEdgeLengths(nmem,2,2)
  integer, intent(in) ::  groupIntPtr(ngroup,2)
  double precision, intent(in) ::  groupInts(nint)

  !Output
  integer, intent(out) ::  groupSplitCount(ngroup)

  !Working
  integer isurf, imem, igroup, i, j, n, k, i1, i2
  double precision len
  logical done(ngroup)
  double precision, allocatable, dimension(:) ::  t

  done = .False.
  groupSplitCount = 0
  do isurf=1,nsurf
     do i=1,2
        do j=1,2
           igroup = abs(surf_group(isurf,i,j))
           if (.not. done(igroup)) then
              done(igroup) = .True.
              len = surfEdgeLengths(isurf,i,j)
              i1 = groupIntPtr(igroup,1)
              i2 = groupIntPtr(igroup,2)
              n = i2 - i1 + 1
              allocate(t(n+2))
              t(1) = 0.
              if (n .gt. 0) then
                 t(2:n+1) = groupInts(i1:i2)
              end if
              t(n+2) = 1.
              do k=1,n+1
                 groupSplitCount(igroup) = groupSplitCount(igroup) + &
                      floor(len*(t(k+1)-t(k))/maxL)
              end do
              deallocate(t)
           end if
        end do
     end do
  end do
  do imem=1,nmem
     do i=1,2
        do j=1,2
           igroup = abs(mem_group(imem,i,j))
           if (.not. done(igroup)) then
              done(igroup) = .True.
              len = memEdgeLengths(imem,i,j)
              i1 = groupIntPtr(igroup,1)
              i2 = groupIntPtr(igroup,2)
              n = i2 - i1 + 1
              allocate(t(n+2))
              t(1) = 0.
              if (n .gt. 0) then
                 t(2:n+1) = groupInts(i1:i2)
              end if
              t(n+2) = 1.
              do k=1,n+1
                 groupSplitCount(igroup) = groupSplitCount(igroup) + &
                      floor(len*(t(k+1)-t(k))/maxL)
              end do
              deallocate(t)
           end if
        end do
     end do
  end do

end subroutine countGroupSplits
