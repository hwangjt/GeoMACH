subroutine computeQuads(nvert, nadj, nquad, adjPtr, adjMap, quads)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvert, nadj, nquad, adjPtr, adjMap
  !f2py intent(out) quads
  !f2py depend(nvert) adjPtr
  !f2py depend(nadj) adjMap
  !f2py depend(nquad) quads

  !Input
  integer, intent(in) ::  nvert, nadj, nquad, adjPtr(nvert,2), adjMap(nadj)

  !Output
  integer, intent(out) ::  quads(nquad,4)

  !Working
  integer iquad, ivert0, ivert1, ivert2, ivert3, ivert4, i, i1, i2, i3, i4
  logical found

  iquad = 0
  do ivert0=1,nvert
     do i1=adjPtr(ivert0,1),adjPtr(ivert0,2)
        ivert1 = adjMap(i1)
        do i2=adjPtr(ivert1,1),adjPtr(ivert1,2)
           ivert2 = adjMap(i2)
           do i3=adjPtr(ivert2,1),adjPtr(ivert2,2)
              ivert3 = adjMap(i3)
              do i4=adjPtr(ivert3,1),adjPtr(ivert3,2)
                 ivert4 = adjMap(i4)
                 if ((ivert0 .eq. ivert4) .and. (ivert1 .ne. ivert3) &
                      .and. (ivert0 .ne. ivert2)) then
                    found = .False.
                    do i=adjPtr(ivert0,1),adjPtr(ivert0,2)
                       if (ivert2 .eq. adjMap(i)) then
                          found = .True.
                       end if
                    end do
                    do i=adjPtr(ivert1,1),adjPtr(ivert1,2)
                       if (ivert3 .eq. adjMap(i)) then
                          found = .True.
                       end if
                    end do
                    if (.not. found) then
                       iquad = iquad + 1
                       quads(iquad,1) = ivert1
                       quads(iquad,2) = ivert2
                       quads(iquad,3) = ivert3
                       quads(iquad,4) = ivert4
                    end if
                 end if
              end do
           end do
        end do
     end do
  end do

end subroutine computeQuads




subroutine countQuads(nvert, nadj, adjPtr, adjMap, nquad)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvert, nadj, adjPtr, adjMap
  !f2py intent(out) nquad
  !f2py depend(nvert) adjPtr
  !f2py depend(nadj) adjMap
  
  !Input
  integer, intent(in) ::  nvert, nadj, adjPtr(nvert,2), adjMap(nadj)

  !Output
  integer, intent(out) ::  nquad

  !Working
  integer ivert0, ivert1, ivert2, ivert3, ivert4, i1, i2, i3, i4, i
  logical found

  nquad = 0
  do ivert0=1,nvert
     do i1=adjPtr(ivert0,1),adjPtr(ivert0,2)
        ivert1 = adjMap(i1)
        do i2=adjPtr(ivert1,1),adjPtr(ivert1,2)
           ivert2 = adjMap(i2)
           do i3=adjPtr(ivert2,1),adjPtr(ivert2,2)
              ivert3 = adjMap(i3)
              do i4=adjPtr(ivert3,1),adjPtr(ivert3,2)
                 ivert4 = adjMap(i4)
                 if ((ivert0 .eq. ivert4) .and. (ivert1 .ne. ivert3) &
                      .and. (ivert0 .ne. ivert2)) then
                    found = .False.
                    do i=adjPtr(ivert0,1),adjPtr(ivert0,2)
                       if (ivert2 .eq. adjMap(i)) then
                          found = .True.
                       end if
                    end do
                    do i=adjPtr(ivert1,1),adjPtr(ivert1,2)
                       if (ivert3 .eq. adjMap(i)) then
                          found = .True.
                       end if
                    end do
                    if (.not. found) then
                       nquad = nquad + 1
                    end if
                 end if
              end do
           end do
        end do
     end do
  end do

end subroutine countQuads




subroutine rotateQuads(nvert, nquad, verts, quads)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvert, nquad, verts
  !f2py intent(inout) quads
  !f2py depend(nvert) verts
  !f2py depend(nquad) quads

  !Input
  integer, intent(in) ::  nvert, nquad
  double precision, intent(in) ::  verts(nvert,2)

  !Output
  integer, intent(inout) ::  quads(nquad,4)

  !Working
  double precision v1(2), v2(2)
  integer iquad

  do iquad=1,nquad
     v1 = verts(quads(iquad,2),:) - verts(quads(iquad,1),:)
     v2 = verts(quads(iquad,3),:) - verts(quads(iquad,2),:)
     if ((v1(1)*v2(2) - v1(2)*v2(1)) .gt. 0) then
        quads(iquad,:) = quads(iquad,4:1:-1)
     end if
  end do

end subroutine rotateQuads
