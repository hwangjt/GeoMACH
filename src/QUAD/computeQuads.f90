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
  integer iquad, ivert0, ivert1, ivert2, ivert3, ivert4, i1, i2, i3, i4

  iquad = 1
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
                    quads(iquad,1) = ivert1
                    quads(iquad,2) = ivert2
                    quads(iquad,3) = ivert3
                    quads(iquad,4) = ivert4
                    iquad = iquad + 1
                 end if
              end do
           end do
        end do
     end do
  end do

end subroutine computeQuads
