subroutine computeTriangles(nvert, nadj, ntri, adjPtr, adjMap, triangles)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvert, nadj, ntri, adjPtr, adjMap
  !f2py intent(out) triangles
  !f2py depend(nvert) adjPtr
  !f2py depend(nadj) adjMap
  !f2py depend(ntri) triangles

  !Input
  integer, intent(in) ::  nvert, nadj, ntri, adjPtr(nvert,2), adjMap(nadj)

  !Output
  integer, intent(out) ::  triangles(ntri,3)

  !Working
  integer itri, ivert0, ivert1, ivert2, ivert3, i1, i2, i3

  itri = 1
  do ivert0=1,nvert
     do i1=adjPtr(ivert0,1),adjPtr(ivert0,2)
        ivert1 = adjMap(i1)
        do i2=adjPtr(ivert1,1),adjPtr(ivert1,2)
           ivert2 = adjMap(i2)
           do i3=adjPtr(ivert2,1),adjPtr(ivert2,2)
              ivert3 = adjMap(i3)
              if (ivert0 .eq. ivert3) then
                 triangles(itri,1) = ivert1
                 triangles(itri,2) = ivert2
                 triangles(itri,3) = ivert3
                 itri = itri + 1
              end if
           end do
        end do
     end do
  end do

end subroutine computeTriangles
