subroutine deleteDuplicateTriangles(ntri0, ntri, triangles0, triangles)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) ntri0, ntri, triangles0
  !f2py intent(out) triangles
  !f2py depend(ntri0) triangles0
  !f2py depend(ntri) triangles

  !Input
  integer, intent(in) ::  ntri0, ntri, triangles0(ntri0,3)

  !Output
  integer, intent(out) ::  triangles(ntri,3)

  !Working
  logical sameTri, unused(ntri0)
  integer i1, i2, itri

  unused(:) = .True.

  itri = 1
  do i1=1,ntri0
     if (unused(i1)) then
        triangles(itri,:) = triangles0(i1,:)
        itri = itri + 1
        do i2=i1+1,ntri0
           if (unused(i2)) then
              if (sameTri(triangles0(i1,:),triangles0(i2,:))) then
                 unused(i2) = .False.
              end if
           end if
        end do
     end if
  end do

end subroutine deleteDuplicateTriangles




function sameTri(A, B)

  implicit none

  !Input
  integer, intent(in) ::  A(3), B(3)

  !Output
  logical sameTri

  !Working
  integer minA, maxA, minB, maxB

  minA = minloc(A,1)
  maxA = maxloc(A,1)
  minB = minloc(B,1)
  maxB = maxloc(B,1)

  if ((A(minA) .eq. B(minB)) .and. &
       (A(6-minA-maxA) .eq. B(6-minB-maxB)) &
       .and. (A(maxA) .eq. B(maxB))) then
     sameTri = .True.
  else
     sameTri = .False.
  end if

end function sameTri
