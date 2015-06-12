subroutine removeDuplicateQuads(nquad0, nquad, quads0, quads)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nquad0, nquad, quads0
  !f2py intent(out) quads
  !f2py depend(nquad0) quads0
  !f2py depend(nquad) quads

  !Input
  integer, intent(in) ::  nquad0, nquad, quads0(nquad0,4)

  !Output
  integer, intent(out) ::  quads(nquad,4)

  !Working
  logical sameQuad, unused(nquad0)
  integer i1, i2, iquad

  unused(:) = .True.

  iquad = 0
  do i1=1,nquad0
     if (unused(i1)) then
        iquad = iquad + 1
        quads(iquad,:) = quads0(i1,:)
        do i2=i1+1,nquad0
           if (unused(i2)) then
              if (sameQuad(quads0(i1,:),quads0(i2,:))) then
                 unused(i2) = .False.
              end if
           end if
        end do
     end if
  end do
  if (iquad .ne. nquad) then
     print *, 'Error in deleteDuplicateQuads', iquad, nquad
  end if

end subroutine removeDuplicateQuads




function sameQuad(A, B)

  implicit none

  !Input
  integer, intent(in) ::  A(4), B(4)

  !Output
  logical sameQuad

  if ((A(1) .eq. B(1)) .and. (A(3) .eq. B(3))) then
     sameQuad = .True.
  else if ((A(1) .eq. B(2)) .and. (A(3) .eq. B(4))) then
     sameQuad = .True.
  else if ((A(1) .eq. B(3)) .and. (A(3) .eq. B(1))) then
     sameQuad = .True.
  else if ((A(1) .eq. B(4)) .and. (A(3) .eq. B(2))) then
     sameQuad = .True.
  else
     sameQuad = .False.
  end if

end function sameQuad
