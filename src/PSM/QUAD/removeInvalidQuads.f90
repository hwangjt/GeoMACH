subroutine removeInvalidQuads(nquad0, nquad, quads0, quads)

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
  integer iquad0, iquad

  iquad = 0
  do iquad0=1,nquad0
     if (&
          (quads0(iquad0,1) .ne. quads0(iquad0,2)) .and. &
          (quads0(iquad0,1) .ne. quads0(iquad0,3)) .and. &
          (quads0(iquad0,1) .ne. quads0(iquad0,4)) .and. &
          (quads0(iquad0,2) .ne. quads0(iquad0,3)) .and. &
          (quads0(iquad0,2) .ne. quads0(iquad0,4)) .and. &
          (quads0(iquad0,3) .ne. quads0(iquad0,4))) then
        iquad = iquad + 1
        quads(iquad,:) = quads0(iquad0,:)
     end if
  end do

end subroutine removeInvalidQuads



subroutine countInvalidQuads(nquad, quads, ninv)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nquad, quads
  !f2py intent(out) ninv
  !f2py depend(nquad) quads

  !Input
  integer, intent(in) ::  nquad, quads(nquad,4)

  !Output
  integer, intent(out) ::  ninv

  !Working
  integer iquad

  ninv = 0
  do iquad=1,nquad
     if (&
          (quads(iquad,1) .eq. quads(iquad,2)) .or. &
          (quads(iquad,1) .eq. quads(iquad,3)) .or. &
          (quads(iquad,1) .eq. quads(iquad,4)) .or. &
          (quads(iquad,2) .eq. quads(iquad,3)) .or. &
          (quads(iquad,2) .eq. quads(iquad,4)) .or. &
          (quads(iquad,3) .eq. quads(iquad,4))) then
        ninv = ninv + 1
     end if
  end do

end subroutine countInvalidQuads
