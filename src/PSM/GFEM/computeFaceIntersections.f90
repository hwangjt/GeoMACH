subroutine computeFaceIntersections(nmem, nint, faces, coords, intFaces, intCoords)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nmem, nint, faces, coords
  !f2py intent(out) intFaces, intCoords
  !f2py depend(nmem) faces, coords
  !f2py depend(nint) intFaces, intCoords

  !Input
  integer, intent(in) ::  nmem, nint, faces(nmem,4,2)
  double precision, intent(in) ::  coords(nmem,4,2,2,3)

  !Output
  integer, intent(out) ::  intFaces(nint,5)
  double precision, intent(out) ::  intCoords(nint,4)

  !Working
  integer imem, int, contr, intersectsFace

  int = 1
  do imem=1,nmem
     if (intersectsFace(imem,1,1,1,2,nmem,coords) .ne. 0) then
        contr = intersectsFace(imem,1,1,1,2,nmem,coords)
        intFaces(int,1:2) = faces(imem,contr,:)
        intFaces(int,3:5) = (/ imem , 2 , 1 /)
        intCoords(int,1:2) = coords(imem,contr,1,1,1:2)
        intCoords(int,3:4) = coords(imem,contr,1,2,1:2)
        int = int + 1
     end if
     if (intersectsFace(imem,2,1,2,2,nmem,coords) .ne. 0) then
        contr = intersectsFace(imem,2,1,2,2,nmem,coords)
        intFaces(int,1:2) = faces(imem,contr,:)
        intFaces(int,3:5) = (/ imem , 2 , 2 /)
        intCoords(int,1:2) = coords(imem,contr,2,1,1:2)
        intCoords(int,3:4) = coords(imem,contr,2,2,1:2)
        int = int + 1
     end if
     if (intersectsFace(imem,1,1,2,1,nmem,coords) .ne. 0) then
        contr = intersectsFace(imem,1,1,2,1,nmem,coords)
        intFaces(int,1:2) = faces(imem,contr,:)
        intFaces(int,3:5) = (/ imem , 1 , 1 /)
        intCoords(int,1:2) = coords(imem,contr,1,1,1:2)
        intCoords(int,3:4) = coords(imem,contr,2,1,1:2)
        int = int + 1
     end if
     if (intersectsFace(imem,1,2,2,2,nmem,coords) .ne. 0) then
        contr = intersectsFace(imem,1,2,2,2,nmem,coords)
        intFaces(int,1:2) = faces(imem,contr,:)
        intFaces(int,3:5) = (/ imem , 1 , 2 /)
        intCoords(int,1:2) = coords(imem,contr,1,2,1:2)
        intCoords(int,3:4) = coords(imem,contr,2,2,1:2)
        int = int + 1
     end if
  end do

end subroutine computeFaceIntersections




subroutine countFaceIntersections(nmem, coords, nint)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nmem, coords
  !f2py intent(out) nint
  !f2py depend(nmem) coords

  !Input
  integer, intent(in) ::  nmem
  double precision, intent(in) ::  coords(nmem,4,2,2,3)

  !Output
  integer, intent(out) ::  nint

  !Working
  integer imem, intersectsFace

  nint = 0
  do imem=1,nmem
     if (intersectsFace(imem,1,1,1,2,nmem,coords) .ne. 0) then
        nint = nint + 1
     end if
     if (intersectsFace(imem,2,1,2,2,nmem,coords) .ne. 0) then
        nint = nint + 1
     end if
     if (intersectsFace(imem,1,1,2,1,nmem,coords) .ne. 0) then
        nint = nint + 1
     end if
     if (intersectsFace(imem,1,2,2,2,nmem,coords) .ne. 0) then
        nint = nint + 1
     end if
  end do

end subroutine countFaceIntersections




function intersectsFace(imem, i1, j1, i2, j2, nmem, coords)

  implicit none

  !Input
  integer, intent(in) ::  imem, i1, j1, i2, j2, nmem
  double precision, intent(in) ::  coords(nmem,4,2,2,3)

  !Output
  integer intersectsFace

  !Working
  integer i
  double precision w1(4), w2(4)
  logical isOne

  w1 = coords(imem,:,i1,j1,3)
  w2 = coords(imem,:,i2,j2,3)

  intersectsFace = 0
  do i=1,4
     if (isOne(w1(i)) .and. isOne(sum(w1)) .and. isOne(w2(i)) .and. isOne(sum(w2))) then
        intersectsFace = i
     end if
  end do

end function intersectsFace




function isOne(val)

  implicit none

  !Input
  double precision, intent(in) ::  val

  !Output
  logical isOne

  isOne = abs(val - 1) .lt. 1e-12

end function isOne
