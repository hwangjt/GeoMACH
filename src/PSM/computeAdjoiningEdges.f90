subroutine computeAdjoiningEdges(nmem, nadj, membersInt, membersFlt, adjoiningInt, adjoiningFlt)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nmem, nadj, membersInt, membersFlt
  !f2py intent(out) adjoiningInt, adjoiningFlt
  !f2py depend(nmem) membersInt, membersFlt
  !f2py depend(nadj) adjoiningInt, adjoiningFlt

  !Input
  integer, intent(in) ::  nmem, nadj, membersInt(nmem,4,2)
  double precision, intent(in) ::  membersFlt(nmem,4,2,2,6)

  !Output
  integer, intent(out) ::  adjoiningInt(nadj,5)
  double precision, intent(out) ::  adjoiningFlt(nadj,4)

  !Working
  integer imem, iadj, isrc, intersectsFace

  iadj = 1
  do imem=1,nmem
     isrc = intersectsFace(imem,1,1,2,1,nmem,membersFlt)
     if (isrc .ne. 0) then
        adjoiningInt(iadj,1:2) = membersInt(imem,isrc,:)
        adjoiningInt(iadj,3:5) = (/ imem , 1 , 1 /)
        adjoiningFlt(iadj,1:2) = membersFlt(imem,isrc,1,1,2:3)
        adjoiningFlt(iadj,3:4) = membersFlt(imem,isrc,2,1,2:3)
        iadj = iadj + 1
     end if
     isrc = intersectsFace(imem,1,2,2,2,nmem,membersFlt)
     if (isrc .ne. 0) then
        adjoiningInt(iadj,1:2) = membersInt(imem,isrc,:)
        adjoiningInt(iadj,3:5) = (/ imem , 1 , 2 /)
        adjoiningFlt(iadj,1:2) = membersFlt(imem,isrc,1,2,2:3)
        adjoiningFlt(iadj,3:4) = membersFlt(imem,isrc,2,2,2:3)
        iadj = iadj + 1
     end if
     isrc = intersectsFace(imem,1,1,1,2,nmem,membersFlt)
     if (isrc .ne. 0) then
        adjoiningInt(iadj,1:2) = membersInt(imem,isrc,:)
        adjoiningInt(iadj,3:5) = (/ imem , 2 , 1 /)
        adjoiningFlt(iadj,1:2) = membersFlt(imem,isrc,1,1,2:3)
        adjoiningFlt(iadj,3:4) = membersFlt(imem,isrc,1,2,2:3)
        iadj = iadj + 1
     end if
     isrc = intersectsFace(imem,2,1,2,2,nmem,membersFlt)
     if (isrc .ne. 0) then
        adjoiningInt(iadj,1:2) = membersInt(imem,isrc,:)
        adjoiningInt(iadj,3:5) = (/ imem , 2 , 2 /)
        adjoiningFlt(iadj,1:2) = membersFlt(imem,isrc,2,1,2:3)
        adjoiningFlt(iadj,3:4) = membersFlt(imem,isrc,2,2,2:3)
        iadj = iadj + 1
     end if
  end do

end subroutine computeAdjoiningEdges




subroutine countAdjoiningEdges(nmem, membersFlt, nadj)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nmem, membersFlt
  !f2py intent(out) nadj
  !f2py depend(nmem) membersFlt

  !Input
  integer, intent(in) ::  nmem
  double precision, intent(in) ::  membersFlt(nmem,4,2,2,6)

  !Output
  integer, intent(out) ::  nadj

  !Working
  integer imem, intersectsFace

  nadj = 0
  do imem=1,nmem
     if (intersectsFace(imem,1,1,2,1,nmem,membersFlt) .ne. 0) then
        nadj = nadj + 1
     end if
     if (intersectsFace(imem,1,2,2,2,nmem,membersFlt) .ne. 0) then
        nadj = nadj + 1
     end if
     if (intersectsFace(imem,1,1,1,2,nmem,membersFlt) .ne. 0) then
        nadj = nadj + 1
     end if
     if (intersectsFace(imem,2,1,2,2,nmem,membersFlt) .ne. 0) then
        nadj = nadj + 1
     end if
  end do

end subroutine countAdjoiningEdges




function intersectsFace(imem, i1, j1, i2, j2, nmem, membersFlt)

  implicit none

  !Input
  integer, intent(in) ::  imem, i1, j1, i2, j2, nmem
  double precision, intent(in) ::  membersFlt(nmem,4,2,2,6)

  !Output
  integer intersectsFace

  !Working
  integer i
  double precision w1(4), w2(4)
  logical isOne

  w1 = membersFlt(imem,:,i1,j1,1)
  w2 = membersFlt(imem,:,i2,j2,1)

  intersectsFace = 0
  do i=1,4
     if (isOne(w1(i)) .and. isOne(sum(w1)) .and. &
          isOne(w2(i)) .and. isOne(sum(w2))) then
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
