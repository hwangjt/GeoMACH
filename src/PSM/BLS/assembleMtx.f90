subroutine assembleRHS1(nvert, nB, verts, vertCon, B)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvert, nB, verts, vertCon
  !f2py intent(out) B
  !f2py depend(nvert) verts, vertCon
  !f2py depend(nB) B

  !Input
  integer, intent(in) ::  nvert, nB
  double precision, intent(in) ::  verts(nvert,2)
  logical, intent(in) ::  vertCon(nvert)

  !Output
  double precision, intent(out) ::  B(nB,2)

  !Working
  integer iB, ivert

  B(:,:) = 0.0
  iB = nB
  do ivert=nvert,1,-1
     if (vertCon(ivert)) then
        B(iB,:) = verts(ivert,:)
        iB = iB - 1
     end if
  end do

end subroutine assembleRHS1




subroutine assembleConMtx1(nvert, nC, vertCon, Ca, Ci, Cj)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvert, nC, vertCon
  !f2py intent(out) Ca, Ci, Cj
  !f2py depend(nvert) vertCon
  !f2py depend(nC) Ca, Ci, Cj

  !Input
  integer, intent(in) ::  nvert, nC
  logical, intent(in) ::  vertCon(nvert)

  !Output
  double precision, intent(out) ::  Ca(nC)
  integer, intent(out) ::  Ci(nC), Cj(nC)

  !Working
  integer iC, ivert, n

  Ca(:) = 1.0

  iC = 0
  do ivert=1,nvert
     if (vertCon(ivert)) then
        iC = iC + 1
        Ci(iC) = nvert + iC
        Cj(iC) = ivert
     end if
  end do
  if (2*iC .ne. nC) then
     print *, 'Error in assembleConMtx1', 2*iC, nC
     call exit(1)
  end if

  n = int(nC/2)

  Cj(n+1:nC) = Ci(1:n)
  Ci(n+1:nC) = Cj(1:n)

end subroutine assembleConMtx1




subroutine assembleMtx1(nquad, nM, quads, Ma, Mi, Mj)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nquad, nM, quads
  !f2py intent(out) Ma, Mi, Mj
  !f2py depend(nquad) quads
  !f2py depend(nM) Ma, Mi, Mj

  !Input
  integer, intent(in) ::  nquad, nM, quads(nquad,4)

  !Output
  double precision, intent(out) ::  Ma(nM)
  integer, intent(out) ::  Mi(nM), Mj(nM)

  !Working
  integer getDOFindex1, iquad, iM, i, j, k, l
  double precision F(2,2), B(2,2)

  F(1,:) = (/  1. , -1. /)
  F(2,:) = (/ -1. ,  1. /)

  B(1,:) = (/  2. ,  1. /)
  B(2,:) = (/  1. ,  2. /)
  
  iM = 0
  do iquad=1,nquad
     do i=1,2
        do j=1,2
           do k=1,2
              do l=1,2
                 iM = iM + 1
                 Ma(iM) = F(i,k)*B(j,l) + F(j,l)*B(i,k)
                 Mi(iM) = getDOFindex1(iquad, i, j, nquad, quads)
                 Mj(iM) = getDOFindex1(iquad, k, l, nquad, quads)
              end do
           end do
        end do
     end do
  end do
  if (iM .ne. nM) then
     print *, 'Error in assembleMtx1', iM, nM
     call exit(1)
  end if

end subroutine assembleMtx1



function getDOFindex1(iquad, i, j, nquad, quads)

  implicit none
  integer, intent(in) ::  iquad, i, j, nquad, quads(nquad,4)
  integer getDOFindex1

  if ((i .eq. 1) .and. (j .eq. 1)) then
     getDOFindex1 = quads(iquad,1)
  else if ((i .eq. 2) .and. (j .eq. 1)) then
     getDOFindex1 = quads(iquad,2)
  else if ((i .eq. 2) .and. (j .eq. 2)) then
     getDOFindex1 = quads(iquad,3)
  else if ((i .eq. 1) .and. (j .eq. 2)) then
     getDOFindex1 = quads(iquad,4)
  else
     getDOFindex1 = 0
     print *, 'Error in getDOFindex1'
     call exit(1)
  end if

end function getDOFindex1




subroutine assembleRHS2(nvert, nedge, nB, &
     verts, edges, vertCon, edgeCon, B)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvert, nedge, nB, verts, edges, vertCon, edgeCon
  !f2py intent(out) B
  !f2py depend(nvert) verts, vertCon
  !f2py depend(nedge) edges, edgeCon
  !f2py depend(nB) B

  !Input
  integer, intent(in) ::  nvert, nedge, nB
  double precision, intent(in) ::  verts(nvert,2)
  integer, intent(in) ::  edges(nedge,2)
  logical, intent(in) ::  vertCon(nvert), edgeCon(nedge)

  !Output
  double precision, intent(out) ::  B(nB,2)

  !Working
  integer iB, ivert, iedge

  B(:,:) = 0.0
  iB = nB
  do iedge=nedge,1,-1
     if (edgeCon(iedge)) then
        B(iB,:) = 0.5*verts(edges(iedge,1),:) + 0.5*verts(edges(iedge,2),:)
        iB = iB - 1
     end if
  end do
  do ivert=nvert,1,-1
     if (vertCon(ivert)) then
        B(iB,:) = verts(ivert,:)
        iB = iB - 1
     end if
  end do

end subroutine assembleRHS2




subroutine assembleConMtx2(nvert, nedge, nquad, nC, &
     vertCon, edgeCon, Ca, Ci, Cj)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvert, nedge, nquad, nC, vertCon, edgeCon
  !f2py intent(out) Ca, Ci, Cj
  !f2py depend(nvert) vertCon
  !f2py depend(nedge) edgeCon
  !f2py depend(nC) Ca, Ci, Cj

  !Input
  integer, intent(in) ::  nvert, nedge, nquad, nC
  logical, intent(in) ::  vertCon(nvert), edgeCon(nedge)

  !Output
  double precision, intent(out) ::  Ca(nC)
  integer, intent(out) ::  Ci(nC), Cj(nC)

  !Working
  integer iC, ivert, iedge, nDOF, n

  Ca(:) = 1.0

  nDOF = nvert + nedge + nquad
  iC = 0
  do ivert=1,nvert
     if (vertCon(ivert)) then
        iC = iC + 1
        Ci(iC) = nDOF + iC
        Cj(iC) = ivert
     end if
  end do
  do iedge=1,nedge
     if (edgeCon(iedge)) then
        iC = iC + 1
        Ci(iC) = nDOF + iC
        Cj(iC) = nvert + iedge
     end if
  end do
  if (2*iC .ne. nC) then
     print *, 'Error in assembleConMtx2', 2*iC, nC
     call exit(1)
  end if

  n = int(nC/2)

  Cj(n+1:nC) = Ci(1:n)
  Ci(n+1:nC) = Cj(1:n)

end subroutine assembleConMtx2




subroutine assembleMtx2(nvert, nedge, nquad, nM, &
     quads, quad2edge, Ma, Mi, Mj)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) nvert, nedge, nquad, nM, quads, quad2edge
  !f2py intent(out) Ma, Mi, Mj
  !f2py depend(nquad) quads, quad2edge
  !f2py depend(nM) Ma, Mi, Mj

  !Input
  integer, intent(in) ::  nvert, nedge, nquad, nM
  integer, intent(in) ::  quads(nquad,4), quad2edge(nquad,2,2)

  !Output
  double precision, intent(out) ::  Ma(nM)
  integer, intent(out) ::  Mi(nM), Mj(nM)

  !Working
  integer getDOFindex2, iquad, iM, i, j, k, l
  double precision F(3,3), B(3,3)

  F(1,:) = (/  7. , -8. ,  1. /)
  F(2,:) = (/ -8. , 16. , -8. /)
  F(3,:) = (/  1. , -8. ,  7. /)
  F = F/3.

  B(1,:) = (/  4. ,  2. , -1. /)
  B(2,:) = (/  2. , 16. ,  2. /)
  B(3,:) = (/ -1. ,  2. ,  4. /)
  B = B/30.
  
  iM = 0
  do iquad=1,nquad
     do i=1,3
        do j=1,3
           do k=1,3
              do l=1,3
                 iM = iM + 1
                 Ma(iM) = F(i,k)*B(j,l) + F(j,l)*B(i,k)
                 Mi(iM) = getDOFindex2(iquad, i, j, nvert, nedge, nquad, quads, quad2edge)
                 Mj(iM) = getDOFindex2(iquad, k, l, nvert, nedge, nquad, quads, quad2edge)
              end do
           end do
        end do
     end do
  end do
  if (iM .ne. nM) then
     print *, 'Error in assembleMtx2', iM, nM
     call exit(1)
  end if

end subroutine assembleMtx2



function getDOFindex2(iquad, i, j, nvert, nedge, nquad, &
     quads, quad2edge)

  implicit none
  integer, intent(in) ::  iquad, i, j, nvert, nedge, nquad
  integer, intent(in) ::  quads(nquad,4), quad2edge(nquad,2,2)
  integer getDOFindex2

  if ((i .eq. 1) .and. (j .eq. 1)) then
     getDOFindex2 = quads(iquad,1)
  else if ((i .eq. 3) .and. (j .eq. 1)) then
     getDOFindex2 = quads(iquad,2)
  else if ((i .eq. 3) .and. (j .eq. 3)) then
     getDOFindex2 = quads(iquad,3)
  else if ((i .eq. 1) .and. (j .eq. 3)) then
     getDOFindex2 = quads(iquad,4)
  else if ((i .eq. 2) .and. (j .eq. 1)) then
     getDOFindex2 = nvert + abs(quad2edge(iquad,1,1))
  else if ((i .eq. 2) .and. (j .eq. 3)) then
     getDOFindex2 = nvert + abs(quad2edge(iquad,1,2))
  else if ((i .eq. 1) .and. (j .eq. 2)) then
     getDOFindex2 = nvert + abs(quad2edge(iquad,2,1))
  else if ((i .eq. 3) .and. (j .eq. 2)) then
     getDOFindex2 = nvert + abs(quad2edge(iquad,2,2))
  else if ((i .eq. 2) .and. (j .eq. 2)) then
     getDOFindex2 = nvert + nedge + iquad
  else
     getDOFindex2 = 0
     print *, 'Error in getDOFindex2'
     call exit(1)
  end if

end function getDOFindex2
