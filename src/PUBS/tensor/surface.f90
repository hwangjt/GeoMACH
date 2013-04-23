subroutine evalsurface(k1, k2, m1, m2, n1, n2, B1, B2, i0, j0, C, P)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) k1,k2,m1,m2,n1,n2,B1,B2,i1,i2,C
  !f2py intent(out) P
  !f2py depend(n1,n2,k1) B1
  !f2py depend(n1,n2,k2) B2
  !f2py depend(n1,n2) i0
  !f2py depend(n1,n2) j0
  !f2py depend(m1,m2) C
  !f2py depend(n1,n2) P

  !Input
  integer, intent(in) ::  k1, k2, m1, m2, n1, n2
  double precision, intent(in) ::  B1(n1,n2,k1), B2(n1,n2,k2)
  integer, intent(in) ::  i0(n1,n2), j0(n1,n2)
  double precision, intent(in) ::  C(m1,m2)

  !Output
  double precision, intent(out) ::  P(n1,n2)

  !Working
  integer i, j, u, v

  P(:,:) = 0
  do j=1,k2
     do i=1,k1
        do v=1,n2
           do u=1,n1
              P(u,v) = P(u,v) + B1(u,v,i)*B2(u,v,j)*C(i0(u,v)+i,j0(u,v)+j)
           end do
        end do
     end do
  end do

end subroutine evalsurface


subroutine surfacejacob(k1, k2, m1, m2, n1, n2, B1, B2, i0, j0, dPdC)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) k1,k2,m1,m2,n1,n2,B1,B2,i0,j0
  !f2py intent(out) dPdC
  !f2py depend(n1,n2,k1) B1
  !f2py depend(n1,n2,k2) B2
  !f2py depend(n1,n2) i0
  !f2py depend(n1,n2) j0
  !f2py depend(n1,n2,m1,m2) dPdC

  !Input
  integer, intent(in) ::  k1, k2, m1, m2, n1, n2
  double precision, intent(in) ::  B1(n1,n2,k1), B2(n1,n2,k2)
  integer, intent(in) ::  i0(n1,n2), j0(n1,n2)

  !Output
  double precision, intent(out) ::  dPdC(n1*n2,m1*m2)

  !Working
  integer i, j, u, v, uv

  dPdC(:,:) = 0
  
  do j=1,k2
     do i=1,k1
        uv = 1
        do v=1,n2
           do u=1,n1
              dPdC(uv,i0(u,v)+i+(j0(u,v)+j-1)*m1) = B1(u,v,i)*B2(u,v,j)
              uv = uv + 1
           end do
        end do
     end do
  end do
  
end subroutine surfacejacob


subroutine surfacefit(k1, k2, m1, m2, n1, n2, B1, B2, i0, j0, P0, C, dPdC, P)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) k1,k2,m1,m2,n1,n2,B1,B2,i0,j0,P0,C
  !f2py intent(out) dPdC,P
  !f2py depend(n1,n2,k1) B1
  !f2py depend(n1,n2,k2) B2
  !f2py depend(n1,n2) i0
  !f2py depend(n1,n2) j0
  !f2py depend(n1,n2) P0
  !f2py depend(m1,m2) C
  !f2py depend(n1,n2,m1,m2) dPdC
  !f2py depend(n1,n2) P

  !Input
  integer, intent(in) ::  k1, k2, m1, m2, n1, n2
  double precision, intent(in) ::  B1(n1,n2,k1), B2(n1,n2,k2)
  integer, intent(in) ::  i0(n1,n2), j0(n1,n2)
  double precision, intent(in) ::  P0(n1,n2,3)
  double precision, intent(in) ::  C(m1,m2,3)

  !Output
  double precision, intent(out) ::  dPdC((n1-2)*(n2-2),(m1-2)*(m2-2))
  double precision, intent(out) ::  P((n1-2)*(n2-2),3)

  !Working
  integer i, j, u, v, l
  double precision dPdC0(n1*n2,m1*m2)
  double precision sum

  dPdC0(:,:) = 0  
  do j=1,k2
     do i=1,k1
        do v=1,n2
           do u=1,n1
              dPdC0(u+(v-1)*n1,i0(u,v)+i+(j0(u,v)+j-1)*m1) = B1(u,v,i)*B2(u,v,j)
           end do
        end do
     end do
  end do
  do j=1,m2-2
     do i=1,m1-2
        do v=1,n2-2
           do u=1,n1-2
              dPdC(u+(v-1)*(n1-2),i+(j-1)*(m1-2)) = dPdC0(u+1+v*n1,i+1+j*m1)
           end do
        end do
     end do
  end do
  do l=1,3
     do v=1,n2-2
        do u=1,n1-2
           sum = 0.0
           do i=1,m1
              sum = sum &
                   + dPdC0(u+1+v*n1,i)*C(i,1,l) &
                   + dPdC0(u+1+v*n1,i+(m2-1)*m1)*C(i,m2,l)
           end do
           do j=2,m2-1
              sum = sum &
                   + dPdC0(u+1+v*n1,1+(j-1)*m1)*C(1,j,l) &
                   + dPdC0(u+1+v*n1,m1+(j-1)*m1)*C(m1,j,l)
           end do
           P(u+(v-1)*(n1-2),l) = P0(u+1,v+1,l) - sum
        end do
     end do
  end do
           
end subroutine surfacefit


subroutine expand2d(m1, m2, X, C0, C)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) m1,m2,X,C0
  !f2py intent(out) C
  !f2py depend(m1,m2) X
  !f2py depend(m1,m2) C0
  !f2py depend(m1,m2) C

  !Input
  integer, intent(in) ::  m1, m2
  double precision, intent(in) ::  X((m1-2)*(m2-2))
  double precision, intent(in) ::  C0(m1,m2)

  !Output
  double precision, intent(out) ::  C(m1,m2)

  !Working
  integer i,j,ij

  C(:,:) = C0(:,:)

  ij = 1
  do j=1,m2-2
     do i=1,m1-2
        C(i+1,j+1) = X(i+(j-1)*(m1-2)) !ij)
        ij = ij + 1
     end do
  end do
  
end subroutine expand2d


subroutine flatten2d(n1, n2, P, Q)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) n1,n2,P
  !f2py intent(out) Q
  !f2py depend(n1,n2) P
  !f2py depend(n1,n2) Q

  !Input
  integer, intent(in) ::  n1, n2
  double precision, intent(in) ::  P(n1,n2)

  !Output
  double precision, intent(out) ::  Q(n1*n2)

  !Working
  integer u,v,uv

  Q(:) = 0

  uv = 1
  do v=1,n2
     do u=1,n1
        Q(uv) = P(u,v)
        uv = uv + 1
     end do
  end do
  
end subroutine flatten2d


subroutine expand2d0(m1, m2, D, C)

  implicit none

  !Fortran-python interface directives
  !f2py intent(in) m1,m2,D
  !f2py intent(out) C
  !f2py depend(m1,m2) D
  !f2py depend(m1,m2) C

  !Input
  integer, intent(in) ::  m1, m2
  double precision, intent(in) ::  D(m1*m2)

  !Output
  double precision, intent(out) ::  C(m1,m2)

  !Working
  integer i,j,ij

  C(:,:) = 0

  ij = 1
  do j=1,m2
     do i=1,m1
        C(i,j) = D(ij)
        ij = ij + 1
     end do
  end do
  
end subroutine expand2d0
