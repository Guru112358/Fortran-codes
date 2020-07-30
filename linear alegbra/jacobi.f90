program Jacobi
    !This program computes the solution to a system of equations using the Iterative Jacobi method which works when the system of equations is strongly diagonally dominant
    implicit none.
    integer,parameter::n=3
    real,dimension(n,n)::A,L,U,D,Dinv
    real,dimension(n)::b,xn,xnp1
    integer::i,j,nsteps,t
    real::tol,residual

    write(*,*)"enter the matrix A"
    read *,((A(i,j),j=1,n),i=1,n)

do i=1,n
    do j=1,n
write(*,'(f8.3,t3)',advance='no')A(i,j)
end do
write(*,*)"    "
    end do

write(*,*)"enter the matrix b"
read *,(b(j),j=1,n)

do i=1,n
    write(*,*)b(i)
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!splitting the matrices into L,U ,D and Dinv!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
L=0
U=0
D=0
Dinv=0
do i=1,n
    do j=1,n
        !!!!!!splitting diagonal!!!!!!!!!
      if(i==j)then
        D(i,j)=A(i,j)
        Dinv(i,j)=1/D(i,j)
      end if
      !!!!!!!!splitting L!!!!!!!!!!!
      if(i>j)then
        L(i,j)=-A(i,j)
      end if
    !!!!!!!!!!splitting U!!!!!!!!!!!
    if(i<j)then
        U(i,j)=-A(i,j)
    end if
    end do
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(*,*)"The Coefficient matrix A ="
do i=1,n
    do j=1,n
write(*,'(f8.3,t3)',advance='no')A(i,j)
end do
write(*,*)"    "
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

write(*,*)"The Lower triangular matrix L ="
do i=1,n
    do j=1,n
write(*,'(f8.3,t3)',advance='no')L(i,j)
end do
write(*,*)"    "
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(*,*)"The upper triangular matrix U ="
do i=1,n
    do j=1,n
write(*,'(f8.3,t3)',advance='no')U(i,j)
end do
write(*,*)"    "
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(*,*)"The Diagnonal matrix D ="
do i=1,n
    do j=1,n
write(*,'(f8.3,t3)',advance='no')D(i,j)
end do
write(*,*)"    "
end do


!!!!!!!!!!!!!!!!!Jacobi method implementation!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!initial conditions
nsteps=100
tol=0.01
xn=0


do t=1,nsteps
xnp1=MATMUL(Dinv,b+MATMUL(L+U,xn))
residual=MAXVAL(ABS(xnp1-xn))
!defining the residual
if(residual<tol)exit
xn=xnp1
if(t>nsteps)then
    write(*,*)"max number of iterations exceeded"
end if

end do

write(*,*)"The solution matrix is  X="

do j=1,n
write(*,*)xn(j)
end do


end program
