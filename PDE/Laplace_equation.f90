
program Laplace
!-----------------------------------------------------------------------------------------------------------!
!!!!!!!!!!!!!This program solves the Laplace equation(uxx+uyy=0) by simultaneous relaxation !!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------------------------------------------!
implicit none
integer,parameter::nx=2,ny=2
integer::i,j,nsteps,k
real::a,b,dx,dy
real,dimension(0:nx+1,0:ny+1)::phi,phin,x,y
real,parameter::pi=4*ATAN(1.0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
a=1.0 ! length of domain
b=1.0  !breadth of domain
dx=b/(nx)
dy=a/(ny)
!defining boundary conditions
phi=0.0 !relaxing inner points as zero
!!!!!!!!!!!!!!!!!!!!!!!! Setting Boundary conditions,selecting a simple grid to verify results quickly!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i=0,nx+1
do j=0,ny+1
x(i,j)=i*dx
y(i,j)=j*dy
!top boundary
phi(0,j)=1
!bottom boundary
phi(nx+1,j)=2
!left  boundary
phi(i,0)=2
!right boundary
phi(i,ny+1)=1
end do
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(*,*)"Initial Boundary values are: "
write(*,*)
do i=0,nx+1
 do j=0,ny+1
write(*,'(f8.3,t5)',advance='no')phi(i,j)
end do
write(*,*)
end do
!-----------------------------------------------------------------------------------------------------------!
!!!!!!!!!!!!!solving the equations for the above boundary conditions by Jacobi iterations!!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------------------------------------------!
nsteps=100

write(*,*)
write(*,*)"the final solution is:"
write(*,*)

do k=1,nsteps
do i=1,nx
do j=1,ny
phin(i,j)=0.25*(phi(i+1,j)+phi(i,j+1)+phi(i,j-1)+phi(i-1,j))
phi(i,j)=phin(i,j)
end do
end do
end do

!printing the final solution
do i=1,nx
do j=1,ny
write(*,'(f8.3,t5)',advance='no')phin(i,j)
end do
write(*,*)
end do



end program
