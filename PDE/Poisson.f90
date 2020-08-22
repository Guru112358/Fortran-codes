program Laplace
!-----------------------------------------------------------------------------------------------------------
!!!!!!!!!!!!!This program solves the Poisson's equation(uxx+uyy=rho(x,y)) by simultaneous relaxation in 2D
!These problems are typically encountered in electrostatics in order to find the potential of a given charge distribution (rho) in space.
!-----------------------------------------------------------------------------------------------------------
implicit none
integer,parameter::nx=50,ny=50
integer::i,j,nsteps,k
real::a,b,dx,dy,h
real,dimension(0:nx+1,0:ny+1)::phi,phin,x,y,rho
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
a=1.0 ! length of domain
b=1.0  !breadth of domain
dx=b/(nx)
dy=a/(ny)
!defining boundary conditions
phi=0.0 !relaxing inner points as zero
!!!!!!!!!!!!!!!!!!!!!!!!Defining the charge distrbiution!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
open(1,file='xyrho.dat',status='replace')
open(2,file='xyrho.plt',status='replace')
do i=0,nx+1
do j=0,ny+1
x(i,j)=i*dx
y(i,j)=j*dy
rho(i,j)=x(i,j)**2+y(i,j)**2
write(1,*)x(i,j),y(i,j),rho(i,j)
end do
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
h=1/(1+nx) !spacing parameter
!!!!!!!!!plotting charge distribution!!!!!!!!!!
write(2,*)"set xlabel 'x'"
write(2,*)"set ylabel 'y'"
write(2,*)"set zlabel 'rho'"
write(2,*)"splot 'xyrho.dat'"
CALL SYSTEM('gnuplot -p xyrho.plt')
close(1)
close(2)

!!!!!!!!!!!!!!!!!!!!!!!! Setting Boundary conditions!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i=0,nx+1
do j=0,ny+1
x(i,j)=i*dx
y(i,j)=j*dy
!top boundary
phi(0,j)=1
!bottom boundary
phi(nx+1,j)=2
!left  boundary
phi(i,0)=3
!right boundary
phi(i,ny+1)=4
end do
end do
!-----------------------------------------------------------------------------------------------------------!
!!!!!!!!!!!!!solving the Poisson equation for the above boundary conditions by Jacobi iterations!!!!!!!!!!!!!
!-----------------------------------------------------------------------------------------------------------!
nsteps=200

write(*,*)
write(*,*)"the final solution is:"
write(*,*)
open(3,file="phixy.dat",status="replace")
open(4,file="phixy1.plt",status="replace")
do k=1,nsteps
do i=1,nx
do j=1,ny
x(i,j)=i*dx
y(i,j)=j*dy
phin(i,j)=0.25*(phi(i+1,j)+phi(i,j+1)+phi(i,j-1)+phi(i-1,j)-((4*h**2)*rho(i,j)))
phi(i,j)=phin(i,j)
write(3,*)x(i,j),y(i,j),phi(i,j)
end do
end do
end do
!plotting the obtained Potential for the given charge dsitribution
write(4,*)"set xlabel 'x'"
write(4,*)"set ylabel 'y'"
write(4,*)"set zlabel 'phi'"
write(4,*)"splot 'phixy.dat'"
CALL SYSTEM('gnuplot -p phixy1.plt')
close(3)
close(4)
end program
