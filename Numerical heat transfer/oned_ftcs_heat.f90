program oneDheatequation
    !this program solves the 1 D heat equation with a givwen initial heat distribution with a simple explicit Finite difference method
    !the scheme used is a FTCS scheme
    !warning::(alpha*dt)/(dx**2) should be lesser than or equal to one for a stable solution
integer::i,n,Tp,nsteps,time
real,allocatable,dimension(:)::T,x,Te
real::dx,alpha,dt,C
real,parameter::pi=4*ATAN(1.0)

alpha=0.001 !thermal diffusivity of the  material
n=100       !number of grid points
dt=0.01     !time step
dx=1.0/n    !grid spacing
Tp=200      !time period
allocate(T(n),Te(n),x(n))
nsteps=INT(Tp/dt)
C=(alpha*dt)/(dx**2)

open(1,file="temperature.dat",status='replace')
open(2,file="tempdata.plt",status='replace')
open(3,file="temptime.dat",status='replace')
open(4,file="temptime.plt",status='replace')
open(5,file='temptime3d.dat',status='replace')
open(6,file='temptime3d.plt',status='replace')

!defining heat dsitribution at time t-0,for simplicity,a sinusoidal distribution is used.However other distributions can be specified.

do i=0,n
x(i)=i*dx
T(i)=500*SIN(pi*x(i))
end do


do time=1,nsteps !time loop

do i=1,n-1       !space loop

Te(i)=T(i)+C*((T(i+1)-2*T(i)+T(i-1)))  !discretised time stepping
if(MOD(time,100)==0)then
write(1,*)x(i),T(i)
write(3,*)time*dt,T(i)
write(5,*)time*dt,x(i),T(i)
end if
end do
T=Te

end do

!post processing data files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(2,*)'set xlabel "X"'
write(2,*)'set ylabel "Temperature"'
write(2,*)'plot "temperature.dat" '
CALL SYSTEM('gnuplot -p tempdata.plt')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(4,*)'set xlabel "time"'
write(4,*)'set ylabel "Temperature"'
write(4,*)'set title "3D heat evolution"'
write(4,*)'plot "temptime.dat" '
CALL SYSTEM('gnuplot -p temptime.plt')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(6,*)'set xlabel "time"'
write(6,*)'set ylabel "x"'
write(6,*)'set zlabel "Temperature"'
write(6,*)'set title "temperature time plot"'
write(6,*)'splot "temptime3d.dat" '
CALL SYSTEM('gnuplot -p temptime3d.plt')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

close(1)
close(2)
close(3)
close(4)
close(5)
close(6)
deallocate(x,T,Te)
end program oneDheatequation
