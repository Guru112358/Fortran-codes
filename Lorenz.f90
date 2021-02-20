program Lorenz
implicit none

! This is a simple program to solve the Lorenz series of 3 ODE's by taking the constants rho,sigma and beta such that chaotic behaviour is observed from the literature.
!the numerical method employed is a simple midpoint method.
!visualisation using gnuplot.

real,parameter::s=10.0 !sigma
real,parameter::b=8/3  !beta
real,parameter::r=28   !rho
real::x0,y0,z0,xnp1,ynp1,znp1,Tp,dt,xbar,ybar,zbar
real,allocatable,dimension(:)::x,y,z,time
integer::i,nsteps


Tp=30
dt=0.0001
x0=0.1
y0=0.0
z0=0.0
nsteps=int(Tp/dt)


open(1,file='xyz.dat',status='replace')
open(2,file='plot1.plt',status='replace')
open(3,file='plot2.plt',status='replace')
open(4,file='plot3.plt',status='replace')
open(5,file='plot4.plt',status='replace')

allocate(x(nsteps),y(nsteps),z(nsteps),time(nsteps))


do i=1,nsteps

xbar=x0+0.5*dt*dxdt(x0,y0)
ybar=y0+0.5*dt*dydt(x0,y0,z0)
zbar=z0+0.5*dt*dzdt(x0,y0,z0)

xnp1=x0+dt*dxdt(xbar,ybar)
ynp1=y0+dt*dydt(xbar,ybar,zbar)
znp1=z0+dt*dzdt(xbar,ybar,zbar)


x0=xnp1
y0=ynp1
z0=znp1

x(i)=xnp1
y(i)=ynp1
z(i)=znp1
time(i)=i*dt
write(1,*)x(i),y(i),z(i),time(i)
end do

write(2,*)"set xlabel 'x'"
write(2,*)"set ylabel 'y'"
write(2,*)"set zlabel 'z'"
write(2,*)"splot 'xyz.dat'  with line "
write(3,*)"set xlabel 'time'"
write(3,*)"set grid"
write(3,*)"plot 'xyz.dat' using 4:1 with line lt rgb 'red' title 'X'"
write(3,*)"set xlabel 'time'"
write(4,*)"set grid"
write(4,*)"set xlabel 'time'"
write(4,*)"plot 'xyz.dat' using 4:2 with line lt rgb 'blue' title 'y'"
write(5,*)"set grid"
write(5,*)"set xlabel 'time'"
write(5,*)"plot 'xyz.dat' using 4:3 with line lt rgb 'dark-violet' title 'Z'"

call system('gnuplot -p plot1.plt')
call system('gnuplot -p plot2.plt')
call system('gnuplot -p plot3.plt')
call system('gnuplot -p plot4.plt')

deallocate(x,y,z,time)

 close(1)
 close(2)
 close(3)
 close(4)
 close(5)

contains

function dxdt(x,y)result(xdot)
real::x,y,xdot
xdot=s*(y-x)
return
end function

function dydt(x,y,z)result(ydot)
real::x,y,z,ydot
ydot=x*(r-z)-y
return
end function


function dzdt(x,y,z)result(zdot)
real::x,y,z,zdot
zdot=(x*y)-(b*z)
return
end function


function Jacobian(x,y,z)result(J)
real::x,y,z
real,dimension(3,3)::J
J= reshape((/-s,r-z,y,s,-1.00,x,0.00,-x,-b /), (/3,3/))
return
end function

end program
