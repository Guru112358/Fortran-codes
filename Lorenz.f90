program Lorenz
implicit none
! This is a simple program to solve the Lorenz series of 3 ODE's by taking the constants rho,sigma and beta such that chaotic behaviour is observed from the literature.
!simple 4 th order Runge kutta method is utiliased to solve the system of ODEs
!visualisation using gnuplot.

real(8),parameter::s=10.0 !sigma
real(8),parameter::b=2.666666667  !beta
real(8),parameter::r=28  !rho
real(8)::x0,y0,z0,xnp1,ynp1,znp1,Tp,dt,xbar,ybar,zbar
real(8)::k1,k2,k3,k4,l1,l2,l3,l4,m1,m2,m3,m4
real(8),allocatable,dimension(:)::x,y,z,time
integer::i,nsteps


Tp=50
dt=0.0001
x0=0.0
y0=1.0
z0=0.0
nsteps=int(Tp/dt)


open(1,file='xyz.dat',status='replace')
open(2,file='plot1.plt',status='replace')
open(3,file='plot2.plt',status='replace')
open(4,file='plot3.plt',status='replace')
open(5,file='plot4.plt',status='replace')

allocate(x(nsteps),y(nsteps),z(nsteps),time(nsteps))


do i=1,nsteps


k1=dxdt(x0,y0)
l1=dydt(x0,y0,z0)
m1=dzdt(x0,y0,z0)

k2=dxdt(x0+0.5*k1*dt,y0+0.5*l1*dt)
l2=dydt(x0+0.5*k1*dt,y0+0.5*l1*dt,z0+0.5*m1*dt)
m2=dzdt(x0+0.5*k1*dt,y0+0.5*l1*dt,z0+0.5*m1*dt)

k3=dxdt(x0+0.5*k2*dt,y0+0.5*l2*dt)
l3=dydt(x0+0.5*k2*dt,y0+0.5*l2*dt,z0+0.5*m2*dt)
m3=dzdt(x0+0.5*k2*dt,y0+0.5*l2*dt,z0+0.5*m2*dt)

k4=dxdt(x0+k3*dt,y0+l3*dt)
l4=dydt(x0+k3*dt,y0+l3*dt,z0+m3*dt)
m4=dzdt(x0+k3*dt,y0+k3*dt,z0+m3*dt)

xnp1=x0+(dt/6)*(k1+(2*k2)+(2*k3)+k4)
ynp1=y0+(dt/6)*(l1+(2*l2)+(2*l3)+l4)
znp1=z0+(dt/6)*(m1+(2*m2)+(2*m3)+m4)

x0=xnp1
y0=ynp1
z0=znp1


x(i)=xnp1
y(i)=ynp1
z(i)=znp1

time(i)=i*dt

write(1,*)x(i),y(i),z(i),time(i)
!write(*,*)x(i),y(i),z(i),time(i)
!write(*,*)Jacobian(xnp1,ynp1,znp1)


end do

write(2,*)"set xlabel 'x'"
write(2,*)"set ylabel 'y'"
write(2,*)"set zlabel 'z'"
write(2,*)"splot 'xyz.dat'  with line lt rgb 'red' title 'Phase Plot' "
write(3,*)"set xlabel 'time'"
write(3,*)"set grid"
write(3,*)"set mouse"
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
real(8)::x,y,xdot
xdot=s*(y-x)
return
end function

function dydt(x,y,z)result(ydot)
real(8)::x,y,z,ydot
ydot=x*(r-z)-y
return
end function


function dzdt(x,y,z)result(zdot)
real(8)::x,y,z,zdot
zdot=(x*y)-(b*z)
return
end function


function Jacobian(x,y,z)result(J)
real(8)::x,y,z
real(8),dimension(3,3)::J
J= reshape((/-s,r-z,y,s,-1.00d0,x,0.00d0,-x,-b /), (/3,3/))
return
end function

end program

