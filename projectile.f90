program Projectile
    !this program solves for the trajecttory of a projectile with quadratic air resistance using a second order modified Euler Scheme.
    !RK4 scheme gives a more accurate trajecory but this is just a simple exercise for a bored mind hell bent on world domination during the lockdown of 2020.
    implicit none
    double precision,parameter::pi=4*ATAN(1.0)
    double precision::v0,vx0,vy0,vx,vy,theta,Tp,x,y,x0,y0,dt,vxhalf,vyhalf
    integer::i,nsteps
    v0=100
    theta=45*(pi/180)
    vx0=v0*COS(theta)
    vy0=v0*SIN(theta)
    x0=0.0
    y0=0.0
    Tp=1000
    dt=0.0001
    nsteps=int(Tp/dt)
    open(1,file='xyproj.dat',status='replace')
    open(2,file='xyproj.plt',status='replace')
    do i=1,nsteps
        vxhalf=vx0+(dvx(vx0,vy0))*(dt/2)
        vyhalf=vy0+(dvy(vx0,vy0))*(dt/2)
        vx=vx0+(dvx(vxhalf,vyhalf))*(dt)
        vy=vy0+(dvy(vxhalf,vyhalf))*(dt)
        x=x0+(vx*dt/2)
        y=y0+(vy*dt/2)
        vx0=vxhalf
        vy0=vyhalf
        x0=x
        y0=y
!condition to halt the program after the projectile touches the ground(y=0)
        if(y<=0)then
        exit
        end if
        write(*,*)x,y
        write(1,*)x,y
    end do
    !writing gnuplot file to visualise x-y trajectory.
    write(2,*)'set xlabel "X"'
    write(2,*)'set ylabel "Y"'
    write(2,*)'set title "projectile trajectory"'
    write(2,*)'set size ratio 2'
    write(2,*)'plot "xyproj.dat" with line'
    CALL SYSTEM('gnuplot -p xyproj.plt')
close(1)
close(2)
contains

!functions for variation of x and y velocity components wrt time.

double precision function dvx(vx,vy)
double precision::vx,vy,c,m,g
c=0.05
g=9.81
m=150
dvx=(-c/m)*(SQRT((vx**2)+(vy**2))*vx)
end function


double precision function dvy(vx,vy)
double precision::vx,vy,c,m,g
c=0.05
g=9.81
m=150
dvy=(-g)-((c/m))*(SQRT((vx**2)+(vy**2))*vy)
end function


end program projectile
