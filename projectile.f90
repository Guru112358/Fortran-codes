
program Projectile
    !this program solves for the trajecttory of a projectile with quadratic air resistance using a second order modified Euler Scheme.
    implicit none
    double precision,parameter::pi=4*ATAN(1.0)
    double precision::v0,vx0,vy0,vx,vy,theta,x,y,x0,y0,dt,vxhalf,vyhalf
    double precision::m  !mass of the projectile
    double precision::c   !drag coefficient
    double precision,parameter::g=9.81 !acceleration due to gravity
    integer::i,switch1,t

    write(*,*)"Enter the mass of the projectile(kg):- "
    read *,m
    write(*,*)
    write(*,*)"Enter the Drag Coefficent of the projectile:- "
    read *,c
    write(*,*)
    write(*,*)"Enter the Coordinates of the projectile launcher(x,y):- "
    read *,x0,y0
    write(*,*)
    write(*,*)"Enter the angle of elevation in degrees:-"
    read *,theta
    write(*,*)
    write(*,*)"Enter The Projection velocity(m/s):-"
    read *,v0
    write(*,*)
    write(*,*)"Enter the time step for integration:- "
    read *,dt
    write(*,*)

    vx0=v0*COS(theta*(pi/180))  !x velocity
    vy0=v0*SIN(theta*(pi/180))  ! y velocity
    t=0
    switch1=1

    open(1,file='xyproj.dat',status='replace')
    open(2,file='xyproj.plt',status='replace')

write(*,*)"please wait,computing trajectory......"
write(*,*)

do while(switch1.GT.0)
        t=t+1  
        vxhalf=vx0+(dvx(vx0,vy0,m,c,g))*(dt/2)
        vyhalf=vy0+(dvy(vx0,vy0,m,c,g))*(dt/2)
        vx=vx0+(dvx(vxhalf,vyhalf,m,c,g))*(dt)
        vy=vy0+(dvy(vxhalf,vyhalf,m,c,g))*(dt)
        x=x0+(vx*dt/2)
        y=y0+(vy*dt/2)
        vx0=vxhalf
        vy0=vyhalf
        x0=x
        y0=y
        if(y<=0)then  !condition to halt the program after the projectile touches the ground(y=0)
        switch1=-1
         write(*,*)"Integration complete!,plotting results"
         write(*,*)"Time of flight:",t*dt,"seconds"
        else
        !write(*,*)x,y
        write(1,*)x,y
      
       end if   
 end do

    !writing gnuplot file to visualise x-y trajectory.
    write(2,*)'set xlabel "Range(m)"'
    write(2,*)'set ylabel "Height(m)"'
    write(2,*)'set grid'
    write(2,*)'set autoscale xy'
    write(2,*)'plot "xyproj.dat" with line lt rgb "red" title "Trajectory"'
    CALL SYSTEM('gnuplot -p xyproj.plt')
close(1)
close(2)

contains

!functions for variation of x and y velocity components wrt time.

double precision function dvx(vx,vy,m,c,g)
implicit none
double precision::vx,vy,m,c,g

dvx=(-c/m)*(SQRT((vx**2)+(vy**2))*vx)

return
end function


double precision function dvy(vx,vy,m,c,g)
implicit none
double precision::vx,vy,m,c,g

dvy=(-g)-((c/m))*(SQRT((vx**2)+(vy**2))*vy)

return
end function


end program projectile




