program Projectile
    !this program solves for the trajecttory of a projectile with quadratic air resistance  Euler Scheme.
    implicit none
    double precision,parameter::pi=4*ATAN(1.0)
    double precision::v0,vx0,vy0,vx,vy,theta,x,y,x0,y0,dt,vxbar,vybar,vxnp1,vynp1,xnp1,ynp1,xbar,ybar
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
    switch1=0

    open(1,file='xyproj.dat',status='replace')
    open(2,file='xyproj.plt',status='replace')

write(*,*)"please wait,computing trajectory......"
write(*,*)

do while(switch1.EQ.0)
        t=t+1  
        vxnp1=vx0+dvx(vx0,vy0,m,c,g)*(dt)
        vynp1=vy0+dvy(vx0,vy0,m,c,g)*(dt)
        xnp1=x0+(vx0*dt)+(0.5*dt**2)*(dvx(vx0,vy0,m,c,g))
        ynp1=y0+(vy0*dt)+(0.5*dt**2)*(dvy(vx0,vy0,m,c,g))
        x0=xnp1
        y0=ynp1
        vx0=vxnp1
        vy0=vynp1

        if(ynp1<=0)then  !condition to halt the program after the projectile touches the ground(y=0)
        switch1=1
         write(*,*)"Integration complete!,plotting results"
         write(*,*)"Time of flight:",t*dt,"seconds"
        else
        write(*,*)xnp1,ynp1
        write(1,*)xnp1,ynp1
      
       end if   
 end do

    !writing gnuplot file to visualise x-y trajectory.
    write(2,*)'set xlabel "Range(m)"'
    write(2,*)'set ylabel "Height(m)"'
    write(2,*)'set grid'
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




