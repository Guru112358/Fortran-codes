program heat2D
    !this program solves the 2D heat equation with a given initial heat distribution with a simple explicit Finite difference method.
    !the scheme used is a FTCS scheme.
    !warning::(alpha*dt)/(dx**2)+*(beta*dt)/(dy**2) should be lesser than or equal to half(0.5) for a stable solution.
    !the domain is a square plate of dimensions 1 x 1
    implicit none
    integer::i,Tp,nsteps,time,j,nx,ny
    real,dimension(50,80)::T,Te,x,y
    real::dx,D,alpha,dy,beta,dt
    D=0.0000001
    dt=0.1
    nx=50 !grid points in x direction
    ny=80 !grid points in y direction
    dx=1.0/nx
    dy=1.0/ny
    Tp=500
    nsteps=INT(Tp/dt)
    alpha=(D*dt)/(dx**2)
    beta=(D*dt)/(dy**2)
    open(1, file = 'initial.dat', status = 'replace')
    open(2,file='tempdata.plt',status='replace')
    open(3,file = 'final.dat', status = 'replace')
    open(4,file = 'initial1.plt', status = 'replace')
  !Defining initial heat distribution with a function.

    do i=0,nx
        do j=0,ny
            x(i,j)=i*dx
            y(i,j)=j*dy
            T(i,j)=500*SIN(x(i,j)+y(i,j))
            write(1,*)x(i,j),y(i,j),T(i,j)
        end do
    end do

    do time=1,nsteps !time loop
            do i=1,nx-1 !x traversal loop
            do j=1,ny-1 !y traversal loop
                Te(i,j)=alpha*(T(i+1,j)-T(i-1,j))+beta*(T(i,j+1)-T(i,j-1))+(1-2*alpha-2*beta)*T(i,j)!2D FTCS time stepping
                if(time==nsteps)then
                write(3,*)x(i,j),y(i,j),Te(i,j)
                else if(Te(i,j)<0.0)Then
                    exit
                end if
	        end do
            end do
	        T=Te
        end do
!Post processing the data files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            write(4,*)'set xlabel "x"'
            write(4,*)'set ylabel "y"'
            write(4,*)'set zlabel "T"'
            write(4,*)'set grid'
            write(4,*)'set title "Initial temperature distribution"'
            write(4,*)'plot "initial.dat" u 1:2:3 with image'
            CALL SYSTEM('gnuplot -p initial1.plt')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            write(2,*)'set xlabel "x"'
            write(2,*)'set ylabel "y"'
            write(2,*)'set zlabel "T"'
            write(2,*)'set grid'
            write(2,*)'set title "final temperature distribution"'
            write(2,*)'plot "final.dat" u 1:2:3 with image'
            CALL SYSTEM('gnuplot -p tempdata.plt')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    close(1)
    close(2)
    close(3)
    close(4)

end program heat2D
