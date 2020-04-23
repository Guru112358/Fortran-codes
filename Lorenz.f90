program Lorenz
    ! This is a simple program to solve the Lorenz series of 3 ODE's by taking the constants rho,sigma and beta such that chaotic behaviour is observed from the literature.
    !the numerical method employed is a simple midpoint method.
    !visualisation using gnuplot.
    implicit none
    real::x,y,z,xbar,ybar,zbar,x__,y__,z__,dt,Tp,sigma,rho,beta
    integer::i,nsteps
    Tp=100
    dt=0.01
    nsteps=INT(Tp/dt)
    x=1
    y=2
    z=5
    sigma=10
    rho=28
    beta=8/3
    open(1,file='lorenz.dat',status='replace')
    open(2,file='phaseplotlorenz.plt',status='replace')
    do i=1,nsteps
        xbar=x+(0.5*dt)*(sigma*(y-x))
        ybar=y+(0.5*dt)*(x*(rho-z)-y)
        zbar=z+(0.5*dt)*((x*y)-(beta*z))

        x__=x+(dt)*(sigma*(ybar-xbar))
        y__=y+(dt)*(xbar*(rho-zbar)-ybar)
        z__=z+(dt)*((xbar*ybar)-(beta*zbar))

        x=x__
        y=y__
        z=z__

        write(1,*)x,y,z
    end do

        write(2,*)'set xlabel "x"'
        write(2,*)'set ylabel "y"'
        write(2,*)'set zlabel "z"'
        write(2,*)'set title "3D phase plot"'
        write(2,*)'splot "lorenz.dat" with line'
        CALL SYSTEM('gnuplot -p phaseplotlorenz.plt')


close(1)
close(2)

end program
