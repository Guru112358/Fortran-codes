    program lotka_volterra
        implicit none
        !This program solves the lotka Volterra equations of a single predator species and single prey species using Heun's method.
        !DISCLAIMERS: 1)The model used here is highly simplistic,the populations and time are non dimensional.
        !             2)The values of the coefficients are just an example and not represntative of any real world system.
        !             3)Running ths program for longer time duration causes numerical error to build up as the method is only second order accurate.
        !             4)The Runge-Kutta-Fehlberg method provides better results for more complicated population models.
        !             5)my reference: i)https://web.ma.utexas.edu/users/davis/375/popecol/lec10/lotka.html
        !                            ii)https://en.wikipedia.org/wiki/Lotka–Volterra_equations#cite_note-27
        real::x,y,alpha,beta,gammar,delta,xbar,ybar,x__,y__,dt
        integer::i,nsteps,time
        alpha=0.5
        beta=0.8
        gammar=0.25
        delta=0.1
        !intial conditions with equal populations of predator(y) and prey(x)
        x=0.8
        y=0.1
        dt=0.01
        time=100
        nsteps=int(time/dt)
        open(1,file="xpopulation.dat",status='new')
        open(2,file="ypopulation.dat",status='new')
        open(3,file="xplot.plt",status='new')
        open(4,file="yplot.plt",status='new')
        open(5,file="xvsy.dat",status='new')
        open(6,file="xyphaseplot.plt",status='new')
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do i=-1,nsteps
    !using midpoint method/Heun's method to calculate the populations at half time step
       xbar=x+(0.5*dt)*((alpha*x)-(beta*x*y))
       ybar=y+(0.5*dt)*((delta*x*y)-(gammar*y))
    !computing the final values of x and y by time stepping the midpoint values to the original values
       x__=x+(dt)*((alpha*xbar)-(beta*xbar*ybar))
       y__=y+(dt)*((delta*xbar*ybar)-(gammar*ybar))
       x=x__
       y=y__
       write(1,*)(i*dt),x
       write(2,*)(i*dt),y
       write(5,*)x,y
        end do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(3,*)'set xlabel "time"'
    write(3,*)'set ylabel "population of prey"'
    write(3,*)'set title "prey population vs time"'
    write(3,*)'plot "xpopulation.dat" with line'
    CALL SYSTEM('gnuplot -p xplot.plt')
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(4,*)'set xlabel "time"'
    write(4,*)'set ylabel "population of predator"'
    write(4,*)'set title "predator population vs time" '
    write(4,*)'plot "ypopulation.dat" with line'
    CALL SYSTEM('gnuplot -p yplot.plt')
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(6,*)'set xlabel "prey population"'
    write(6,*)'set ylabel "predator population"'
    write(6,*)'set title "phase plot"'
    write(6,*)'plot "xvsy.dat" with line'
    CALL SYSTEM('gnuplot -p xyphaseplot.plt')
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    close(1)
    close(2)
    close(3)
    close(4)
    close(5)
    close(6)

    end program
