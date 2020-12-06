    program lotka_volterra
            implicit none
            !This program solves the lotka Volterra equations of a apex predator,middle predator and one prey species and single prey species using Heun's method.
            !DISCLAIMERS: 1)The model used here is highly simplistic,the populations and time are non dimensional.
            !             2)The values of the coefficients are just an example and not represntative of any real world system.
            !             3)Running ths program for longer time duration causes numerical error to build up as the method is only second order accurate.
            !             4)The Runge-Kutta-Fehlberg method provides better results for more complicated population models.
            !             5)e represents the effect of predation on species y by species z,
            !               f represents the natural death rate of the predator z in the absence of prey,
            !               g represents the efficiency and propagation rate of the predator z in the presence of prey

            !             6)my reference: i)https://web.ma.utexas.edu/users/davis/375/popecol/lec10/lotka.html
            !                            ii)https://en.wikipedia.org/wiki/LotkaVolterra_equations#cite_note-27
            !                            iii)A Lotka-Volterra Three-Species Food Chain Author(s): Erica Chauvet, Joseph E. Paullet, Joseph P. Previte and Zac WallsSource: Mathematics Magazine, Vol. 75, No. 4 (Oct., 2002), pp. 243-255Published
            !
            real::x,y,alpha,beta,gammar,delta,xbar,ybar,x__,y__,dt,z,z__,e,f,g,zbar
            integer::i,nsteps,time
            alpha=0.8
            beta=0.65
            gammar=0.45
            delta=0.60
            e=0.7
            f=0.5
            g=0.6
            !intial conditions with equal populations of apex predator(Z), middle predator(y) and prey(x)
            x=1.5
            y=0.8
            z=1.0
            dt=0.001
            time=30
            nsteps=int(time/dt)
            open(1,file="xpopulation3.dat",status='replace')
            open(2,file="ypopulation3.dat",status='replace')
            open(6,file="xyphaseplot3.plt",status='replace')
            open(7,file="zpopulation3.dat",status='replace')
            open(8,file='xyzplot3.plt',status='replace')
            open(9,file='xyzphaseplot3.plt',status='replace')
            open(10,file='xyz3.dat',status='replace')

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            do i=1,nsteps
        !using midpoint method/Heun's method to calculate the populations at half time step
           xbar=x+(0.5*dt)*((alpha*x)-(beta*x*y))
           ybar=y+(0.5*dt)*((delta*x*y)-(gammar*y)-(e*y*z))
           zbar=z+(0.5*dt)*((-f*z)+(g*y*z))
        !computing the final values of x,y and z by time stepping the midpoint values to the original values
           x__=x+(dt)*((alpha*xbar)-(beta*xbar*ybar))
           y__=y+(dt)*((delta*xbar*ybar)-(gammar*ybar)-(e*ybar*zbar))
           z__=z+(dt)*((-f*zbar)+(g*ybar*zbar))
           x=x__
           y=y__
           z=z__
           write(1,*)(i*dt),x
           write(2,*)(i*dt),y
           write(7,*)(i*dt),z
           write(10,*)x,y,z
            end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        write(8,*)'set xlabel "time"'
        write(8,*)'set ylabel "population"'
        write(8,*)'set title "Species population vs time"'
        write(8,*)'plot "xpopulation3.dat" with line,"zpopulation3.dat" with line,\'
        write(8,*)'"ypopulation3.dat" with line'
        CALL SYSTEM('gnuplot -p xyzplot3.plt')
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        write(9,*)'set xlabel "prey population"'
        write(9,*)'set ylabel "middle predator population"'
        write(9,*)'set zlabel "apex predator population"'
        write(9,*)'set title "3D phase plot"'
        write(9,*)'splot "xyz3.dat" with line'
        CALL SYSTEM('gnuplot -p xyzphaseplot3.plt')
        close(1)
        close(2)
        close(6)
        close(7)
        close(8)
        close(9)
        close(10)
    end program

