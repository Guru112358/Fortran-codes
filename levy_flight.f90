program Levy_flight

implicit none

real,parameter::pi=4*ATAN(1.0)

integer::i,nsteps

real:: x0,y0,z0,p,theta,step,q,phi


open(1,file='levy_flight.dat',status='replace')

open(2,file='levy_flight_plot.plt',status='replace')

nsteps=1000

x0=0
y0=0
z0=0

do i=1,nsteps


call random_number(p)
call random_number(q)
call cauchy_sample(step)

theta=2*pi*(1-p)
phi=2*pi*(1-q)

 x0=x0+step*COS(phi)*SIN(theta)
 y0=y0+step*SIN(theta)*SIN(phi)
 z0=z0+step*COS(theta)



write(1,*)x0,y0,z0

end do


write(2,*)"set xlabel 'x'"
write(2,*)"set ylabel 'y'"
write(2,*)"set zlabel 'z'"
write(2,*)"set title 'Levy flight'"
write(2,*)"splot 'levy_flight.dat'  with line lt rgb 'blue' title 'trajectory'" 
call system('gnuplot -p levy_flight_plot.plt')

close(1)
close(2)



contains

subroutine random_stduniform(u)
        implicit none
        real,intent(out) :: u
        real :: r
        call random_number(r)
        u = 1 - r
end subroutine random_stduniform


subroutine random_stdnormal(x)
        implicit none
        real,intent(out) :: x
        real,parameter :: pi=4*ATAN(1.0)
        real :: u1,u2
        call random_stduniform(u1)
        call random_stduniform(u2)
        x = sqrt(-2*log(u1))*cos(2*pi*u2)
end subroutine random_stdnormal


subroutine cauchy_sample(k)
        implicit none
        real,intent(out)::k
	real::n1,n2
	call random_stdnormal(n1)
	call random_stdnormal(n2)
	k=n1/n2
end subroutine cauchy_sample



end program Levy_flight
