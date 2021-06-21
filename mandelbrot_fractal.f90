program Mandelbrot

!This program  is a bare bones implementation  to plot the mandelbrot fratcal using gnuplot to plot the results

implicit none
double precision,parameter::l=4.0000000d0
double precision,parameter::h=2.0000000d0
double precision:: dx,dy,xc,yc,x,y
integer,parameter::nx=1200  
integer,parameter::ny=800
integer,parameter::niter=100
integer::i,j,counter
complex::Z,Znp1,C

xc=l/2
yc=h/2

dx=l/nx
dy=h/ny

Z=0
C=0

open(1,file='xy_mandelbrot.dat',status='replace')
open(2,file='xyplot.plt',status='replace')
         
do i=0,nx
do j=0,ny

x=xc-dx*i
y=yc-dy*j

do counter=1,niter
			
	C=complex(x,y)
			
	Znp1=(Z**2)+C
			
	Z=Znp1
	
end do

if(ABS(Znp1).LE.4.0)then
write(1,*)x,y
end if

z=0

end do
end do

write(2,*)'set title "Mandelbrot Set"'
write(2,*)'set xlabel "Re"'
write(2,*)'set ylabel "Im"'
write(2,*)'set size ratio 1'
write(2,*)'set term png'
write(2,*)'set output "mandelbrot_set.png" '
write(2,*)'plot "xy_mandelbrot.dat" with points pt 0.2 lc rgb "black" title " " '

call system('gnuplot -p xyplot.plt')


end program Mandelbrot