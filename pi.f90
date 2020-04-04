program monte
!this program computes the value of Pi approximately by employing a Monte carlo method.
!the random__number module is called to produce random points between 0 and 1 in x and y coordinates of a square region of space 1 X 1.
!if they fall within a unit circle,the point is counted.
!increasing number of random points increases accuracy.
!this method is limted by the quality of the pseudorandom number generator and number of points the memory can handle.
implicit none
real::x,y,r,pi
integer::n,i,c
c=0 !counter
n=10000!number of random points
open(1,file="pi12.dat",status="old")!data file
open(2,file="montepi.plt",status="old")!plot file
do i=1,n
call random_number(x)
call random_number(y)
r=(x**2+y**2)**0.5
if(r<=1)then
c=c+1 !incrementing each time the random point falls within the circle
pi=4.0*real(c)/(i)
print *,i,x,y,r,pi
write(1,*)x,y
end if
end do
write(2,*)'set xlabel "X"'
write(2,*)'set ylabel "Y"'
write(2,*)'plot "pi12.dat"'
CALL SYSTEM('gnuplot -p montepi.plt') !calling gnuplot using the call system command

close(1)
close(2)
end program
