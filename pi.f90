program monte
implicit none
real::x,y,r,pi
integer::n,i,c
c=0 !counter
n=10000000
open(1,file="pi12.dat",status="old")
do i=1,n
call random_number(x)
call random_number(y)
r=(x**2+y**2)**0.5
if(r<=1)then
c=c+1
pi=4.0*real(c)/(i)
print *,i,pi
write(1,*)i,x,y,r,pi
end if
end do
close(1)
end program
