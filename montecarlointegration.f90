program Montecarlointegration
!This program computes the value of an integral by using the Monte Carlo method
!any integral can be be broken up into a summation with infinite terms as the step tends to zero according to the fundamental theorem of Integral calculus.
!the function in the integral sign is repeatedly calculated for a series of random points within the range of the limits.
!for this particular example, integral(exp(x)) is taken as a benchmark.
!a-lower limit
!b-upper limit
implicit none
real::a,b,res,av,sum1,k1
real,dimension(100000)::x,f,k
integer::i,n,count1
n=100000 !number of points
print *,"enter the value of a and b: "
read *,a,b
count1=0
sum1=0.0
k1=0.01
do i=1,n
!building a simple filter
call random_number(k)
!note: the above module only calls pesudorandom numbers between 0 and 1,for other limits,an appropriate scaling factor will have to be introduced.
x(i)=k(i)
if(x(i)>=a)then
if(x(i)<=b)then
count1=count1+1
f(i)=exp(x(i))
sum1=Sum1+f(i)
print *,count1,x(i),f(i)
end if
end if
end do
av=sum1/real(count1)
res=(b-a)*av
print *,"the sum is:  ",sum1
print *,"the average is:  ",av
print *,"the Integral is: ",res
end program

