program monte_carlo_integration_parallel
use omp_lib

implicit none

integer,parameter::npoints=100000000
integer::i
real(8)::a,b,x,sum,I_term,result

a=0.0
b=5.0
sum=0;

!$omp parallel do private(x,I_term) reduction(+:sum)

do i =1,npoints

call random_uniform(a,b,x)
I_term=f(x)
sum=sum+I_term
end do
!$omp end parallel do

result=((b-a)/npoints)*sum

write(*,*)'the value of the integral is: ',result



contains

real function f(x) result(func)
real(8)::x
func=(x**4)*EXP(-x)


end function 


subroutine random_stduniform(u)
   implicit none
   real(8),intent(out) :: u
   real(8) :: r
   call random_number(r)
   u = 1 - r
end subroutine random_stduniform



subroutine random_uniform(a,b,x)
   implicit none
   real(8),intent(in) :: a,b
   real(8),intent(out) :: x
   real(8) :: u
   call random_stduniform(u)
   x = (b-a)*u + a
end subroutine random_uniform



end program monte_carlo_integration_parallel
