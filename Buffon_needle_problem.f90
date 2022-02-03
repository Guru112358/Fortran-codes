program needle
use omp_lib
implicit none 

integer,parameter::ndrops=10000000

real(8),parameter::pi=4*ATAN(1.0000000d0)

real(8)::y,theta,space,pi_estimate

integer::i,hits

hits=0

!$omp parallel do private(y,theta,space) reduction(+:hits)

do i=1,ndrops

call random_uniform(0.00000d0,pi,theta)
call random_uniform(0.000d0,1.0000d0,y)

space=0.5d0*sin(theta)

if((y+space>1.0000d0).OR.(y-space<0))then

hits=hits+1

end if

end do

pi_estimate=(2*ndrops)/real(hits)


write(*,*)"The value of Pi estimated by",ndrops,"needles is",pi_estimate
write(*,*)'Error of estimate is',abs(pi-pi_estimate)*100,"percent"


end program

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

