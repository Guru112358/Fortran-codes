program Thomas_algorithm
    !this is a very simple program that solves a tri-diagonal system of equations by the Thomas algorithm by taking the diagonals and RHS matrix as input.
    !this can easily be turned into a subroutine to solve problems in heat transfer,CFD etc
    implicit none
    integer,parameter::n=5
    integer::i
    real,dimension(n)::a,b,c,d,x
    real,dimension(n)::gammar,beta

!!!!!!!!!taking diagonals as input!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

write(*,*) "enter the diagonal elements(b)"
read *,(b(i),i=1,n)
write(*,*)'b= '
write(*,*)(b(i),i=1,n)

write(*,*) "enter the lower diagonal elements(a)"
read *,(a(i),i=1,n)
write(*,*)'a= '
write(*,*)(a(i),i=1,n)

write(*,*) "enter the upper diagonal elements(c)"
read *,(c(i),i=1,n)
write(*,*)'c= '
write(*,*)(c(i),i=1,n)

write(*,*)"enter the RHS matrix(d)"
read *,(d(i),i=1,n)
write(*,*)'d= '
write(*,*)(d(i),i=1,n)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Algorithm !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!forward sweep
 gammar(1) =c(1)/b(1)
 beta(1)= d(1)/b(1)

do i=2,n

    gammar(i)=c(i)/(b(i)-(a(i)*gammar(i-1)))
    beta(i)=(d(i)-a(i)*beta(i-1))/(b(i)-a(i)*gammar(i-1))

end do
!backward sweep
x(n)=beta(n)
do i=(n-1),1,-1
    x(i)=beta(i)-gammar(i)*x(i+1)
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(*,*),'the solution matrix x is: '
write(*,*),(x(i),i=1,n)


end program Thomas_algorithm
