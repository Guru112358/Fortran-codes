program RK4secondorder
    implicit none
    !solving a ordinary differential equation of the form y''=f(x,y,y')
    !y'=dy/dx=z,dz/dx=f(x,y,z)
    !The current problem is taken for validation from "higher engineering mathematics" by B.S Grewal(43rd edition,page 1033)
    real::x,y,h,k1,k2,k3,k4,f,xp,l1,l2,l3,l4,l,k,z,phi
    integer::i,n
    write(*,*)"enter values of x and y"
    read *,x,y
    write(*,*)"input value of x at which y is required"
    read *,xp
    write(*,*)"enter the value of z/(dy/dx) at x and y"
    read *,z
    write(*,*)"enter h"
    read *,h
    n=int((xp-x)/h)
    do i=1,n
        l1=h*phi(x,y,z)
        l2=h*phi((x+(h/2)),(y+(k1/2)),(z+l1/2))
        l3=h*phi((x+(h/2)),(y+(k2/2)),(z+l2/2))
        l4=h*phi((x+h),(y+k3),(z+l3))
        l=((1.0/6.0)*(l1+(2*l2)+(2*l3)+l4))
        k1=h*f(z)
        k2=h*f(z+l1/2)
        k3=h*f(z+l2/2)
        k4=h*f(z+l3)
        k=((1.0/6.0)*(k1+(2*k2)+(2*k3)+k4))
        x=x+h
        y=y+k
        z=z+l
        write(*,*)x,y,z
    end do
end program
!note:It is also possible to create one more general function for performing the numerical action of this method and call it twice for computing l and k.
!There are two functions needed,one for each simultaneous coupled Ordinary Differential Equations
real function f(z)
    real::z
    f=z
end function

real function phi(x,y,z)
    real::x,y,z
    phi=(x*z**2)-y**2

end function
