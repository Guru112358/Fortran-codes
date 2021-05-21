program Newton_Raphson
IMPLICIT NONE    


call root_finder()


contains


real function f(x)
implicit none
real::x
f=x**3-5
return
end function 


subroutine root_finder()
!This Subroutine finds the roots of the function f(x) using a Newton Raphson iterative method and approximating the first derivative using the secant method.
!may not always work,has many limitations,user be warned,see Wikipedia of examples of stuff that does not work.
real::x0,x1,xnp1,residual,tol
integer::niter,i

write(*,*)"Enter the two initial guesses x0 and x1: "
write(*,*)
read *,x0,x1
write(*,*)
write(*,*)"enter the number of max iterations: "
write(*,*)
read *,niter
write(*,*)
write(*,*)"enter the desired order of accuracy: "
write(*,*)
read *,tol

write(*,*)
write(*,*)"Solution loop starting.......... "
write(*,*)

do i=1,niter
    xnp1=x1-f(x1)*((x1-x0)/(f(x1)-f(x0)))
    x0=x1
    x1=xnp1
    write(*,*)i,xnp1
    residual=abs(f(xnp1))

    if(residual<tol)then
    write(*,*)"solution successfully converged to a tolerance of",tol
    write(*,*)
    write(*,*)"The root is x=",xnp1
    write(*,*)
    EXIT

    else if(i==niter)then
    write(*,*) 
    write(*,*)i,"WARNING: max iterations reached!"
    write(*,*)
    end if
end do
return

end subroutine root_finder





end program Newton_Raphson
