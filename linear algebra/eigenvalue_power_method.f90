program  power_method
implicit none


real,dimension(5,5) :: a = reshape((/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25/),(/5,5/))


call power_iteration(a,100,0.00000001)


contains


subroutine power_iteration(A,niter,tol)  !A is a square n x n matrix,niter is the max number of iterations,tol is the error tolerance before exiting
!This sunroutine implements a scheme to iteratively find the largest eigenvalue and its associated eigenvector for a n x n matrix A
!This method is also attributed to Richard Von Mises
implicit none
real,dimension(:,:)::A
real,allocatable,dimension(:)::bnp1,b0,indices,eigenvector
integer:: niter,imax,jmax,i,j,count
real::lambdamax,p,tol,res


indices=shape(A)
imax=indices(1)
jmax=indices(2)

allocate(b0(imax),bnp1(imax))

write(*,*)"Beginning power iterations......"
write(*,*)

lambdamax=0
!using a vector of random elements to intitialsie the procedure

do i=1,imax
    call random_number(p)
    b0(i)=p
end do

do count=1,niter 

    bnp1=MATMUL(A,b0)/(NORM2(MATMUL(A,b0)))
    eigenvector=(1/MAXVAL(bnp1))*bnp1
    lambdamax=DOT_PRODUCT(bnp1,MATMUL(A,bnp1))/(DOT_PRODUCT(bnp1,bnp1))
    res= ABS(MAXVAL(b0-bnp1))
    b0=bnp1
    write(*,*)count,res
    write(*,*)

    if(res.LE.tol)then
    write(*,*)"solution converged to a tolerance of",tol,"after",count,"Iterations"
    write(*,*)
    write(*,*)"Exiting the interations loop..."
    EXIT
    else if(count.EQ.niter)then
    write(*,*)"MAX iterations reached!!"
    end if

end do

write(*,*)
write(*,*)
write(*,*)'The largest Eigenvalue is: ',lambdamax
write(*,*)
write(*,*)
write(*,*)'The associated Eigenvector is:  [',eigenvector,']'

deallocate(b0,bnp1)


end subroutine  power_iteration



end program power_method
