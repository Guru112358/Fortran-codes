program rwtest 

implicit none

call random_walk()

contains

subroutine random_walk()
implicit none

integer::i,nsteps,ntrials,t
real,parameter::pi=4*ATAN(1.0)
real::phi,theta,k,k1,k2,u1,u2,p
real::x0,y0,z0,r,step,mean_x,mean_y,mean_z
real,allocatable,dimension(:)::x,y,z,rsquared
real,allocatable,dimension(:):: sigmax,sigmay,sigmaz,distance
real::sumx,sumy,sumz,rbar,time1,time2


write(*,*)"Enter the number of steps for a random walk: "
write(*,*)
read *,nsteps
write(*,*)
write(*,*)"Enter the number of trials: "
write(*,*)
read *,ntrials

allocate(sigmax(ntrials),sigmay(ntrials),sigmaz(ntrials))
allocate(distance(ntrials),rsquared(ntrials))
allocate(x(nsteps),y(nsteps),z(nsteps))

write(*,*)
write(*,*)"Beginning random walk simulation of",nsteps,"steps and",ntrials,"trials" 
write(*,*)
write(*,*)"please wait...."
write(*,*)

call CPU_TIME(time1)

do t=1,ntrials

x0=0
y0=0
z0=0

do i =1,nsteps

    call RANDOM_NUMBER(k)
    call RANDOM_NUMBER(k1)
    call RANDOM_NUMBER(k2)

    call RANDOM_NUMBER(p)

    u1=1-k1
    u2=1-k2
    !step=1
    step = SQRT(-2*LOG(u1))*COS(2*pi*u2)
    theta=2*pi*(1-k)
    phi=2*pi*(1-p)

        x0=x0+step*COS(phi)*SIN(theta)
        y0=y0+step*SIN(theta)*SIN(phi)
        z0=z0+step*COS(theta)
        
        x(i)=x0
        y(i)=y0
        z(i)=z0

end do

distance(t)=SQRT(x0**2+y0**2+z0**2)

rsquared(t)=distance(t)**2

mean_x=SUM(x)/nsteps
mean_y=SUM(y)/nsteps
mean_z=SUM(z)/nsteps

sigmax(t)=mean_x
sigmay(t)=mean_y
sigmaz(t)=mean_z

end do

sumx=SUM(sigmax)/ntrials
sumx=SUM(sigmay)/ntrials
sumx=SUM(sigmaz)/ntrials
rbar=SQRT((SUM(rsquared))/(ntrials))

call CPU_TIME(time2)


write(*,*)"Simulation complete!!"
write(*,*)
write(*,*)"Time taken:",time2-time1,"seconds"
write(*,*)
write(*,*)"+++++++++++++++++++++++++++++++"
write(*,*)"!-----------RESULTS-----------!"
write(*,*)"+++++++++++++++++++++++++++++++"
write(*,*)
write(*,*)
write(*,*)"Average x distance travelled by the particle:"
write(*,*)sumx
write(*,*)
write(*,*)"Average y distance travelled by the particle:"
write(*,*)sumy
write(*,*)
write(*,*)"Average z distance travelled by the particle:"
write(*,*)sumz
write(*,*)
write(*,*)"Average distance from origin travelled by the particle:"
write(*,*)rbar
write(*,*)

deallocate(sigmax,sigmay,sigmaz,distance)
deallocate(x,y,z,rsquared)

end subroutine random_walk


end program rwtest



















