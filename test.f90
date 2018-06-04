program test
implicit none
integer n
integer,allocatable :: t(:),w(:,:)
n=10
print *,n
allocate(t(n))
allocate(w(n,n))
t(1)=10
t(2)=11
t(3:9)=19
t(10)=20
print *,mod(t,n)
end program 
