program test
implicit none
integer n
real*8,allocatable :: t(:,:),w(:,:)
n=2**17
print *,n
allocate(t(n,n))
allocate(w(n,n))
t=1
w=1
print *,t(n,n)
end program 
