program test_idecom
use utils
use idecomp

implicit double precision (a-h,o-z)

! double precision, allocatable :: a(:,:),a0(:,:),r(:,:),b(:,:),r2(:,:)
! integer, allocatable          :: ipivs(:)

! n = 100
! m = 200
! allocate(a(n,m),a0(n,m))

! do i=1,n
! do j=1,m
! a(i,j) = cos(i*j/1000.0d0)
! end do
! end do
! a0 = a
! eps = 1.0d-15
! call idecomp_construct(eps,a0,krank,ipivs,r)

! deallocate(a0)
! allocate(a0(n,krank))
! a0 = a(:,ipivs(1:krank))

! print *,norm2(a - matmul(a0,r))


double complex, allocatable   :: a(:,:),a0(:,:),r(:,:),b(:,:),r2(:,:)
integer, allocatable          :: ipivs(:)
double complex                :: ima

ima = (0.0d0,1.0d0)
n = 100
m = 200
allocate(a(n,m),a0(n,m))

do i=1,n
do j=1,m 
a(i,j) = cos(i*j/100.0d0) + ima * 17*sin(i*j/50.0d0)
end do
end do
a0 = a
eps = 1.0d-14
call cidecomp_construct(eps,a0,krank,ipivs,r)


deallocate(a0)
allocate(a0(n,krank))
a0 = a(:,ipivs(1:krank))

print *,norm2(abs(a - matmul(a0,r)))

end program
