program test_jacobi_transform

use utils
use chebyshev
use idecomp
use jacobi_exp
use jacobi_transform

implicit double precision (a-h,o-z)

double precision, allocatable :: vals(:,:),ts(:),whts(:)

type(jacobi_expansion_data)   :: expdata
type(jacobi_transform_data)   :: jacdata
double precision, allocatable :: x(:),y(:),x0(:),y0(:)

double precision, allocatable :: times(:,:),errs(:),dsizes(:)
integer, allocatable          :: kranks(:)


pi       = acos(-1.0d0)
ima      = (0.0d0,1.0d0)

da       = -0.25d0
db       =  0.00d0

i1       = 5
i2       = 20

allocate(times(3,i1:i2),errs(i1:i2),dsizes(i1:i2),kranks(i1:i2))
times =0 
errs  = 0
dsizes = 0
kranks  =0 

do i=i1,i2

eps      = 1.0d-12/2
n        = 2**i
dmax     = n
iffactor = 1
ifeval   = 0

call prini("n = ",n)
call elapsed(t1)
call jacobi_expansion(eps,iffactor,dmax,da,db,expdata)
call jacobi_transform_prepare(expdata,n,jacdata)
call elapsed(t2)


times(1,i) = t2-t1
kranks(i)  = expdata%krank

call jacobi_expansion_size(expdata,dsize1)
call jacobi_transform_size(jacdata,dsize2)
dsize = dsize1+dsize2
call prin2("jacobi data size (MB) = ",dsize)
dsizes(i) = dsize

call prini("krank = ",expdata%krank)
call prin2("prepare time = ",t2-t1)


! construct an input vector and allocate memory for the output
allocate(x(n),y0(n),y(n),x0(n))
x      = 0

do j=1,n
call random_number(x(j))
end do
x = x / norm2(x)

call elapsed(t1)
call jacobi_transform_forward(jacdata,x,y)
call elapsed(t2)
call prin2("forward apply time = ",t2-t1)
times(2,i) = t2-t1

! call jacobi_transform_forward_bf(da,db,n,n,jacdata%ts,jacdata%whts,x,y)


call elapsed(t1)
call jacobi_transform_backward(jacdata,y,x0)
call elapsed(t2)
call prin2("backward  apply time = ",t2-t1)
times(3,i) = t2-t1

! call jacobi_transform_backward_bf(da,db,n,n,jacdata%ts,jacdata%whts,y,x0)


print *,"composition error           =",maxval(abs(x-x0))
errs(i) = maxval(abs(x-x0))


!times(3,i) = t2-t1


deallocate(x,y,x0,y0)
call jacobi_transform_destroy(jacdata)
end do

print *,""
print *,""


! print *,"n     krank    prepare time   forward time    backward time     error        data size"
! print *,"---------------------------------------------------------------------------------------"
! do i=i1,i2
! write (*,"(1X,I2.2,5X,I3.3,5X,D8.3,8X,D8.3,8X,D8.3,8X,D8.3,7X,D8.3)")  &
!   i,kranks(i),times(1,i),times(2,i),times(3,i),errs(i),dsizes(i)
!end do

print *,"n     krank    prepare time   forward time    backward time     error   "
print *,"---------------------------------------------------------------------------------------"
do i=i1,i2
write (*,"(1X,I2.2,5X,I3.3,5X,D8.3,8X,D8.3,8X,D8.3,8X,D8.3)")  &
  i,kranks(i),times(1,i),times(2,i),times(3,i),errs(i)

end do
print *,""
print *,""




end program



