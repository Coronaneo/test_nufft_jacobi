program test_jacobi_exp
use utils
use chebyshev
use jacobi_exp
implicit double precision (a-h,o-z)
type(jacobi_expansion_data)  :: expdata, expdata2

!
!  Build and test an expansion
!

pi       = acos(-1.0d0)
da       =  0.00d0
db       =  0.25d0
eps      = 1.0d-12
iffactor = 1
dmax     = 2**16

call elapsed(t1)
call jacobi_expansion(eps,iffactor,dmax,da,db,expdata)
call elapsed(t2)
call prin2("expansion time = ",t2-t1)

call jacobi_expansion_size(expdata,dsize)
call prin2("expansion size (MB) = ",dsize)

call prini("n = ",expdata%n)
call prini("m = ",expdata%m)
call prini("krank = ",expdata%krank)

call jacobi_expansion_test(expdata)

! test the reading and writing of the expansion

! iw = 20
! open(iw,FILE='jacobi_expansion.txt')
! call jacobi_expansion_write(iw,expdata)
! close(iw)

! iw = 20
! open(iw,FILE='jacobi_expansion.txt',STATUS='OLD')
! call jacobi_expansion_read(iw,expdata2)
! close(iw)

! call jacobi_expansion_test(expdata2)

end program
