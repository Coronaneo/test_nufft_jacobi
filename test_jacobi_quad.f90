program test_jacobi_quad
use utils
use chebyshev
implicit double precision (a-h,o-z)

double precision, allocatable :: xs(:),xwhts(:)
double precision, allocatable :: ts(:),twhts(:)
double precision, allocatable :: vals(:)

double precision, allocatable :: psivals1(:,:),psivals2(:,:),avals1(:,:),avals2(:,:)
double precision, allocatable :: vals1(:),vals2(:),ab(:,:)
double precision, allocatable :: avals0(:),psivals0(:)
type(chebexps_data)           :: chebdata


pi = acos(-1.0d0)
n  = 2**16+1
da = 0.25d0
db = 0.25d0
allocate(xs(n),xwhts(n))
call elapsed(t1)
call jacobi_quad(n,da,db,xs,xwhts)
call elapsed(t2)

call prini("n = ",n)
call prin2("time = ",t2-t1)
call prin2("time per node = ",(t2-t1)/n)

! iw = 20
! open(iw,FILE='quad.txt')

! do i=1,n
! read(iw,*) x0,t0,wht0
! !if (n/2-10 .lt. i .AND. i .lt. n/2+10) then
! !if (i .gt. n-10) then
! !if (i .lt. 10) then
! !print *,i,xs(i)-x0,(xwhts(i)-wht0)/wht0
! !endif

! end do

! close(iw)


!
!  Make the modified quadrature rule
!
allocate(ts(n),twhts(n))
call jacobi_quad_mod(n,da,db,ts,twhts)
print *,""

! call jacobi_quad_mod2(n,da,db,ts,twhts)
! stop
sum1 = 0
nu1  = 102
nu2  = 102
dnu1 = nu1
dnu2 = nu2

k = 30
call chebexps(k,chebdata)
nints = 50

allocate(ab(2,nints))
allocate(avals1(k,nints),psivals1(k,nints),vals1(n))
allocate(avals2(k,nints),psivals2(k,nints),vals2(n))
allocate(avals0(n),psivals0(n))

call jacobi_phase_disc(nints,ab)
call jacobi_phase(chebdata,dnu1,da,db,nints,ab,avals1,psivals1)
call jacobi_phase(chebdata,dnu2,db,da,nints,ab,avals2,psivals2)

call jacobi_phase_eval(chebdata,dnu1,da,db,nints,ab,avals1,psivals1,n,ts,avals0,psivals0)
vals1 = cos(psivals0)*avals0

call jacobi_phase_eval(chebdata,dnu2,da,db,nints,ab,avals2,psivals2,n,ts,avals0,psivals0)
vals2 = cos(psivals0)*avals0


vals1 = vals1 * sqrt(twhts)
vals2 = vals2 * sqrt(twhts)

print *,sum(vals1*vals2)


end program
