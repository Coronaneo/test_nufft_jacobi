
#include "fintrf.h"

! standard header for mex functions
subroutine mexfunction(nlhs, plhs, nrhs, prhs)
use utils
use chebyshev
implicit double precision (a-j,o-z)

! number of input arguments, number of output arguments
integer      :: nlhs, nrhs  
! pointer to inputs and outputs
mwPointer    :: plhs(*), prhs(*) 
! get some of the matlab mex functions
mwPointer    :: mxGetPr, mxGetPi, mxCreateDoubleMatrix 
! define a size integer so that we can get its type
mwSize       :: n

type(chebexps_data)           :: chebdata
real*8, allocatable :: c(:),twhts(:),ab(:,:)
complex*16, allocatable :: r(:),ier(:)
real*8, allocatable :: psivals(:),avals(:)
real*8, allocatable :: ts(:),avals0(:),psivals0(:)
real*8, allocatable :: psival(:,:),aval(:,:)
integer*4 k,ii,jj,kk
integer*4 it,i,j
real*8 da,db
complex*16 a


allocate(ier(5))
ier=0
n = mxGetM(prhs(1))

if (n .lt. (2**12-1)) then
it = 27
else 
it = 9
end if

allocate(c(n-it),r(n),ts(n),twhts(n))
allocate(avals0(n),psivals0(n))

call mxCopyPtrToReal8(mxGetPr(prhs(2)),c,n-it)
call mxCopyPtrToReal8(mxGetPr(prhs(3)),da,1)
call mxCopyPtrToReal8(mxGetPr(prhs(4)),db,1)

call jacobi_quad_mod(n,da,db,ts,twhts)
!ier(1)=1
k  = 16
call chebexps(k,chebdata)

dd      = n*1.0d0
dd      = min(0.10d0,1/dd)
dd      = log(dd)/log(2.0d0)
nints   = ceiling(-dd)+1
nints   = 2*nints
allocate(ab(2,nints))

call jacobi_phase_disc(nints,ab)
!ier(2)=1
allocate(psivals(k*nints),avals(k*nints))
allocate(psival(k*nints,n-it),aval(k*nints,n-it))

do i=it,n-1
dnu = i
call jacobi_phase(chebdata,dnu,da,db,nints,ab,avals,psivals)
psival(:,i-it+1) = psivals
aval(:,i-it+1) = avals
end do

r=0
do i=it,n-1
dnu = i
call jacobi_phase_eval(chebdata,dnu,da,db,nints,ab,aval(:,i-it+1),psival(:,i-it+1),n,ts,avals0,psivals0)
r = avals0*exp(dcmplx(0,1)*psivals0)*c(i-it+1)+r
end do
!ier(4)=1
!ier=c(1:5)
plhs(1) = mxCreateDoubleMatrix(n, 1, 1)
plhs(2) = mxCreateDoubleMatrix(5,1,1)
plhs(3) = mxCreateDoubleMatrix(n,1,0)
call mxCopyComplex16ToPtr(r, mxGetPr(plhs(1)),mxGetPi(plhs(1)),n)
!ier(5)=1
call mxCopyComplex16ToPtr(ier,mxGetPr(plhs(2)),mxGetPi(plhs(2)),5)
call mxCopyReal8ToPtr(ts,mxGetPr(plhs(3)),n)
deallocate(c,twhts,ts,ab,r,psivals,avals,psivals0,avals0,psival,aval)


end subroutine
