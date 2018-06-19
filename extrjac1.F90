#include "fintrf.h"

! standard header for mex functions
subroutine mexfunction(nlhs, plhs, nrhs, prhs)
use utils
use chebyshev
implicit double precision (a-k,o-z)

! number of input arguments, number of output arguments
integer      :: nlhs, nrhs  
! pointer to inputs and outputs
mwPointer    :: plhs(*), prhs(*) 
! get some of the matlab mex functions
mwPointer    :: mxGetPr, mxGetPi, mxCreateDoubleMatrix 
! define a size integer so that we can get its type
mwSize       :: n,n1,n2,m1,m2,nts,nnu

type(chebexps_data)           :: chebdata
real*8, allocatable :: twhts(:),ab(:,:),t(:),k(:)
complex*16, allocatable :: r(:),m(:,:),ier(:)
real*8, allocatable :: psivals(:),avals(:),ts(:),nu(:)
real*8, allocatable :: ts(:),avals0(:),psivals0(:),xs(:)
integer*4, allocatable :: t1(:),k1(:)
integer*4 it,kk,ii,i,jj
real*8 da,db,flag,pi


pi=acos(-1.0d0)


allocate(ier(5))
n = mxGetM(prhs(1))
if (n .lt. (2**12-1)) then
it = 27
else 
it = 9
end if

m1 = mxGetN(prhs(2))
m2 = mxGetN(prhs(3))
n1 = mxGetM(prhs(2))
n2 = mxGetM(prhs(3))
n1 = n1*m1
n2 = n2*m2
nts = mxGetM(prhs(7))
nnu = mxGetM(prhs(8))

allocate(t(n1),k(n2),t1(n1),k1(n2),m(n1,n2),r(n1),ts(nts),nu(nnu))
call mxCopyPtrToReal8(mxGetPr(prhs(2)),t,n1)
call mxCopyPtrToReal8(mxGetPr(prhs(3)),k,n2)
call mxCopyPtrToReal8(mxGetPr(prhs(4)),flag,1)
call mxCopyPtrToReal8(mxGetPr(prhs(5)),da,1)
call mxCopyPtrToReal8(mxGetPr(prhs(6)),db,1)
call mxCopyPtrToReal8(mxGetPr(prhs(7)),ts,nts)
call mxCopyPtrToReal8(mxGetPr(prhs(8)),nu,nnu)
t1=int(t+0.5)
k1=int(k+0.5)
ier(1)=flag
allocate(xs(nts))
allocate(avals0(n1),psivals0(n1))
ier(1:5)=k1(1:5)
!ier(3:5)=k1(1:3)
!call jacobi_quad_mod(n,da,db,ts,twhts)
xs = int(ts/2/pi*(n+0.0d0))*2*pi/n

kk  = 16
call chebexps(kk,chebdata)

dd      = n*1.0d0
dd      = min(0.10d0,1/dd)
dd      = log(dd)/log(2.0d0)
nints   = ceiling(-dd)+1
nints   = 2*nints
allocate(ab(2,nints))

call jacobi_phase_disc(nints,ab)

allocate(psivals(kk*nints),avals(kk*nints))


m=0
do i=1,n2
if (k1(i) .gt. it) then
    dnu = nu(k1(i)-it)
    call jacobi_phase(chebdata,dnu,da,db,nints,ab,avals,psivals)
    call jacobi_phase_eval(chebdata,dnu,da,db,nints,ab,avals,psivals,n,ts(t1),avals0,psivals0)

    if (flag .lt. 0) then
       r = avals0*exp(dcmplx(0,1)*(psivals0-dnu*ts(t1)))
    else
       r = avals0*exp(dcmplx(0,1)*(psivals0-dnu*xs(t1)))
    end if
    
    m(:,i) = r
end if
end do

plhs(1)=mxCreateDoubleMatrix(n1, n2, 1)
plhs(2)=mxCreateDoubleMatrix(5,1,1)
call mxCopyComplex16ToPtr(m, mxGetPr(plhs(1)),mxGetPi(plhs(1)),n1*n2)
call mxCopyComplex16ToPtr(ier,mxGetPr(plhs(2)),mxGetPi(plhs(2)),5)
deallocate(ab,t,k,t1,k1,m,r,ts,nu,xs,avals,psivals,avals0,psivals0)

end subroutine
