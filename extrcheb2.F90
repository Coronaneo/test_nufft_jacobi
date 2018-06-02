#include "fintrf.h"

! standard header for mex functions
subroutine mexfunction(nlhs, plhs, nrhs, prhs)
use utils
use chebyshev
implicit double precision (a-h,o-z)

! number of input arguments, number of output arguments
integer      :: nlhs, nrhs  
! pointer to inputs and outputs
mwPointer    :: plhs(*), prhs(*) 
! get some of the matlab mex functions
mwPointer    :: mxGetPr, mxGetPi, mxCreateDoubleMatrix 
! define a size integer so that we can get its type
mwSize       :: n,n1,n2

type(chebexps_data)           :: chebdata
real*8, allocatable :: twhts(:),ab(:,:),t(:),k(:)
complex*16, allocatable :: r(:),m(:,:)
real*8, allocatable :: psivals(:),avals(:)
real*8, allocatable :: ts(:),avals0(:),psivals0(:),xs(:)
real*8, allocatable :: avals1(:),psivals1(:)
integer, allocatable :: t1(:),k1(:)
integer*4 it,kk,a
real*8 da,db,flag

da=-0.5d0
db=-0.5d0

n = mxGetM(prhs(1))
if (n .lt. (2**12-1)) then
it = 9
else 
it = 27
end if

n1 = mxGetM(prhs(2))
n2 = mxGetM(prhs(3))
allocate(t(n1),k(n2),t1(n1),k1(n2),m(n1,n2),r(n**2))
call mxCopyPtrToReal8(mxGetPr(prhs(2)),t,n1)
call mxCopyPtrToReal8(mxGetPr(prhs(3)),k,n2)
call mxCopyPtrToReal8(mxGetPr(prhs(4)),flag,1)
t1=int(t)
k1=int(k)

allocate(ts(n),twhts(n),xs(n))
allocate(avals0(n),psivals0(n),avals1(n),psivals1(n))

call jacobi_quad_mod(n,da,db,ts,twhts)
xs = int(ts/2/pi*(n+0.0d0))*2*pi/n

kk  = 16
call chebexps(kk,chebdata)

dd      = nts*1.0d0
dd      = min(0.10d0,1/dd)
dd      = log(dd)/log(2.0d0)
nints   = ceiling(-dd)+1
nints   = 2*nints
allocate(ab(2,nints))

call jacobi_phase_disc(nints,ab)

allocate(psivals(kk*nints),avals(kk*nints))

m=0
do i=1,n2
dnu = mod(k1(i)+it*n,n)
if (dnu .ge. it) then
    r=0
    j = int(k1(i)/n)+it
    call jacobi_phase(chebdata,j,da,db,nints,ab,avals,psivals)
    call jacobi_phase_eval(chebdata,j,da,db,nints,ab,avals,psivals,n,ts,avals0,psivals0)
    call jacobi_phase(chebdata,dnu,da,db,nints,ab,avals,psivals)
    call jacobi_phase_eval(chebdata,dnu,da,db,nints,ab,avals,psivals,n,ts,avals1,psivals1)
    do a=1,n
       r((a-1)*n+1:a*n)=avals0(a)*exp(dcmplx(0,1)*(psivals0(a)-j*ts(a)-j*flag*(xs(a)-ts(a))))*avals1*exp(dcmplx(0,1)   &
        *(psivals1-dnu*ts-dnu*flag*(xs-ts)))+r((a-1)*n+1:a*n)
    end do
    m(:,i)=r(t1)
end if
end do

plhs(1)=mxCreateDoubleMatrix(n1, n2, 1)
call mxCopyComplex16ToPtr(m, mxGetPr(plhs(1)),mxGetPi(plhs(1)),n1*n2)

deallocate(twhts,ab,t,k,t1,k1,m,r,ts,xs,avals,psivals,avals0,psivals0,avals1,psivals1)

end subroutine
