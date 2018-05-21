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
  mwSize       :: nj

type(chebexps_data)           :: chebdata
double precision, allocatable :: ab(:,:)
double precision, allocatable :: psivals(:),avals(:)
double precision, allocatable :: ts(:),avals0(:),psivals0(:),polvals(:),polvals0(:)
double precision, allocatable :: xx(:),xs(:),twhts(:)
complex*16,allocatable :: jacobi1(:,:),jacobi2(:,:)
integer*4 nts,k
double precision pi

pi  = acos(-1.0d0)
!call mxCopyPtrToInteger4(mxGetPr(prhs(1)),nts,1)
nts = 2**9
allocate(jacobi1(nts,nts-27),jacobi2(nts,nts-27))
allocate(ts(nts),xs(nts),twhts(nts),avals0(nts),psivals0(nts),polvals(nts),polvals0(nts))				  

da=0.25d0
db=0.25d0
call jacobi_quad_mod(nts,da,db,ts,twhts)

do i=1,nts
xs(i) = (mod(floor(ts(i)/2/pi*nts+0.5),nts))*2*pi/nts
end do
call quicksort(nts,ts)
call quicksort(nts,xs)


do i = 27,nts-1

k  = 24
call chebexps(k,chebdata)
if (i.ge.28) then
   deallocate(xx,psivals,avals,ab)
end if


! this must be integer because we use the recurrence relations to test accuracy
nu  = i

dnu = nu
da  = 0.25d0
db  = 0.25d0
p   = dnu + (da+db+1)/2
allocate(xx(0:nu))

dd      = nu
dd      = 1/(dd) * 2/pi
dd      = log(dd)/log(2.0d0)
nints   = ceiling(-dd)*2
nints   = max(nints,20)
allocate(ab(2,nints))

call jacobi_phase_disc(nints,ab)
!call prini("nints = ",nints)
!call prin2("before jacobi_phase, ab = ",ab)

allocate(psivals(k*nints),avals(k*nints))

call jacobi_phase(chebdata,dnu,da,db,nints,ab,avals,psivals)
!call prin2("time to construct phase = ",(t2-t1))


call jacobi_phase_eval(chebdata,dnu,da,db,nints,ab,avals,psivals,nts,ts,avals0,psivals0)
!call prin2("average eval time = ",(t2-t1)/nts)
!polvals0 = cos(psivals0)*avals0
!print *,size(avals0),size(psivals0)
jacobi1(:,i-26)=avals0*exp(dcmplx(0,1)*(psivals0-i*ts))
jacobi2(:,i-26)=avals0*exp(dcmplx(0,1)*(psivals0-i*xs))
end do

plhs(1) = mxCreateDoubleMatrix(nts, 1, 0)
plhs(2) = mxCreateDoubleMatrix(nts, nts-27, 1)
plhs(3) = mxCreateDoubleMatrix(nts, nts-27, 1)

call mxCopyReal8ToPtr(ts, mxGetPr(plhs(1)),nts) 
call mxCopyComplex16ToPtr(jacobi1, mxGetPr(plhs(2)),mxGetPi(plhs(2)),nts*(nts-27)) 
call mxCopyComplex16ToPtr(jacobi2, mxGetPr(plhs(3)),mxGetPi(plhs(3)),nts*(nts-27)) 

deallocate(jacobi1,jacobi2,ts,xs,twhts,avals0,psivals0,polvals,polvals0)
end subroutine mexfunction
