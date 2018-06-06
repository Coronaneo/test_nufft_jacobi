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
complex*16, allocatable :: r(:),rr(:),ier(:),rrr(:)
real*8, allocatable :: psivals(:),avals(:),rd(:),rd1(:),rd2(:)
real*8, allocatable :: ts(:),avals0(:),psivals0(:)
real*8, allocatable :: psival(:,:),aval(:,:),avals1(:),psivals1(:)
integer*4 k,nn
integer*4 it,i,j
real*8 da,db
complex*16 a

!da=-0.50d0
!db=-0.50d0
allocate(ier(5))
ier=0
n = mxGetM(prhs(1))
nn = int(log(real(nn))/log(2.0d0)+0.01)

if (n .lt. (2**12-1)) then
it = 27
else 
it = 9
end if

allocate(c((n-it)**2),r(n**2),rr(n),rrr(n),ts(n),twhts(n))
allocate(avals0(n),psivals0(n),avals1(n),psivals1(n),rd(nn),rd1(nn),rd2(nn))

call mxCopyPtrToReal8(mxGetPr(prhs(2)),c,(n-it)**2)
call mxCopyPtrToReal8(mxGetPr(prhs(3)),da,1)
call mxCopyPtrToReal8(mxGetPr(prhs(4)),db,1)
call mxCopyPtrToReal8(mxGetPr(prhs(5)),rd,nn)
rd1 = int(rd+0.01)
call mxCopyPtrToReal8(mxGetPr(prhs(6)),rd,nn)
rd2 = int(rd+0.01)

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
!ier(1)=aval(1,2)
!ier(2)=psival(1,2)
r=0
do i=1,nn
dnu = rd1(i)
call jacobi_phase_eval(chebdata,dnu,da,db,nints,ab,aval(:,rd1(i)-it+1),psival(:,rd1(i)-it+1),n,ts,avals0,psivals0)
rrr = avals0*exp(dcmplx(0,1)*psivals0)
   do j=1,nn
      dnu1 = rd2(j)
      call jacobi_phase_eval(chebdata,dnu1,da,db,nints,ab,aval(:,rd2(j)-it+1),psival(:,rd2(j)-it+1),n,ts,avals1,psivals1)

      rr=avals1*exp(dcmplx(0,1)*psivals1)
      do k=1,n
         r((k-1)*n+1:k*n)=c(rd2(j)-it+1+(rd1(i)-it)*(n-it))*rrr(k)*rr+r((k-1)*n+1:k*n)
      end do
   end do
end do
!ier(4)=1
!ier=c(1:5)
plhs(1) = mxCreateDoubleMatrix(n**2, 1, 1)
plhs(2) = mxCreateDoubleMatrix(5,1,1)
plhs(3) = mxCreateDoubleMatrix(n,1,0)
call mxCopyComplex16ToPtr(r, mxGetPr(plhs(1)),mxGetPi(plhs(1)),n**2)
!ier(5)=1
call mxCopyComplex16ToPtr(ier,mxGetPr(plhs(2)),mxGetPi(plhs(2)),5)
call mxCopyReal8ToPtr(ts,mxGetPr(plhs(3)),n)
deallocate(c,twhts,ts,ab,r,psivals,avals,psivals0,avals0,psival,aval,avals1,psivals1)


end subroutine
