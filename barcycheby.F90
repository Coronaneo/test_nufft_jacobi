
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
mwSize       :: nt,nx,nw

real*8, allocatable :: t(:),x(:),w(:),S(:,:),p(:),c(:)
real*8, allocatable :: medi(:)
integer*4 i,j

nt = mxGetM(prhs(1))
nx = mxGetM(prhs(2))
nw = mxGetM(prhs(3))

allocate(x(nx),t(nt),w(nw),S(nt,nx),p(nt),c(nx),medi(nx))

call mxCopyPtrToReal8(mxGetPr(prhs(1)),x,nx)
call mxCopyPtrToReal8(mxGetPr(prhs(2)),t,nt)
call mxCopyPtrToReal8(mxGetPr(prhs(3)),w,nw)

c = 1.0d0
do i = 1,nt
   medi = c*t(i)-x
   medi = 1.0d0/medi
   p(i) = sum(w*medi)
enddo

do i = 1,nt
   do j = 1,nx
      S(i,j) = w(j)/(t(i)-x(j))/p(i)
   enddo
enddo

plhs(1) = mxCreateDoubleMatrix(nt, nx, 0)
call mxCopyReal8ToPtr(S, mxGetPr(plhs(1)),nt*nx)

deallocate(x,t,w,S,p,c,medi)
end subroutine
