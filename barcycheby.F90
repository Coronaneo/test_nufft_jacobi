
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

real*8, allocatable :: t(:),x(:),w(:),S(:,:)
integer*4 i,j,temp,flag

nt = mxGetM(prhs(1))
nx = mxGetM(prhs(2))
nw = mxGetM(prhs(3))

allocate(x(nx),t(nt),w(nw),S(nt,nx))

call mxCopyPtrToReal8(mxGetPr(prhs(1)),t,nt)
call mxCopyPtrToReal8(mxGetPr(prhs(2)),x,nx)
call mxCopyPtrToReal8(mxGetPr(prhs(3)),w,nw)


do i = 1,nt
   flag = 0
   do j = 1,nx
      if (abs(t(i)-x(j)) .lt. 1.0d-16) then
         temp = j
         flag = 1
      endif
   enddo
   if (abs(flag) .lt. 1.0d-12) then
      do j = 1,nx
         S(i,j) = w(j)/(t(i)-x(j))
      enddo
      S(i,:) = S(i,:)/sum(S(i,:))
   else
      S(i,:) = 0
      S(i,temp) = 1
   endif
enddo

plhs(1) = mxCreateDoubleMatrix(nt, nx, 0)
call mxCopyReal8ToPtr(S, mxGetPr(plhs(1)),nt*nx)

deallocate(x,t,w,S)
end subroutine
