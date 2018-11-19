
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
mwSize       :: n,kk,xx


real*8, allocatable :: ts(:),x(:),w(:),ll(:,:),S(:,:),omega(:),c(:),cc(:)
integer*4 ii,i,j,flag
real*8 ome

kk = mxGetM(prhs(1))
xx = mxGetM(prhs(2))

allocate(ts(kk),x(xx),ll(kk-1,kk),w(kk-1),c(kk-1),cc(kk))
call mxCopyPtrToReal8(mxGetPr(prhs(1)),ts,kk)
call mxCopyPtrToReal8(mxGetPr(prhs(2)),x,xx)

do j = 1,kk
   if (j == 1) then
       w = ts(2:kk)
   end if
   if (j == kk) then
       w = ts(1:kk-1)
   end if
   if (1<j .AND. j<kk) then
       w(1:j-1) = ts(1:j-1)
       w(j:kk-1) = ts(j+1:kk)
   end if
   ll(:,j) = ts(j)*c-w;
end do
    
c = 1
cc = 1
allocate(S(xx,kk),omega(kk))
do j = 1,xx
   omega = x(j)*cc-ts;
   ome = 1
   do i = 1,kk
      if (abs(x(j)-ts(i))<1d-15) then
         ome = 0
         flag = i
      else
         ome = ome*(x(j)-ts(i))
      end if
   end do     
   if (abs(ome)>1d-15) then
       do i = 1,kk
          S(j,i) = ome/omega(i)
          do ii = 1,kk-1
             S(j,i) = S(j,i)/ll(ii,i)
          end do
       end do
   else
       S(j,:) = 0.0d0
       S(j,flag) = 1.0d0
   end if
end do
       
plhs(1) = mxCreateDoubleMatrix(xx, kk, 0)
call mxCopyReal8ToPtr(S, mxGetPr(plhs(1)),xx*kk)

deallocate(ts,x,ll,w,S,omega,c,cc)

end subroutine
