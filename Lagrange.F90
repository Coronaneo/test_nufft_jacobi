
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
mwSize       :: n,kk,


real*8, allocatable :: ts(:),x(:),w(:),ll(:,:),S(:,:),omega(:)
integer*4 ii,i,jj,nint
real*8 da,db,flag,pi,dd

kk = mxGetM(prhs(1))
xx = mxGetM(prhs(2))

allocate(ts(kk),x(xx),ll(kk-1,kk),w(kk-1))
call mxCopyPtrToReal8(mxGetPr(prhs(1)),ts,kk)
call mxCopyPtrToReal8(mxGetPr(prhs(2)),x,xx)

do j = 1,kk
   if (j == 1) then
       w = ts(2:kk)
   end if
   if (j == kk) then
       w = ts(1:kk-1)
   end if
   if (1<j & j<kk) then
       w(1:j-1) = ts(1:j-1)
       w(j:kk-1) = ts(j+1:kk)
   end if
   ll(:,j) = ts(j)-w;
end do
    
allocate(S(xx,kk),omega(kk))
do j = 1,xx
   omega = x(j)-ts;
   flag = find(abs(omega) <= eps);
   if  isempty(flag)
            
            omega1 = prod(omega);
            ww = (omega1*ones(kk,1)./omega);
	        for jj = 1:kk
	     	    for ii = 1:kk-1 
	             	ww(jj) = ww(jj)*ll(ii,jj);
                    end
	        end
        else
            ww = zeros(kk,1);
            ww(flag,1) = 1;
        end
        S(j,:) = ww.';
    end

end subroutine