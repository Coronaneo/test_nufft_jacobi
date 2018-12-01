
#include "fintrf.h"

! standard header for mex functions
subroutine mexfunction(nlhs, plhs, nrhs, prhs)
use utils
use chebyshev
use jacobi_exp
implicit double precision (a-k,o-z)

! number of input arguments, number of output arguments
integer      :: nlhs, nrhs  
! pointer to inputs and outputs
mwPointer    :: plhs(*), prhs(*) 
! get some of the matlab mex functions
mwPointer    :: mxGetPr, mxGetPi, mxCreateDoubleMatrix 
! define a size integer so that we can get its type
mwSize       :: nv,mv,nr,mr

type(jacobi_expansion_data)  :: expdata
real*8, allocatable :: cosvals(:,:),sinvals(:,:),r(:,:)
real*8 eps,iffactor,dmax,da,db
integer*4 iffactor1,i,j

call mxCopyPtrToReal8(mxGetPr(prhs(1)),eps,1)
call mxCopyPtrToReal8(mxGetPr(prhs(2)),iffactor,1)
call mxCopyPtrToReal8(mxGetPr(prhs(3)),dmax,1)
call mxCopyPtrToReal8(mxGetPr(prhs(4)),da,1)
call mxCopyPtrToReal8(mxGetPr(prhs(5)),db,1)
iffactor1 = int(iffactor)
call jacobi_expansion(eps,iffactor1,dmax,da,db,expdata)


if (iffactor .eq. 1.0d0) then
   nv = size(expdata%cosvals0,1)
   mv = size(expdata%cosvals0,2)
   nr = size(expdata%r,1)
   mr = size(expdata%r,2)
else
   nv = size(expdata%cosvals,1)
   mv = size(expdata%cosvals,2)
   nr = 1
   mr = 1
end if
allocate(cosvals(nv,mv),sinvals(nv,mv),r(nr,mr))

!cosvals = expdata%cosvals0
if (iffactor .eq. 1.0d0) then
   do i = 1,nv
      do j = 1,mv
         cosvals(i,j) = expdata%cosvals0(i,j)
         sinvals(i,j) = expdata%sinvals0(i,j)
      end do
   end do
   do i = 1,nr
      do j = 1,mr
         r(i,j) = expdata%r(i,j)
      end do
   end do
else
   do i = 1,nv
      do j = 1,mv
         cosvals(i,j) = expdata%cosvals(i,j)
         sinvals(i,j) = expdata%sinvals(i,j)
      end do
   end do
   r = 1
end if

plhs(1) = mxCreateDoubleMatrix(nv, mv, 0)
plhs(2) = mxCreateDoubleMatrix(nv, mv, 0)
plhs(3) = mxCreateDoubleMatrix(nr, mr, 0)
call mxCopyReal8ToPtr(cosvals, mxGetPr(plhs(1)),nv*mv)
call mxCopyReal8ToPtr(sinvals, mxGetPr(plhs(2)),nv*mv)
call mxCopyReal8ToPtr(r, mxGetPr(plhs(3)),nr*mr)
!plhs(1) = mxCreateDoubleMatrix(1, 1, 0)
!call mxCopyReal8ToPtr(sinvals(nv,mv), mxGetPr(plhs(1)),1)

deallocate(cosvals,sinvals,r)

end subroutine
