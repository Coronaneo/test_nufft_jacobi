#include "fintrf.h"

! standard header for mex functions
subroutine mexfunction(nlhs, plhs, nrhs, prhs)

use utils
use chebyshev
use idecomp
use jacobi_exp
use jacobi_transform

  implicit double precision (a-h,o-z)

  ! number of input arguments, number of output arguments
  integer      :: nlhs, nrhs  
  ! pointer to inputs and outputs
  mwPointer    :: plhs(*), prhs(*) 
  ! get some of the matlab mex functions
  mwPointer    :: mxGetPr, mxGetPi, mxCreateDoubleMatrix 
  ! define a size integer so that we can get its type
  mwSize       :: dmax


type(jacobi_expansion_data)   :: expdata
type(jacobi_transform_data)   :: jacdata

real*8 :: da,db,eps


dmax = mxGetM(prhs(1))
call mxCopyPtrToReal8(mxGetPr(prhs(2)),da,1)
call mxCopyPtrToReal8(mxGetPr(prhs(3)),db,1)
call mxCopyPtrToReal8(mxGetPr(prhs(4)),eps,1)

call jacobi_expansion(eps,1,dmax,da,db,expdata)
call jacobi_transform_prepare(expdata,dmax,jacdata)

plhs(2) = mxCreateDoubleMatrix(dmax, jacdata%krank, 0)
plhs(3) = mxCreateDoubleMatrix(dmax, jacdata%krank, 1)

call mxCopyInteger2ToPtr(jacdata%krank,mxGetPr(plhs(1)),1)
call mxCopyReal8ToPtr(jacdata%r,mxGetPr(plhs(2)),dmax*jacdata%krank)
call mxCopyComplex16ToPtr(jacdata%expvals,mxGetPr(plhs(3)),mxGetPi(plhs(3)),dmax*jacdata%krank)

end subroutine mexfunction
