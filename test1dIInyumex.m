format long
nj=10;
ms=10;
iflag=1
eps=1e-12
xj=(pi*[0:nj-1]/nj)'
fk=[1:ms]'
ex=exp(1i*nj/2*xj);
[cj] = nufft1dIInyumex(xj,iflag,eps,fk).*ex
nufftc = nufftII(nj*xj/2/pi,iflag,ms,60,eps);
fftc = nufftc(fk)*nj
%fftc = fftshift(fftc)
