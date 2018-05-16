function nufftfun = nufftII(x,iflag,ns,rt,tol)
% NUFFTII: Computation of nonuniform FFT in R^1 - Type II.
%   gfun = NUFFTII(x,iflag,ns,rt,tol) provides the fast algorithm for nufft
%   of type II. x is the location of sources on interval [0,ns-1], iflag
%   determines the sign of FFT, and ns is the number of Fourier modes
%   computed. The function returns a function that fast evaluate nufft.
%
%   NUFFT of type II is defined as follows.
%                 ns
%     f(x(j)) = w sum g(k) exp(+/-i 2pi * k x(j) / n) 
%                 k=1
%   for 1 <= k <= ns and 1 <= x(j) <= ns.
%   If (iflag .ge.0) the + sign is used in the exponential and a = 1/ns.
%   If (iflag .lt.0) the - sign is used in the exponential and a = 1.
%
%   See also FFT, NUFFTI, NUFFTIII.


if nargin < 2
    iflag = 1;
else
    iflag = sign(iflag);
end

if nargin < 3
    ns = ceil(max(x))+1;
end

if nargin < 4
    rt = 15;
end

if nargin < 5
    tol = 1e-12;
end

n = size(x,1);
k = (0:ns-1)';

fftconst = iflag*1i/ns*2*pi;

ratiofun = @(x,k)exp(fftconst*(x-round(x))*k');
[U,V] = lowrank(x,k,ratiofun,tol,rt,rt);

xsub = mod(round(x),ns)+1;
r = size(V,2);

nufftfun = @(c)nufftIIfun(c);

    function fftc = nufftIIfun(c)
        ncol = size(c,2);

        c = repmat(conj(V),1,ncol).*reshape(repmat(c,r,1),ns,r*ncol);
        if iflag < 0
            fftc = fft(c);
        else
            fftc = ifft(c);
        end
        fftc = fftc(xsub,:);
        fftc = squeeze(sum(reshape(repmat(U,1,ncol).*fftc,n,r,ncol),2));
    end

end