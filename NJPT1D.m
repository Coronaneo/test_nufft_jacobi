function [fun,rank] = NJPT1D(nts,ts,da,db,tR,mR,tol)

    if nts < 2^12
       it = 27;
    else
       it = 9;
    end
nu = [it:nts-1]';
xs = mod(floor(ts*nts/2/pi),nts)+1;
JTM = @(rs,cs,n,da,db,ts,nu)JTM1D(rs,cs,n,da,db,ts,nu);
[U,V] = lowrank(nts,JTM,da,db,tol,tR,mR,ts,nu);
rank = size(U,2);

fun = @(c)NJacPT1d(c);

    function y = NJacPT1d(c)
        ncol = size(c,2);
        d = repmat(conj(V),1,ncol).*reshape(repmat(c,rank,1),nts,rank*ncol);
        fftc = ifft(d);
        fftc = fftc(xs,:);
        y = nts*squeeze(sum(reshape(repmat(U,1,ncol).*fftc,nts,rank,ncol),2));
        y = real(y);
    end

end