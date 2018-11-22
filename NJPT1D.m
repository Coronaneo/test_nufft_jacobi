function [fun,rank] = NJPT1D(nts,ts,da,db,tR,mR,tol,opt,R_or_N)

    if nts < 2^12
       it = 9;
    else
       it = 27;
    end
nu = [it:nts-1]';
xs = mod(floor(ts*nts/2/pi),nts)+1;
wghts = ones(nts,1);
nt = zeros(nts,1);

if  opt > 0
    JTM = @(ts,nu)interpjac1(nt,ts,nu,da,db,R_or_N);
    [U,V] = lowrank(nts,JTM,ts,nu,tol,tR,mR);
     V = conj(V);
elseif 0 <= opt && opt<1
    [U,V] = ID_Cheby1(nts,ts,nu,da,db,tol,1,R_or_N,tR,mR);
elseif opt < 0
    [U,V] = ID_Cheby1(nts,ts,nu,da,db,tol,-1,R_or_N,tR,mR);
end
end
rank = size(U,2);
V = [zeros(it,rank);V];

if  R_or_N > 0
    fun = @(c)NJacPT1d1(c);
else
    ex = exp(1i*nts/2*ts);
    U = U.*repmat(ex,1,rank);
    fun = @(c)NJacPT1d2(c);
end

    function y = NJacPT1d1(c)
        ncol = size(c,2);
        d = repmat(V,1,ncol).*reshape(repmat(c,rank,1),nts,rank*ncol);
        fftc = ifft(d);
        fftc = fftc(xs,:);
        y = nts*squeeze(sum(reshape(repmat(U,1,ncol).*fftc,nts,rank,ncol),2));
        y = real(y);
    end
    
    function y = NJacPT1d2(c)
        y = zeros(nts,1);
        for i=1:rank
            cj = nufft1dIInyumex(ts,1,tol,V(:,i).*c);
            y = y + U(:,i).*cj;
        end
	y = real(y)./sqrt(wghts);
    end
end
