function [fun,rank] = NJPT1D(nts,ts,da,db,tR,mR,tol,opt)

    if nts < 2^12
       it = 9;
    else
       it = 27;
    end
nu = [it:nts-1]';
xs = mod(floor(ts*nts/2/pi),nts)+1;
wghts = ones(nts,1);

if opt > 0
    JTM = @(rs,cs,n,da,db,ts,nu,wghts)JTM1d(rs,cs,n,da,db,ts,nu,wghts);
   [U,V] = lowrank(nts,JTM,da,db,tol,tR,mR,ts,nu,wghts);
   V = conj(V);
else
    JTM = @(rs,cs,ts,nu)JTM1d(rs,cs,n,da,db,ts,nu,wghts);
    grid = cos(((2*[nts:-1:1]'-1)*pi/2/nts)+1)*pi/2;
    [U,V] = ID_Cheby(JTM,ts,nu,wghts,grid,50,tol,'r',1);
end
rank = size(U,2);

fun = @(c)NJacPT1d(c);

    function y = NJacPT1d(c)
        ncol = size(c,2);
        d = repmat(V,1,ncol).*reshape(repmat(c,rank,1),nts,rank*ncol);
        fftc = ifft(d);
        fftc = fftc(xs,:);
        y = nts*squeeze(sum(reshape(repmat(U,1,ncol).*fftc,nts,rank,ncol),2));
        y = real(y);
    end

end
