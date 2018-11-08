function [fun,rank,ts,wghts] = invJPT1D(nts,da,db,tR,mR,tol)

    if nts < 2^12
       it = 27;
    else
       it = 9;
    end
nt = zeros(nts,1);
[ts,wghts] = getts(nt,da,db);
nu = [it:nts-1]';
xs = mod(floor(ts*nts/2/pi),nts)+1;
JTM = @(rs,cs,n,da,db,ts,nu,wghts)JTM1d(rs,cs,n,da,db,ts,nu,wghts);
[U,V] = lowrank(nts,JTM,da,db,tol,tR,mR,ts,nu,wghts);
rank = size(U,2);

fun = @(c)JacPT1d(c);

    function y = JacPT1d(c)
        ncol = size(c,2);
        c = c./repmat(sqrt(wghts),1,ncol);
        d = repmat(conj(U),1,ncol).*reshape(repmat(c,rank,1),nts,rank*ncol);
        d = d(xs,:);
        fftc = fft(d);
        y = nts*squeeze(sum(reshape(repmat(V,1,ncol).*fftc,nts,rank,ncol),2));
        y = real(y);
    end

end