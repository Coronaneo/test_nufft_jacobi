function [fun,ts,wghts] = JPT1D(nts,da,db,tR,mR,tol)

    if nts < 2^12
       it = 27;
    else
       it = 9;
    end
nt = zeros(nts,1);
[ts,wghts] = getts(nt,da,db);
nu = [it:nts-1]';
xs = mod(floor(ts*nts/2/pi),nts)+1;
[U,V] = lowrank(nts,JTM1d,da,db,tol,tR,mR,ts,nu);
rank = size(U,2);

fun = @(c)JacPT1d(c);

    function y = JacPT1d(c)
        ncol = size(c,2);
        d = repmat(conj(V),1,ncol).*reshape(repmat(c,rank,1),nts,rank*ncol);
        fftc = ifft(d);
        fftc = fftc(xs,:);
        y = nts*squeeze(sum(reshape(repmat(U,1,ncol).*fftc,nts,rank,ncol),2));
        y = real(y)./wghts;
    end

end