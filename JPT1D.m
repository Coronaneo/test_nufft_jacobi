function [fun,rank,ts,wghts] = JPT1D(nts,da,db,tR,mR,tol)

    if nts < 2^12
       it = 9;
    else
       it = 27;
    end
nt = zeros(nts,1);
[ts,wghts] = getts(nt,da,db);
%size(ts)
%vals0 = jacobi_recurrence(ts,da,db,it);
%size(vals0)
nu = [it:nts-1]';
xs = mod(floor(ts*nts/2/pi),nts)+1;

JTM = @(rs,cs,n,da,db,ts,nu,wghts)JTM1d(rs,cs,n,da,db,ts,nu,wghts);

[U,V] = lowrank(nts,JTM,da,db,tol,tR,mR,ts,nu,wghts);
rank = size(U,2);

fun = @(c)JacPT1d(c);

    function y = JacPT1d(c)
        ncol = size(c,2);
        d = repmat(conj(V),1,ncol).*reshape(repmat(c,rank,1),nts,rank*ncol);
        fftc = ifft(d);
        fftc = fftc(xs,:);
        y = nts*squeeze(sum(reshape(repmat(U,1,ncol).*fftc,nts,rank,ncol),2));
        y = real(y)./sqrt(wghts);
 %       y = y + vals0*c(1:it,:);
    end

end
