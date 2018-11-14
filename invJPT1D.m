function [fun,rank,ts,wghts] = invJPT1D(nts,da,db,tR,mR,tol)

    if nts < 2^12
       it = 9;
    else
       it = 27;
    end
nt = zeros(nts,1);
[ts,wghts] = getts(nt,da,db);
nu = [it:nts-1]';
xs = mod(floor(ts*nts/2/pi),nts)+1;
cs = zeros(nts,1);
Y = [1:nts]';
X = zeros(nts,1);
S = ones(nts,1);
for i = 1:nts
    X(i) = xs(i);   
end
P = sparse(X,Y,S,nts,nts,nts);
JTM = @(rs,cs,n,da,db,ts,nu,wghts)JTM1d(rs,cs,n,da,db,ts,nu,wghts);
[U,V] = lowrank(nts,JTM,da,db,tol,tR,mR,ts,nu,wghts);
rank = size(U,2);

fun = @(c)invJacPT1d(c);

    function y = invJacPT1d(c)
        ncol = size(c,2);
        c = c.*repmat(sqrt(wghts),1,ncol);
        d = repmat(U,1,ncol).*reshape(repmat(c,rank,1),nts,rank*ncol);
        d = P*d;
        fftc = conj(fft(conj(d)));
        y = squeeze(sum(reshape(repmat(conj(V),1,ncol).*fftc,nts,rank,ncol),2));
        y = real(y);
    end

end