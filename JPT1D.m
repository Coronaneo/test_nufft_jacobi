function [fun,rank,ts,wghts] = JPT1D(nts,da,db,tR,mR,tol,opt,R_or_N)

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


if opt >= 1
    JTM = @(rs,cs,n,da,db,ts,nu,wghts)JTM1d(rs,cs,n,da,db,ts,nu,wghts,R_or_N);
   [U,V] = lowrank(nts,JTM,da,db,tol,tR,mR,ts,nu,wghts);
   V = conj(V);
elseif 0 <= opt && opt<1
    %JTM = @(rs,cs,ts,nu)JTM1d(rs,cs,nts,da,db,ts,nu,wghts);
    %grid = cos(((2*[nts:-1:1]'-1)*pi/2/nts)+1)*pi/2;
    [U,V] = ID_Cheby1(nts,ts,nu,wghts,da,db,tol,1,R_or_N);
elseif opt < 0
    [U,V] = ID_Cheby1(nts,ts,nu,wghts,da,db,tol,-1,R_or_N);
end
rank = size(U,2);
V = [zeros(it,rank);V];

if R_or_N > 0

   fun = @(c)JacPT1d(c);

    function y = JacPT1d(c)
        ncol = size(c,2);
        d = repmat(V,1,ncol).*reshape(repmat(c,rank,1),nts,rank*ncol);
        fftc = ifft(d);
        fftc = fftc(xs,:);
        y = nts*squeeze(sum(reshape(repmat(U,1,ncol).*fftc,nts,rank,ncol),2));
        y = real(y)./sqrt(wghts);
 %       y = y + vals0*c(1:it,:);
    end
else
    
end
