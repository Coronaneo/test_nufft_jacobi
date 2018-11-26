function w = directinvjac3d(nts,r,s,t,wghtr,wghts,wghtt,n,da,db,c)
if nts < 2^12
   it = 9;
else
   it = 27;
end
valr = jacrecur(nts,r,it-1,da,db);
vals = jacrecur(nts,s,it-1,da,db);
valt = jacrecur(nts,t,it-1,da,db);
m = length(n);
w = zeros(m,1);
nu = [it:nts-1]';
nt = zeros(nts,1);
c = c(:);
c = kron(kron(1./wghtr,1./wghts),1./wghtt).*c;
for i = 1:m
    kk = floor((n(i)-0.5)/nts/nts);
    jj = floor((n(i)-kk*nts*nts-0.5)/nts);
    ii = n(i)-kk*nts*nts-nts*jj;
    if  kk < it
        x = valr(:,kk+1);
        x = sqrt(wghtr).*x;
        x = x.';
    else
        x = interpjac1(nt,r,kk,da,db,-1);
        x = sqrt(wghtr).*x.*exp(1i*r*kk);
        x = x.';
    end
    if  jj < it
        y = valr(:,jj+1);
        y = sqrt(wghts).*y;
        y = y.';
    else
        y = interpjac1(nt,s,jj,da,db,-1);
        y = sqrt(wghts).*y.*exp(1i*s*jj);
        y = y.';
    end
    if  ii < it+1
        z = valr(:,ii);
        z = sqrt(wghtt).*z;
        z = z.';
    else
        z = interpjac1(nt,t,ii-1,da,db,-1);
        z = sqrt(wghtt).*z.*exp(1i*t*(ii-1));
        z = z.';
    end
    w(i) = real(kron(kron(x,y),z)*c);
end
end