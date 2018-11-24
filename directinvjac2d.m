function y = directinvjac2d(nts,t,s,wghtt,wghts,n,da,db,c)
if nts < 2^12
   it = 9;
else
   it = 27;
end
vals = jacrecur(nts,s,it-1,da,db);
valt = jacrecur(nts,t,it-1,da,db);
m = length(n);
y = zeros(m,1);
nu = [it:nts-1]';
nt = zeros(nts,1);
c = c(:);
c = kron(1./sqrt(wghtt),1./sqrt(wghts)).*c;

for i = 1:m
    jj = floor((n(i)-0.5)/nts);
    ii = n(i) - nts*jj;
    if  jj < it
        x = valt(:,jj+1);
        x = sqrt(wghtt).*x;
        x = x.';
    else
        x = interpjac1(nt,t,jj,da,db,-1);
        x = sqrt(wghtt).*x.*exp(1i*t*jj);
        x = x.';
    end
    if  ii < it+1
        z = vals(:,ii);
        z = sqrt(wghts).*z;
        z = z.';
    else
        z = interpjac1(nt,s,ii-1,da,db,-1);
        z = sqrt(wghts).*z.*exp(1i*s*(ii-1));
        z = z.';
    end
    y(i) = real(kron(x,z)*c);
end
end
