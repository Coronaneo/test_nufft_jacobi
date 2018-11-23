function y = directjac2d(nts,s,t,n,da,db,c)
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

for i = 1:m
    jj = floor((n(i)-0.5)/nts);
    ii = n(i) - nts*jj;
    x = interpjac1(nt,s(jj+1),nu,da,db,-1);
    %x1 = exp(1i*2*pi/nts*floor(s(jj+1)*nts/2/pi)*nu');
    x1 = exp(1i*s(jj+1)*nu');
    x = x.*x1;
    z = interpjac1(nt,t(ii),nu,da,db,-1);
    %z1 = exp(1i*2*pi/nts*floor(t(ii)*nts/2/pi)*nu');
    z1 = exp(1i*t(ii)*nu');
    z = z.*z1;
    x = [vals(jj+1,:) x];
    z = [valt(ii,:) z];
    y(i) = real(kron(x,z)*c);
end
end
