function y = directjac3d(nts,r,s,t,n1,da,db,c)
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

for i = 1:m
    kk = floor((n(i)-0.5)/nts/nts);
    jj = floor((n(i)-kk*nts*nts-0.5)/nts);
    ii = n(i)-kk*nts*nts-nts*jj;
    x = interpjac1(nt,r(kk+1),nu,da,db,-1);
    %x1 = exp(1i*2*pi/nts*floor(s(jj+1)*nts/2/pi)*nu');
    x1 = exp(1i*r(jj+1)*nu');
    x = x.*x1;
    y = interpjac1(nt,s(jj+1),nu,da,db,-1);
    %z1 = exp(1i*2*pi/nts*floor(t(ii)*nts/2/pi)*nu');
    y1 = exp(1i*s(jj+1)*nu');
    y = y.*y1;
    z = interpjac1(nt,t(ii),nu,da,db,-1);
    %z1 = exp(1i*2*pi/nts*floor(t(ii)*nts/2/pi)*nu');
    z1 = exp(1i*t(ii)*nu');
    z = z.*z1;
    x = [valr(kk+1,:) x];
    y = [vals(jj+1,:) y];
    z = [valt(ii,:) z];
    w(i) = real(kron(kron(x,y),z)*c);
end
end