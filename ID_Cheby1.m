function [U,V] = ID_Cheby1(n,x,k,wghts,tol,opt)
% Compute decomposition A = U*V.' via ID approximation A(:,rd) ~ A(:,sk)*T. 
% A =fun(rs,cs,x,k)
% The precision is specified by tol and the rank is given by 'rank'; 
% r_or_c specify row ID or column ID
% opt - whether use adaptive rank or fix rank
%       if opt = 0, fix rank rk in the ID no matter what tol is
%       if opt = 1, use a rank less than or equal to rk trying to obtain an
%        accuracy as good as tol; may not be able to achieve tol if rk is
%        too small
%
% Copyright 2018 Haizhao Yang, Qiyuan Pang

kk = 16;
dd      = n;
dd      = min(0.10,1/dd);
dd      = log2(dd);
nints   = ceil(-dd)+1;
nints   = 2*nints;
chebygrid = cos((2*[kk:-1:1]'-1)*pi/2/kk);


nints0  = nints/2;
dd      = 2.0;

ab = zeros(2,nints);
for int=1,nints0
ab(1,int) = dd^(-nints0+int-1);
ab(2,int) = dd^(-nints0+int);
end 

ab = pi/2*ab;

for int=1,nints0
ab(1,nints-int+1) = pi - ab(2,int);
ab(2,nints-int+1) = pi - ab(1,int);
end 

ts = zeros(kk*nints,1);
for i = 1:nints
    ts((i-1)*kk+1:i*kk) = (chebygrid+1)/2*(ab(2,i)-ab(1,i))+ab(1,i);
end

nt = zeros(n,1);
if opt > 0
   nu = k;
else
end
A = interpjac1(nt,ts,nu,1);
[~,R,E] = qr(A',0);
rr = find( abs(diag(R)/R(1)) > tol, 1, 'last');
sk = E(1:rr);
rd = E(rr+1:end);
T = R(1:rr,1:rr)\R(1:rr,rr+1:end);
T = T';
V1 = A(sk,:);
V1 = V1.';
U1 = zeros(kk*nints,rr);
    for i = 1:kk*nints
        flag = find(sk == i);
        if ~isempty(flag)
            U1(i,flag) = 1;
        else
            flag1 = find(rd == i);
            U1(i,:) = T(flag1,:);
        end
    end

binranges = [ab(1,:) ab(2,end)]';
bincounts = histc(x,binranges);

U = zeros(size(x,1),rr);
totalM = 0;
for i = 1:nints
    w = ts((i-1)*kk+1:i*kk);
    ll = zeros(kk,1);
    for j = 1:kk
        w1 = [w(1:j-1) w(j+1:end)];
        ll(j) = prod(ones(1,kk-1)*w(j)-w1);
    end
    S = zeros(bincounts(i),kk);
    for j = 1:bincounts(i)
        omega = ones(1,kk)*x(totalM+j)-w;
        flag = find(omega == 0);
        if  isempty(flag)
            omega1 = prod(omega);
            ww = (omega1*ones(1,kk)./omega)./ll;
        else
            ww = zeros(1,kk);
            ww(1,flag) = 1;
        end
        S(j,:) = ww;
    end
    U(totalM+1:totalM+bincounts(i),:) = S*U1((i-1)*kk+1:i*kk,:);
    totalM = totalM + bincounts(i);
end
sqrtW = diag(sqrt(wghts));
U = sqrtW*U;

if  opt > 0
    V = V1;
else
end
end