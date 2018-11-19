function [U,V] = ID_Cheby1(n,x,k,wghts,da,db,tol,opt,R_or_N,tR,mR)
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


dd      = n;
dd      = min(0.10,1/dd);
dd      = log2(dd);
nints   = ceil(-dd)+1;
nints   = 2*nints;
%5kktol = ceil((log2(sqrt(pi/2)/tol)-0)/(4-log2(pi*exp(1))))*ceil(log2(n));
%kk = max(20*(ceil(log2(n))+1),kktol);
kk = tR;
chebygrid = cos((2*[kk:-1:1]'-1)*pi/2/kk);
%chebygrid = cos([kk-1:-1:0]'/(kk-1)*pi);
%chebygrid = cos((kk-[1:kk]')*pi/(kk-1))
%chebygrid = cos((2*[kk:-1:1]'-1)*pi/2/kk);

nints0  = nints/2;
dd      = 2.0;

ab = zeros(2,nints);
for int=1:nints0
ab(1,int) = dd^(-nints0+int-1);
ab(2,int) = dd^(-nints0+int);
end 

ab = pi/2*ab;

for int=1:nints0
ab(1,nints-int+1) = pi - ab(2,int);
ab(2,nints-int+1) = pi - ab(1,int);
end 

%ts = zeros(kk*nints,1);
%for i = 1:nints
%    ts((i-1)*kk+1:i*kk) = (chebygrid)/2*(ab(2,i)-ab(1,i))+(ab(1,i)+ab(2,i))/2;
%end
ts = (pi/2-1/n)*chebygrid+pi/2;

nt = zeros(n,1);
if  opt > 0
    nu = k;
else
    xx = 5*(ceil(log2(n))+1);
    chebygrid1 = cos((2*[xx:-1:1]'-1)*pi/2/xx);
    nu = (n-k(1)+1)/2*chebygrid1+(n+k(1)-1)/2;
end
[A,ier] = interpjac1(nt,ts,nu,da,db,R_or_N);
%ier
[~,R,E] = qr(A',0);
rr = find( abs(diag(R)/R(1)) > tol, 1, 'last');
rr = min(rr,mR);
sk = E(1:rr);
rd = E(rr+1:end);
T = R(1:rr,1:rr)\R(1:rr,rr+1:end);
T = T';
V1 = A(sk,:);
V1 = V1.';
nn = kk;
U1 = zeros(nn,rr);
    for i = 1:nn
        flag = find(sk == i);
        if ~isempty(flag)
            U1(i,flag) = 1;
        else
            flag1 = find(rd == i);
            U1(i,:) = T(flag1,:);
        end
    end
%norm(A-U1*V1.')
binranges = [ab(1,:) ab(2,end)]';
bincounts = histc(x,binranges);
%bincounts = histc(x,[ts(1) ts(end)])

U = zeros(size(x,1),rr);
nint = 1;
%SS=zeros(size(x,1),kk*nint);
totalM = 0;
totalN = 0;
for i = 1:nint
    w = ts((i-1)*kk+1:i*kk);
    ll = zeros(kk-1,kk);
    for j = 1:kk
        if j == 1
            w1 = w(2:end);
	end
    	if j == kk
            w1 = w(1:end-1);
        end
	if 1<j && j<kk
            w1 = [w(1:j-1); w(j+1:end)];
        end
        ll(:,j) = 1./(ones(kk-1,1)*w(j)-w1);
    end
    
    count = size(x,1);%bincounts(i);
    S = zeros(count,kk);
    for j = 1:count
        omega = ones(kk,1)*x(totalM+j)-w;
        flag = find(abs(omega) <= eps);
        if  isempty(flag)
            
            omega1 = prod(omega);
            ww = (omega1*ones(kk,1)./omega);
	        for jj = 1:kk
	     	    for ii = 1:kk-1 
	             	ww(jj) = ww(jj)*ll(ii,jj);
                    end
	        end
        else
            ww = zeros(kk,1);
            ww(flag,1) = 1;
        end
        S(j,:) = ww.';
    end
    %SS(totalM+1:totalM+count,totalN+1:totalN+kk) = S;
    U(totalM+1:totalM+count,:) = S*U1((i-1)*kk+1:i*kk,:);
    %totalM = totalM + count;
    %totalN = totalN + kk;
end

%norm(SS,2)
sqrtW = diag(sqrt(wghts));
U = sqrtW*U;



if  opt > 0
    V = V1;
else
V = zeros(size(k,1),rr);
nint = 1;
%PP=zeros(size(k,1),xx*nint);
totalM = 0;
totalN = 0;
for i = 1:nint
    w = nu((i-1)*xx+1:i*xx);
    ll = zeros(xx-1,xx);
    for j = 1:xx
        if j == 1
            w1 = w(2:end);
        end
        if j == xx
            w1 = w(1:end-1);
        end
        if 1<j && j<xx
            w1 = [w(1:j-1); w(j+1:end)];
        end
        ll(:,j) = 1./(ones(xx-1,1)*w(j)-w1);
    end
    
    count = size(k,1);%bincounts(i);
    P = zeros(count,xx);
    for j = 1:count
        omega = ones(xx,1)*k(totalM+j)-w;
	
        flag = find(abs(omega) <= eps);
        if  isempty(flag)

            omega1 = prod(omega);
            ww = omega1*ones(xx,1)./omega;
	    
                for jj = 1:xx
                    for ii = 1:xx-1
                        ww(jj) = ww(jj)*ll(ii,jj);
                    end
                end
        else
            ww = zeros(xx,1);
            ww(flag,1) = 1;
        end
        P(j,:) = ww.';
    end
    %PP(totalM+1:totalM+count,totalN+1:totalN+xx) = P;
    V(totalM+1:totalM+count,:) = P*V1((i-1)*xx+1:i*xx,:);
    %totalM = totalM + count;
    %totalN = totalN + kk;
end

end
%[B,ier] = interpjac1(nt,x,k,da,db,R_or_N);
%max(max(abs(SS)))
%max(max(abs(A)))
%max(max(abs(P)))
%sum(sum(isnan(SS*A*PP.')))
%size(SS)
%size(A)
%size(PP)
%norm(B-SS*A*(PP.'))/norm(B)
%norm(sqrtW*B-U*V.')/norm(sqrtW*B)

%C = A.'\B.';
%C(:,1:bincounts(1)).'
end
