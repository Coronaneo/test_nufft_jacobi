function [U,V] = ID_Cheby(n,x,k,da,db,tol,opt,R_or_N,tR,mR)
% Compute decomposition A = U*V.' via ID approximation A(:,rd) ~ A(:,sk)*T. 
% A =fun(rs,cs)
% x -- sample
% k -- degree
%
% tol -- accuracy
% da,db -- parameters of jacobi polynomial
% The precision is specified by tol and the rank is given by 'mR'; 
% 
% opt - whether use adaptive cheb interpolation for both sample dimension and degree dimension
%       if opt = 0, just use cheb interpolation for samples
%       if opt = 1, cheb interpolation for both samples and degrees
%        
% R_or_N -- decides which form of the matrix to be approximated
%       if R_or_N > 0, use the form : M^(da,db)(v,t)exp(i(psi^(a,b)(v,t)-2pi/n*[t*n/2pi]v))
%       if R_or_N <= 0, use the form : M^(da,db)(v,t)exp(i(psi^(a,b)(v,t)-tv))
%
% mR -- maximum rank
% tR -- p*mR, where p should be a oversampling paramter
% Copyright 2018 Haizhao Yang, Qiyuan Pang


%dd      = n;
%dd      = min(0.10,1/dd);
%dd      = log2(dd);
%nints   = ceil(-dd)+1;
%nints   = 2*nints;
%5kktol = ceil((log2(sqrt(pi/2)/tol)-0)/(4-log2(pi*exp(1))))*ceil(log2(n));
%kk = max(20*(ceil(log2(n))+1),kktol);


kk = tR;
chebygrid = cos([kk-1:-1:0]'*pi/(kk-1));

%%%%%%%% chebyshev grids on (0,pi)
ts = (pi/2)*chebygrid+pi/2;

nt = zeros(n,1);
%%%%%%%% decide which cheb method to use 
if  opt > 0
    nu = k;
else
%%%%%%%%%%%%%%  chebyshev grids on (n-26,n) or (n-8,n), depends on n
    xx = 5*(ceil(log2(n))+1);
    chebygrid = cos([xx-1:-1:0]'*pi/(xx-1));
    nu = (n-k(1)-1)/2*chebygrid1+(n+k(1)-1)/2;
end

%%%%%%%%%% consturct lowrank approximation for the matrix obtained via chebyshev grids
    A = interpjac1(nt,ts,nu,da,db,R_or_N);
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

%%%%%%%%% construct right factor U in fun(x,k) = U*V.'
    U = zeros(size(x,1),rr);
    w = [1/2 (-1).^[1:kk-2] 1/2*(-1)^(kk-1)]';
    S = barcycheby(x,ts,w);
    U = S*U1;
    
if  opt > 0
%%%%%%%%%% construct left factor V in fun(x,k) = U*V.'
    V = V1;
else
V = zeros(size(k,1),rr);
w = [1/2 (-1).^[1:xx-2] 1/2*(-1)^(xx-1)]';
P = Lagrange(k,nu,w);
V = P*V1;

end

end
