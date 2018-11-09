function [yfun,fun_pre] = invNJPT1D(nts,ts,da,db,tol)

nt = zeros(nts,1);
nu = [0:nts-1]';
wghts = ones(nts,1);
n1 = [1:nts];
fun = @(x)directjac1(nt,x,da,db,n1,ts,nu,wghts);
rs = [1:nts]';
cs = [1:nts]';
A = extrjac1(nt,rs,cs,-1,da,db,ts,nu,wghts);
[U,S,V] = svd(A,'econ');
u = sqrt(S(1,1))*U(:,1);
v = sqrt(S(1,1))*conj(V(:,1));
N = exp(1i*ts*nu');
M = N'*N;
a = ones(nts,1);
m1 = N'*a;
m2 = N.'*a;
v = zeros(nts,1);
v(1) = m1(1);
for i = 2:nts
    v(i) = ((nts-i+1)*m1(i)+(i-1)*m2(nts-i+2))/nts;
end
X = [[2:nts]';1];
Y = [1:nts]';
W = ones(nts,1);
WS = sparse(X,Y,W,nts,nts);
PRE = zeros(nts,nts);
for i = 1:nts
    PRE(:,i) = (WS^(i-1))*v;
end
fun_pre = @(y)invNJPT1Dpre(nts,ts,nu,da,db,u,v,y,tol,N,M,PRE);
restart = 10;
maxit =  10;
yfun = @(b)invNJacPT1D(b);
    function [x,flag,relres,iter] = invNJacPT1D(b)
             [x,flag,relres,iter] = gmres(fun,b,restart,tol,maxit,fun_pre);
    end
end