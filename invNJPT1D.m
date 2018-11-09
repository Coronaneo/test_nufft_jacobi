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
fun_pre = @(y)invNJPT1Dpre(nts,ts,nu,da,db,u,v,y,tol,N,M);
restart = 5;
maxit =  5;
yfun = @(b)gmres(fun,b,restart,tol,maxit,fun_pre);
end