function [yfun,fun_pre] = invNJPT1D(nts,ts,da,db,tol)

nt = zeros(nts,1);
nu = [0:nts-1];
wghts = ones(nts,1);
n1 = [1:nts];
fun = @(x)directjac1(nt,x,da,db,n1,ts,nu,wghts);
fun_pre = @(y)invNJPT1Dpre(nts,ts,da,db,y,tol);
restart = 5;
maxit =  5;
yfun = @(b)gmres(fun,b,restart,tol,maxit,fun_pre);
end