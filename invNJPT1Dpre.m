function [x] = invNJPT1Dpre(nts,ts,nu,da,db,u,v,y,tol,N,M,PRE)
%nt = zeros(nts,1);
%rs = [1:nts]';
%cs = [1:nts]';
%nu = [0:nts-1]';
%wghts = ones(nts,1);
%A = extrjac1(nt,rs,cs,-1,da,db,ts,nu,wghts);
%[U,S,V] = svd(A,'econ');
%u = sqrt(S(1,1))*U(:,1);
%v = sqrt(S(1,1))*conj(V(:,1));
y1 = y./u;
y2 = (N')*y1;
[x1,flag,relres,iter] = pcg(M,y2,tol,100,PRE);
x = x1./v;
end