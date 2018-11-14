function [U,V] = ID_Cheby(fun,x,k,grid,rank,tol,r_or_c,opt)
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

if nargin < 8, opt = 1; end
if nargin < 7, r_or_c = 'c'; end % sk is column index

rr = length(grid);
xlen = size(x,1);
klen = size(k,1);
switch r_or_c
    case 'c'
        if xlen*klen == 0
            sk = [];
            idx = [];
            rd = 1:klen;
            T = zeros(0,klen);
            return
        end
        
        if rr < xlen
            idxx = round(grid*(xlen-min(xlen,rr)) + (0:min(xlen,rr)-1)')+1;
            idxx= unique(idxx);
        else
            idxx = 1:xlen;
        end
        px = idxx(1:min(xlen,rr));
        pk = [1:klen]';
        Asub = fun(px,pk,x,k);
        
        [~,R,E] = qr(Asub,0);
        if opt > 0
            if xlen*klen > 0
                rr = find( abs(diag(R)/R(1)) > tol, 1, 'last');
                rr = min(rank,rr);
            end
        else
            rr = rank;
        end
        idx = E(1:min(klen,rr));
        rd = E(min(klen,rr)+1:end);
        sk = k(idx,:);
        T = R(1:rr,1:rr)\R(1:rr,rr+1:end);
    case 'r'
        if xlen*klen == 0
            sk = [];
            idx = [];
            rd = 1:xlen;
            T = zeros(xlen,0);
            return
        end
        if rr < klen
            idxk = round(grid*(klen-min(klen,rr)) + (0:min(klen,rr)-1)')+1;
            idxk= unique(idxk);
        else
            idxk = 1:klen;
        end
        pk = idxk(1:min(klen,rr));
        px = [1:xlen]';
        Asub = fun(px,pk,x,k);
        
        [~,R,E] = qr(Asub',0);
        if opt > 0
            if xlen*klen > 0
                rr = find( abs(diag(R)/R(1)) > tol, 1, 'last');
                rr = min(rank,rr);
            end
        else
            rr = rank;
        end
        idx = E(1:min(xlen,rr));
        rd = E(min(xlen,rr)+1:end);
        sk = x(idx,:);
        T = R(1:rr,1:rr)\R(1:rr,rr+1:end);
end
if  r_or_c == 'c'
    U = fun([1:xlen]',idx,x,k);
    V = zeros(klen,min(xlen,rr));
    for i = 1:klen
        flag = find(idx == i);
        if ~isempty(flag)
            V(i,flag) = 1;
        else
            flag1 = find(rd == i);
            V(i,:) = T(:,flag1).';
        end
    end
else
    T = T';
    V = fun(idx,[1:klen]',x,k);
    V = V.';
    U = zeros(xlen,min(klen,rr));
    for i = 1:xlen
        flag = find(idx == i);
        if ~isempty(flag)
            U(i,flag) = 1;
        else
            flag1 = find(rd == i);
            U(i,:) = T(flag1,:);
        end
    end
end
end