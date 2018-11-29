function [fun,rank1,rank2] = invJPT2D(nts,da,db,tR,mR,tol,opt,R_or_N)
%  Return:fun(c)(a function handle computing 1D inverse uniform Jacobi
%  polynomial transform) such that
%    J*fun(c) = c, 
%    where W = diag(wghts);
%      if R_or_N > 0, J(j,k) = M^(da,db)_(ts(j),k-1)*exp(1i*(psi^(da,db)_(t(j),k-1)-2*pi/nts*[ts(j)*nts/2/pi]*(k-1))), 1=<j,k<=nts;               
%      if R_or_N > 0, J(j,k) = M^(da,db)_(ts(j),k-1)*exp(1i*(psi^(da,db)_(t(j),k-1)-t(j)*(k-1))), 1=<j,k<=nts;
%
%  ts    - specific samples for users
%  mR    - maximum rank    
%  tR    - p*mR, where p>5 sould be a oversampling parameter
%  tol   - accuracy
%  opt   - an option decides which lowrank approximation to use
%          opt >= 1, use Randomized sampling SVD
%          0 <= opt < 1, use lowrank approximation based on applying chebyshev grids interpolation for sample dimension                    
%          opt <= 0, use lowrank approximation based on applying chebyshev grids interpolation for both sample and degree dimensions
%         
%  Copyright reserved by Qiyuan Pang, 22/11/2018 

   if nts < 2^12
       it = 10;
    else
       it = 28;
    end
nt = zeros(nts,1);
[ts,wghts] = getts(nt,da,db);
nu = [it:nts-1]';
xs = mod(floor(ts*nts/2/pi),nts)+1;

ts1 = ts;
ts2 = ts;
Y = [1:nts]';
X = zeros(nts,1);
S = ones(nts,1);
for i = 1:nts
    X(i) = xs(i);   
end
P1 = sparse(X,Y,S,nts,nts,nts);
P2 = P1;
PP = kron(P2,P1);

if opt >= 1
    JTM = @(ts,nu)interpjac1(nt,ts,nu,da,db,R_or_N);
    [U,V] = lowrank(nts,JTM,ts,nu,tol,tR,mR);
    U = diag(sqrt(wghts))*U;
    V = conj(V);
elseif 0 <= opt && opt<1
    [U,V] = ID_Cheby(nts,ts,nu,da,db,tol,1,R_or_N,tR,mR);
    U = diag(sqrt(wghts))*U;
elseif opt < 0
    [U,V] = ID_Cheby(nts,ts,nu,da,db,tol,-1,R_or_N,tR,mR);
    U = diag(sqrt(wghts))*U;
end
U1 = U;
U2 = U;
V1 = V;
V2 = V;
rank1 = size(U1,2);
rank2 = size(U2,2);
V1 = [zeros(it,rank1);V1];
V2 = [zeros(it,rank2);V2];
UU = kron(U1,U2);
VV = kron(V1,V2);
rank = rank1*rank2;

vals1 = diag(sqrt(wghts))*jacrecur(nts,ts,it-1,da,db);
vals1 = vals1.';
vals2 = diag(sqrt(wghts))*jacrecur(nts,ts,it-1,da,db);
vals2 = vals2.';
vals = kron(vals1,vals2);

if  R_or_N > 0
    fun = @(c)invJacPT2d1(c);
else
    ex = exp(1i*nts/2*ts);
    U = U.*repmat(ex,1,rank);
    fun = @(c)invJacPT2d2(c);
end

    function y = invJacPT2d1(c)
        c = c(:);
        c = kron(1./sqrt(wghts),1./sqrt(wghts)).*c;
        y = zeros(nts*nts,1);
        ind = reshape(repmat([0:it-1],it,1),it*it,1)*nts+repmat([1:it]',it,1);
        y(ind,:) = vals*c;

        %z1 = kron([vals2;zeros(nts-it,nts)],[vals1;zeros(nts-it,nts)])*c;
	%e1 = norm(z1-y)/norm(z1)

        d = zeros(nts*it,rank1);
        for i = 1:it
            for j = 1:nts
                d((i-1)*nts+1:i*nts,:) = d((i-1)*nts+1:i*nts,:) + vals2(i,j)*U1.*c((j-1)*nts+1:j*nts,:);
            end
        end
        fft2c = zeros(nts*it,rank1);
        for i = 1:it
            d1 = P1*d((i-1)*nts+1:i*nts,:);
            fft2c((i-1)*nts+1:i*nts,:) = conj(fft(conj(d1)));
        end
        y2 = kron(ones(it,1),V1).*fft2c;
        ind = [1:it*nts]';
        y(ind,:) = y(ind,:) + sum(real(y2),2);

        %F = exp(1i*2*pi/nts*(xs-1)*[0:nts-1]).';
	%FF = zeros(nts,nts);
	%for kk = 1:rank1
	%    FF = FF + diag(V1(:,kk))*F*diag(U1(:,kk));
	%end
	%z2 = kron([vals2;zeros(nts-it,nts)],FF)*c;
	%e2 = norm(z2(ind,:)-sum(y2,2))/norm(z2(ind,:))

        d = kron(U2,ones(nts,1)).*repmat(c,1,rank2);
        sl = [1:nts:nts*nts]';
        fft2c = zeros(nts*nts,rank2);
        for i = 1:nts
            d1 = P2*d(sl+(i-1),:);
            fft2c(sl+(i-1),:) = conj(fft(conj(d1)));
        end
        y2 = zeros(nts*it,rank2);
        for i = 1:nts 
            y2((i-1)*it+1:i*it,:) = repmat(V2(i,:),it,1).*(vals1*fft2c((i-1)*nts+1:i*nts,:));
        end
        ind = reshape(repmat([0:nts-1],it,1),it*nts,1)*nts+repmat([1:it]',nts,1);
        y(ind,:) = y(ind,:) + sum(real(y2),2);

        %z3 = kron(FF,[vals1;zeros(nts-it,nts)])*c;
	%e3 = norm(z3(ind,:)-sum(y2,2))/norm(z3(ind,:))

        d = UU.*repmat(c,1,rank1*rank2); 
        d = PP*d;
        d = reshape(d,nts,nts,rank1*rank2);
        fft2c = conj(fft2(conj(d)));
        fft2c = reshape(fft2c,nts^2,rank1*rank2);
        y2 = VV.*fft2c;

        %z4 = kron(FF,FF)*c;
	%e4 = norm(z4-sum(y2,2))/norm(z4)

        y = y + sum(real(y2),2);

    end

end
