function [fun,rank1,rank2] = JPT2D(nts,da,db,tR,mR,tol,opt,R_or_N)
%  Return:(a function handle computing 2D uniform Jacobi polynomial transform)
%    fun(c) = kron(J,J)*c(:), 
%    where W = diag(wghts);
%      if R_or_N > 0, J(j,k) = M^(da,db)_(ts(j),k-1)*exp(1i*(psi^(da,db)_(t(j),k-1)-2*pi/nts*[ts(j)*nts/2/pi]*(k-1))), 1=<j,k<=nts;               
%      if R_or_N > 0, J(j,k) = M^(da,db)_(ts(j),k-1)*exp(1i*(psi^(da,db)_(t(j),k-1)-t(j)*(k-1))), 1=<j,k<=nts;
%  ts    - built-in uniform samples
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
       it = 9;
    else
       it = 27;
    end
nt = zeros(nts,1);
[ts,wghts] = getts(nt,da,db);
nu = [it:nts-1]';
xs = mod(floor(ts*nts/2/pi),nts)+1;


if opt >= 1
    JTM = @(ts,nu)interpjac1(nt,ts,nu,da,db,R_or_N);
    [U,V] = lowrank(nts,JTM,ts,nu,tol,tR,mR);
    %U = diag(sqrt(wghts))*U;
    V = conj(V);
elseif 0 <= opt && opt<1
    [U,V] = ID_Cheby1(nts,ts,nu,da,db,tol,1,R_or_N,tR,mR);
    %U = diag(sqrt(wghts))*U;
elseif opt < 0
    [U,V] = ID_Cheby1(nts,ts,nu,da,db,tol,-1,R_or_N,tR,mR);
    %U = diag(sqrt(wghts))*U;
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

vals1 = jacrecur(nts,ts,it-1,da,db);
vals2 = jacrecur(nts,ts,it-1,da,db);
vals = kron(vals1,vals2);

xs1 = xs;
xs2 = xs;
xsub = zeros(nts*nts,1);
for i = 1:nts
    for j =1:nts
        xsub((i-1)*nts+j) = (xs1(i)-1)*nts + xs2(j);
    end
end

if  R_or_N > 0
    fun = @(c)JacPT2d1(c);
else
    ex = exp(1i*nts/2*ts);
    U = U.*repmat(ex,1,rank);
    fun = @(c)JacPT2d2(c);
end

    function y = JacPT2d1(c)
        c = c(:);
        c1 = c(reshape(repmat([0:it-1],it,1),it*it,1)*nts+repmat([1:it]',it,1));
        y = vals*c1;

        %z1 = kron([vals1 zeros(nts,nts-it)],[vals2 zeros(nts,nts-it)])*c;
	%e1 = norm(z1-y)/norm(z1)

        c1 = c(1:it*nts);
        d = kron(ones(it,1),V2).*repmat(c1,1,rank2);
        fft2c = zeros(it*nts,rank2);
        for i = 1:it
            d1 = ifft(d((i-1)*nts+1:i*nts,:));
            fft2c((i-1)*nts+1:i*nts,:) = d1(xs,:);
        end
        fft2c = nts*fft2c;
        y2 = zeros(nts*nts,rank2);
        for i = 1:nts
            for j = 1:it
                y2((i-1)*nts+1:i*nts,:) = y2((i-1)*nts+1:i*nts,:) + vals1(i,j)*U2.*fft2c((j-1)*nts+1:j*nts,:);
            end
        end

        %F = exp(2*pi*1i/nts*(xs-1)*[0:nts-1]);
        %FF = zeros(nts,nts);
        %for kk = 1:rank2
        %    FF = FF + diag(U1(:,kk))*F*diag(V1(:,kk));
        %end
	%z2 = kron([vals1 zeros(nts,nts-it)],FF)*c;
	%e2 = norm(z2-sum(y2,2))/norm(z2)


        y = y + sum(y2,2);
        c1 = c(reshape(repmat([0:nts-1],it,1),it*nts,1)*nts+repmat([1:it]',nts,1));
        d = kron(V1,ones(it,1)).*repmat(c1,1,rank1);
        sl = [1:it:it*nts]';
        fft2c = zeros(nts*it,rank1);
        for i = 1:it
            d1 = nts*ifft(d(sl+(i-1),:));
             fft2c(sl+(i-1),:) = d1(xs,:);
        end
        y2 = zeros(nts*nts,rank1);
        for i = 1:nts 
            y2((i-1)*nts+1:i*nts,:) = repmat(U1(i,:),nts,1).*(vals2*fft2c((i-1)*it+1:i*it,:));
        end
        
	zz = zeros(nts*nts,1);
	for kk = 1:rank1
	    zz = zz + kron(diag(U(:,kk)),vals2)*fft2c(:,kk);
	end

     
        %z3 = kron(FF,[vals2 zeros(nts,nts-it)])*c;
        %e3 = norm(z3-sum(y2,2))/norm(z3)


        y = y + sum(y2,2);
        d = VV.*repmat(c,1,rank1*rank2);   
        d = reshape(d,nts,nts,rank1*rank2);
        fft2c = ifft2(d);
        fft2c = nts^2*reshape(fft2c,nts^2,rank1*rank2);
        fft2c = fft2c(xsub,:);
        y2 = UU.*fft2c;


	%z4 = kron(FF,FF)*c;
	%e4 = norm(z4-sum(y2,2))/norm(z4)

        y = y + sum(y2,2);
        y = real(y);
	%z = real(z1+z2+z3+z4);
	%e = norm(z-y)/norm(z)
	%y = z;
    end

    function y = JacPT2d2(c)
        y = zeros(nts,1);
        for i=1:rank
            cj = nufft1dIInyumex(ts,1,tol,V(:,i).*c);
            y = y + U(:,i).*cj;
        end
	y = real(y)./sqrt(wghts);
    end



end
