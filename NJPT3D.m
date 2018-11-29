function [fun,rank1,rank2,rank3] = NJPT3D(nts,ts1,ts2,ts3,da,db,tR,mR,tol,opt,R_or_N)
%  Return:(a function handle computing 2D nonuniform Jacobi polynomial transform)
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
       it = 10;
    else
       it = 28;
    end
nt = zeros(nts,1);
%[ts,wghts] = getts(nt,da,db);
nu = [it:nts-1]';
%xs = mod(floor(ts*nts/2/pi),nts)+1;


if opt >= 1
    JTM1 = @(ts,nu)interpjac1(nt,ts,nu,da,db,R_or_N);
    [U1,V1] = lowrank(nts,JTM1,ts1,nu,tol,tR,mR);
    %U = diag(sqrt(wghts))*U;
    V1 = conj(V1);
    JTM2 = @(ts,nu)interpjac1(nt,ts,nu,da,db,R_or_N);
    [U2,V2] = lowrank(nts,JTM2,ts2,nu,tol,tR,mR);
    %U = diag(sqrt(wghts))*U;
    V2 = conj(V2);
    JTM3 = @(ts,nu)interpjac1(nt,ts,nu,da,db,R_or_N);
    [U3,V3] = lowrank(nts,JTM1,ts3,nu,tol,tR,mR);
    %U = diag(sqrt(wghts))*U;
    V3 = conj(V3);
elseif 0 <= opt && opt<1
    [U1,V1] = ID_Cheby(nts,ts1,nu,da,db,tol,1,R_or_N,tR,mR);
    [U2,V2] = ID_Cheby(nts,ts2,nu,da,db,tol,1,R_or_N,tR,mR);
    [U3,V3] = ID_Cheby(nts,ts3,nu,da,db,tol,1,R_or_N,tR,mR);
    %U = diag(sqrt(wghts))*U;
elseif opt < 0
    [U1,V1] = ID_Cheby(nts,ts1,nu,da,db,tol,-1,R_or_N,tR,mR);
    [U2,V2] = ID_Cheby(nts,ts2,nu,da,db,tol,-1,R_or_N,tR,mR);
    [U3,V3] = ID_Cheby(nts,ts3,nu,da,db,tol,-1,R_or_N,tR,mR);
    %U = diag(sqrt(wghts))*U;
end
rank1 = size(U1,2);
rank2 = size(U2,2);
rank3 = size(U3,2);
V1 = [zeros(it,rank1);V1];
V2 = [zeros(it,rank2);V2];
V3 = [zeros(it,rank2);V3];
UU23 = kron(U2,U3);
VV23 = kron(V2,V3);
UU12 = kron(U1,U2);
VV12 = kron(V1,V2);
UUU = kron(kron(U1,U2),U3);
VVV = kron(kron(V1,V2),V3);
rank = rank1*rank2*rank3;

vals1 = jacrecur(nts,ts1,it-1,da,db);
vals2 = jacrecur(nts,ts2,it-1,da,db);
vals3 = jacrecur(nts,ts3,it-1,da,db);
vals12 = kron(vals1,vals2);
vals23 = kron(vals2,vals3);
vals = kron(kron(vals1,vals2),vals3);

Ints = sparse(eye(nts,nts));
Iit = sparse(eye(it,it));
xs1 = mod(floor(ts1*nts/2/pi),nts)+1;
Sub1 = sparse(Ints(xs1,:));
xs2 = mod(floor(ts2*nts/2/pi),nts)+1;
Sub2 = sparse(Ints(xs2,:));
xs3 = mod(floor(ts3*nts/2/pi),nts)+1;
Sub3 = sparse(Ints(xs3,:));
SSS = kron(kron(Sub1,Sub2),Sub3);


xsub23 = zeros(nts*nts,1);
for i = 1:nts
    for j =1:nts
        xsub23((i-1)*nts+j) = (xs2(i)-1)*nts + xs3(j);
    end
end
xsub12 = zeros(nts*nts,1);
for i = 1:nts
    for j =1:nts
        xsub12((i-1)*nts+j) = (xs1(i)-1)*nts + xs2(j);
    end
end


if  R_or_N > 0
    fun = @(c)NJacPT3d1(c);
else
    ex = exp(1i*nts/2*ts);
    U = U.*repmat(ex,1,rank);
    fun = @(c)NJacPT3d2(c);
end

    function y = NJacPT3d1(c)
        c = c(:);
        indc = reshape(repmat([0:it-1],it,1),it*it,1)*nts+repmat([1:it]',it,1);
        c1 = c(reshape(repmat([0:it-1],it*it,1),it*it*it,1)*nts^2+repmat(indc,it,1));
        y = vals*c1;
        
        %z1 = kron(kron([vals1 zeros(nts,nts-it)],[vals2 zeros(nts,nts-it)]),[vals3 zeros(nts,nts-it)])*c;
        %e1 = norm(z1-y)/norm(z1)
        
        indc = [1:it*nts]';
        c1 = c(reshape(repmat([0:it-1],it*nts,1),it*it*nts,1)*nts^2+repmat(indc,it,1));
        d = kron(kron(ones(it,1),ones(it,1)),V3).*repmat(c1,1,rank3);
        fft3c = zeros(nts*it^2,rank3);
        for i = 1:it^2
            d1 = nts*ifft(d((i-1)*nts+1:i*nts,:));
            fft3c((i-1)*nts+1:i*nts,:) = d1(xs3,:);
        end
        y2 = zeros(nts^3,rank3);
        for i = 1:nts^2
            for j = 1:it^2
                y2((i-1)*nts+1:i*nts,:) = y2((i-1)*nts+1:i*nts,:) + vals12(i,j)*U3.*fft3c((j-1)*nts+1:j*nts,:);
            end
        end
        y = y + sum(real(y2),2);
        
        %F3 = exp(1i*2*pi/nts*(xs3-1)*[0:nts-1]);
        %FF3 = zeros(nts,nts);
        %for i = 1:rank3
        %    FF3 = FF3 + diag(U3(:,i))*F3*diag(V3(:,i));
        %end
        %z2 = kron(kron([vals1 zeros(nts,nts-it)],[vals2 zeros(nts,nts-it)]),FF3)*c;
        %e2 = norm(z2-sum(y2,2))/norm(z2)
        
        indc = reshape(repmat([0:nts-1],it,1),it*nts,1)*nts+repmat([1:it]',nts,1);
        c1 = c(reshape(repmat([0:it-1],it*nts,1),it*it*nts,1)*nts^2+repmat(indc,it,1));
        d = kron(kron(ones(it,1),V2),ones(it,1)).*repmat(c1,1,rank2);
        fft3c = zeros(nts*it^2,rank2);
        for i = 1:it
            d1 = d((i-1)*nts*it+1:i*nts*it,:);
            sl = [1:it:it*nts]';
            fft3c1 = zeros(nts*it,rank2);
            for j = 1:it
                d2 = nts*ifft(d1(sl+(j-1),:));
                fft3c1(sl+(j-1),:) = d2(xs2,:);
            end
            fft3c((i-1)*nts*it+1:i*nts*it,:) = fft3c1;
        end
        y2 = zeros(nts^3,rank2);
        for i = 1:nts
            for j = 1:it
                y3 = zeros(nts^2,rank2);
                fft3c1 = fft3c((j-1)*nts*it+1:j*nts*it,:);
                for k = 1:nts
                    y3((k-1)*nts+1:k*nts,:) = repmat(U2(k,:),nts,1).*(vals3*fft3c1((k-1)*it+1:k*it,:));
                end
                y2((i-1)*nts^2+1:i*nts^2,:) = y2((i-1)*nts^2+1:i*nts^2,:) + vals1(i,j)*y3;
            end
        end
        y = y + sum(real(y2),2);
        
        %F2 = exp(1i*2*pi/nts*(xs2-1)*[0:nts-1]);
        %FF2 = zeros(nts,nts);
        %for i = 1:rank2
        %    FF2 = FF2 + diag(U2(:,i))*F2*diag(V2(:,i));
        %end
        %z3 = kron(kron([vals1 zeros(nts,nts-it)],FF2),[vals3 zeros(nts,nts-it)])*c;
        %e3 = norm(z3-sum(y2,2))/norm(z3)
        
        c1 = c(1:it*nts^2);
        d = kron(ones(it,1),VV23).*repmat(c1,1,rank2*rank3);
        fft3c = zeros(it*nts^2,rank2*rank3);
        for i = 1:it
            d2 = reshape(d((i-1)*nts^2+1:i*nts^2,:),nts,nts,rank2*rank3);
            d1 = nts^2*ifft2(d2);
            d1 = reshape(d1,nts^2,rank2*rank3);
            fft3c((i-1)*nts^2+1:i*nts^2,:) = d1(xsub23,:);
        end
        y2 = zeros(nts^3,rank2*rank3);
        for i = 1:nts
            for j = 1:it
                y2((i-1)*nts^2+1:i*nts^2,:) = y2((i-1)*nts^2+1:i*nts^2,:) + vals1(i,j)*UU23.*fft3c((j-1)*nts^2+1:j*nts^2,:);
            end
        end
        y = y + sum(real(y2),2);
        
        %z4 = kron(kron([vals1 zeros(nts,nts-it)],FF2),FF3)*c;
        %e4 = norm(z4-sum(y2,2))/norm(z4)
        
        indc = reshape(repmat([0:it-1],it,1),it*it,1)*nts+repmat([1:it]',it,1);
        c1 = c(reshape(repmat([0:nts-1],it*it,1),it*it*nts,1)*nts^2+repmat(indc,nts,1));
        d = kron(V1,kron(ones(it,1),ones(it,1))).*repmat(c1,1,rank1);
        sl = [1:it^2:it*it*nts]';
        fft3c = zeros(nts*it*it,rank1);
        for i = 1:it*it
            d1 = nts*ifft(d(sl+(i-1),:));
             fft3c(sl+(i-1),:) = d1(xs1,:);
        end
        y2 = zeros(nts^3,rank1);
        for i = 1:nts 
            y2((i-1)*nts^2+1:i*nts^2,:) = repmat(U1(i,:),nts^2,1).*(vals23*fft3c((i-1)*it^2+1:i*it^2,:));
        end
        y = y + sum(y2,2);
        
        %F1 = exp(1i*2*pi/nts*(xs1-1)*[0:nts-1]);
        %FF1 = zeros(nts,nts);
        %for i = 1:rank3
        %    FF1 = FF1 + diag(U1(:,i))*F1*diag(V1(:,i));
        %end
        %z5 = kron(kron(FF1,[vals2 zeros(nts,nts-it)]),[vals3 zeros(nts,nts-it)])*c;
        %e5 = norm(z5-sum(y2,2))/norm(z5)
        
        c1 = c(reshape(repmat([0:nts-1],it*nts,1),it*nts^2,1)*nts^2+repmat([1:it*nts]',nts,1));
        d = kron(V1,kron(ones(it,1),V3)).*repmat(c1,1,rank1*rank3);
        fft3c = zeros(it*nts^2,rank1*rank3);
        for j = 1:nts
            d1 = d((j-1)*it*nts+1:j*it*nts,:);
            fft3c1 = zeros(it*nts,rank1*rank3);
            for k = 1:it
                d2 = nts*ifft(d1((k-1)*nts+1:k*nts,:));
                fft3c1((k-1)*nts+1:k*nts,:) = d2(xs3,:);
            end
            fft3c((j-1)*nts*it+1:j*nts*it,:) = fft3c1;
        end
        sl = [1:it*nts:it*nts^2]';
        for i = 1:it*nts
            d1 = nts*ifft(fft3c(sl+(i-1),:));
            fft3c(sl+(i-1),:) = d1(xs1,:);
        end
        y2 = zeros(nts^3,rank1*rank3);
        for i = 1:rank1
            for k = 1:rank3
                for j = 1:nts
                    d1 = fft3c((j-1)*it*nts+1:j*it*nts,(i-1)*rank3+k);
                    y3 = zeros(nts^2,1);
                    for jj = 1:nts
                        for jjj = 1:it
                            y3((jj-1)*nts+1:jj*nts,:) = y3((jj-1)*nts+1:jj*nts,:) + vals2(jj,jjj)*U3(:,k).*d1((jjj-1)*nts+1:jjj*nts,:);
                        end
                    end
                    y2((j-1)*nts^2+1:j*nts^2,(i-1)*rank3+k) = U1(j,i)*y3;
                end
            end
        end
        y = y + sum(real(y2),2);
        
        %z6 = kron(kron(FF1,[vals2 zeros(nts,nts-it)]),FF3)*c;
        %e6 = norm(z6-sum(y2,2))/norm(z6)
	%y = y + sum(real(z6),2);
        
        c1 = c(reshape(repmat([0:nts^2-1],it,1),it*nts^2,1)*nts+repmat([1:it]',nts^2,1));
        d = kron(VV12,ones(it,1)).*repmat(c1,1,rank1*rank2);
        fft3c = zeros(it*nts^2,rank1*rank2);
        sl = [1:it:it*nts^2]';
        for i = 1:it
            d2 = reshape(d(sl+(i-1),:),nts,nts,rank1*rank2);
            d1 = nts^2*ifft2(d2);
            d1 = reshape(d1,nts^2,rank1*rank2);
            fft3c(sl+(i-1),:) = d1(xsub12,:);
        end
        y2 = zeros(nts^3,rank1*rank2);
        for i = 1:nts^2 
            y2((i-1)*nts+1:i*nts,:) = repmat(UU12(i,:),nts,1).*(vals3*fft3c((i-1)*it+1:i*it,:));
        end
        y = y + sum(real(y2),2);
        
        %z7 = kron(kron(FF1,FF2),[vals3 zeros(nts,nts-it)])*c;
        %e7 = norm(z7-sum(y2,2))/norm(z7)
        
        d = VVV.*repmat(c,1,rank1*rank2*rank3);
        fft3c = nts^3*ifft3(reshape(d,nts,nts,nts,rank1*rank2*rank3));
        fft3c = reshape(fft3c,nts^3,rank1*rank2*rank3);
        fft3c = SSS*fft3c;
        y2 = UUU.*fft3c;
        y = y + sum(real(y2),2);
	    y = real(y);
        
        %z8 = kron(kron(FF1,FF2),FF3)*c;
        %e8 = norm(z8-sum(y2,2))/norm(z8)
    end

    
end
