function [fun,rank1,rank2,rank3] = invJPT3D(nts,da,db,tR,mR,tol,opt,R_or_N)
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
[ts,wghts] = getts(nt,da,db);
nu = [it:nts-1]';
ts1 = ts;
ts2 = ts;
ts3 = ts;
%xs = mod(floor(ts*nts/2/pi),nts)+1;


if opt >= 1
    JTM1 = @(ts,nu)interpjac1(nt,ts,nu,da,db,R_or_N);
    [U1,V1] = lowrank(nts,JTM1,ts1,nu,tol,tR,mR);
    U1 = diag(sqrt(wghts))*U1;
    V1 = conj(V1);
    JTM2 = @(ts,nu)interpjac1(nt,ts,nu,da,db,R_or_N);
    [U2,V2] = lowrank(nts,JTM2,ts2,nu,tol,tR,mR);
    U2 = diag(sqrt(wghts))*U2;
    V2 = conj(V2);
    JTM3 = @(ts,nu)interpjac1(nt,ts,nu,da,db,R_or_N);
    [U3,V3] = lowrank(nts,JTM1,ts3,nu,tol,tR,mR);
    U3 = diag(sqrt(wghts))*U3;
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
V3 = [zeros(it,rank3);V3];
UU23 = kron(U2,U3);
VV23 = kron(V2,V3);
UU12 = kron(U1,U2);
VV12 = kron(V1,V2);
UUU = kron(kron(U1,U2),U3);
VVV = kron(kron(V1,V2),V3);
rank = rank1*rank2*rank3;

vals1 = jacrecur(nts,ts1,it-1,da,db).';
vals1 = vals1*diag(sqrt(wghts));
vals2 = jacrecur(nts,ts2,it-1,da,db).';
vals2 = vals2*diag(sqrt(wghts));
vals3 = jacrecur(nts,ts3,it-1,da,db).';
vals3 = vals3*diag(sqrt(wghts));
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

Y = [1:nts]';
X = zeros(nts,1);
S = ones(nts,1);
for i = 1:nts
    X(i) = xs1(i);   
end
P1 = sparse(X,Y,S,nts,nts,nts);
P2 = P1;
P3 = P1;
PP23 = kron(P2,P3);
PP12 = kron(P1,P2);
PPP = kron(PP12,P3);

if  R_or_N > 0
    fun = @(c)invJacPT3d1(c);
else
    ex = exp(1i*nts/2*ts);
    U = U.*repmat(ex,1,rank);
    fun = @(c)invJacPT3d2(c);
end

    function y = invJacPT3d1(c)
        c = c(:);
        c = kron(kron(1./sqrt(wghts),1./sqrt(wghts)),1./sqrt(wghts)).*c;
        ind1 = reshape(repmat([0:it-1],it,1),it^2,1)*nts+repmat([1:it]',it,1);
        ind = reshape(repmat([0:it-1],it^2,1),it^3,1)*nts^2+repmat(ind1,it,1);
        y = zeros(nts^3,1);
        y(ind,:) = vals*c;
        
        %z1 = kron(kron([vals1;zeros(nts-it,nts)],[vals2;zeros(nts-it,nts)]),[vals3;zeros(nts-it,nts)])*c;
        %e1 = norm(z1-y)/norm(z1)
        
        ind1 = [1:it*nts]';
        ind = reshape(repmat([0:it-1],it*nts,1),it^2*nts,1)*nts^2+repmat(ind1,it,1);
        d = zeros(nts*it^2,rank3);
        C = repmat(c,rank3);
        for i = 1:it^2
            for j = 1:nts^2
                d((i-1)*nts+1:i*nts,:) = d((i-1)*nts+1:i*nts,:) + vals12(i,j)*U3.*C((j-1)*nts+1:j*nts,:);
            end
        end
        fft3c = zeros(nts*it^2,rank3);
        for i = 1:it^2
            d1 = P3*d((i-1)*nts+1:i*nts,:);
            fft3c((i-1)*nts+1:i*nts,:) = conj(fft(conj(d1)));
        end
        y2 = kron(kron(ones(it,1),ones(it,1)),V3).*fft3c;
        y(ind,:) = y(ind,:) + sum(real(y2),2);
        
        %F3 = exp(1i*2*pi/nts*(xs3-1)*[0:nts-1]);
        %FF3 = zeros(nts,nts);
        %for i = 1:rank3
        %    FF3 = FF3 + diag(V3(:,i))*F3.'*diag(U3(:,i));
        %end
        %z2 = kron(kron([vals1;zeros(nts-it,nts)],[vals2;zeros(nts-it,nts)]),FF3)*c;
        %e2 = norm(z2(ind,:)-sum(y2,2))/norm(z2(ind,:))
        
        ind1 = reshape(repmat([0:nts-1],it,1),it*nts,1)*nts+repmat([1:it]',nts,1);
        ind = reshape(repmat([0:it-1],it*nts,1),it*it*nts,1)*nts^2+repmat(ind1,it,1);
        d = zeros(nts*it^2,rank2);
        C = repmat(c,1,rank2);
        for i = 1:it
            for j = 1:nts
                d1 = zeros(nts*it,rank2);
                C1 = C((j-1)*nts*nts+1:j*nts*nts,:);
                for k = 1:nts
                    d1((k-1)*it+1:k*it,:) = repmat(U2(k,:),it,1).*(vals3*C1((k-1)*nts+1:k*nts,:));
                end
                d((i-1)*nts*it+1:i*nts*it,:) = d((i-1)*nts*it+1:i*nts*it,:) + vals1(i,j)*d1;
            end
        end
        fft3c = zeros(nts*it^2,rank2);
        for i = 1:it
            d1 = d((i-1)*nts*it+1:i*nts*it,:);
            sl = [1:it:it*nts]';
            fft3c1 = zeros(nts*it,rank2);
            for j = 1:it
                d2 = P2*d1(sl+(j-1),:);
                fft3c1(sl+(j-1),:) = conj(fft(conj(d2)));
            end
            fft3c((i-1)*nts*it+1:i*nts*it,:) = fft3c1;
        end
        y2 = kron(kron(ones(it,1),V2),ones(it,1)).*fft3c;
        y(ind,:) = y(ind,:) + sum(real(y2),2);
        
        %F2 = exp(1i*2*pi/nts*(xs2-1)*[0:nts-1]);
        %FF2 = zeros(nts,nts);
        %for i = 1:rank2
        %    FF2 = FF2 + diag(V2(:,i))*F2.'*diag(U2(:,i));
        %end
        %z3 = kron(kron([vals1;zeros(nts-it,nts)],FF2),[vals3;zeros(nts-it,nts)])*c;
        %e3 = norm(z3(ind,:)-sum(y2,2))/norm(z3(ind,:))
        
        ind = [1:it*nts^2]';
        d = zeros(nts^2*it,rank2*rank3);
        C = repmat(c,1,rank2*rank3);
        for i = 1:it
            for j = 1:nts
                d((i-1)*nts^2+1:i*nts^2,:) = d((i-1)*nts^2+1:i*nts^2,:) + vals1(i,j)*UU23.*C((j-1)*nts^2+1:j*nts^2,:);
            end
        end 
        fft3c = zeros(it*nts^2,rank2*rank3);
        for i = 1:it
            d2 = PP23*d((i-1)*nts^2+1:i*nts^2,:);
            d2 = reshape(d2,nts,nts,rank2*rank3);
            d1 = conj(fft2(conj(d2)));
            fft3c((i-1)*nts^2+1:i*nts^2,:) = reshape(d1,nts^2,rank2*rank3);
        end
        y2 = kron(ones(it,1),VV23).*fft3c;
        y(ind,:) = y(ind,:) + sum(real(y2),2);
        
        %z4 = kron(kron([vals1;zeros(nts-it,nts)],FF2),FF3)*c;
        %e4 = norm(z4(ind,:)-sum(y2,2))/norm(z4(ind,:))
        
        ind1 = reshape(repmat([0:it-1],it,1),it*it,1)*nts+repmat([1:it]',it,1);
        ind = reshape(repmat([0:nts-1],it*it,1),it*it*nts,1)*nts^2+repmat(ind1,nts,1);
        d = zeros(nts*it^2,rank1);
        C = repmat(c,1,rank1);
        for i = 1:nts 
            d((i-1)*it^2+1:i*it^2,:) = repmat(U1(i,:),it^2,1).*(vals23*C((i-1)*nts^2+1:i*nts^2,:));
        end
        sl = [1:it^2:it*it*nts]';
        fft3c = zeros(nts*it*it,rank1);
        for i = 1:it*it
            d1 = P1*d(sl+(i-1),:);
            fft3c(sl+(i-1),:) = conj(fft(conj(d1)));
        end
        y2 = kron(V1,kron(ones(it,1),ones(it,1))).*fft3c;
        y(ind,:) = y(ind,:) + sum(y2,2);
        
        %F1 = exp(1i*2*pi/nts*(xs1-1)*[0:nts-1]);
        %FF1 = zeros(nts,nts);
        %for i = 1:rank1
        %    FF1 = FF1 + diag(V1(:,i))*F1.'*diag(U1(:,i));
        %end
        %z5 = kron(kron(FF1,[vals2;zeros(nts-it,nts)]),[vals3;zeros(nts-it,nts)])*c;
        %e5 = norm(z5(ind,:)-sum(y2,2))/norm(z5(ind,:))
        
        ind = reshape(repmat([0:nts-1],it*nts,1),it*nts^2,1)*nts^2+repmat([1:it*nts]',nts,1);
        d = zeros(nts^2*it,rank1*rank3);
        C = repmat(c,1,rank1*rank3);
        for i = 1:rank1
            for k = 1:rank3
                for j = 1:nts
                    d1 = C((j-1)*nts*nts+1:j*nts*nts,(i-1)*rank3+k);
                    d2 = zeros(nts*it,1);
                    for jj = 1:it
                        for jjj = 1:nts
                            d2((jj-1)*nts+1:jj*nts,:) = d2((jj-1)*nts+1:jj*nts,:) + vals2(jj,jjj)*U3(:,k).*d1((jjj-1)*nts+1:jjj*nts,:);
                        end
                    end
                    d((j-1)*nts*it+1:j*nts*it,(i-1)*rank3+k) = U1(j,i)*d2;
                end
            end
        end
        fft3c = zeros(it*nts^2,rank1*rank3);
        for j = 1:nts
            d1 = d((j-1)*it*nts+1:j*it*nts,:);
            fft3c1 = zeros(it*nts,rank1*rank3);
            for k = 1:it
                d2 = P3*d1((k-1)*nts+1:k*nts,:);
                fft3c1((k-1)*nts+1:k*nts,:) = conj(fft(conj(d2)));
            end
            fft3c((j-1)*nts*it+1:j*nts*it,:) = fft3c1;
        end
        sl = [1:it*nts:it*nts^2]';
        for i = 1:it*nts
            d1 = P1*fft3c(sl+(i-1),:);
            fft3c(sl+(i-1),:) = conj(fft(conj(d1)));
        end
        y2 = kron(V1,kron(ones(it,1),V3)).*fft3c;
        y(ind,:) = y(ind,:) + sum(real(y2),2);
        
        %z6 = kron(kron(FF1,[vals2;zeros(nts-it,nts)]),FF3)*c;
        %e6 = norm(z6(ind,:)-sum(y2,2))/norm(z6(ind,:))
	%y = y + sum(real(z6),2);
        
        ind = reshape(repmat([0:nts^2-1],it,1),it*nts^2,1)*nts+repmat([1:it]',nts^2,1);
        d = zeros(nts^2*it,rank1*rank2);
        C = repmat(c,1,rank1*rank2);
        for i = 1:nts^2 
            d((i-1)*it+1:i*it,:) = repmat(UU12(i,:),it,1).*(vals3*C((i-1)*nts+1:i*nts,:));
        end
        fft3c = zeros(it*nts^2,rank1*rank2);
        sl = [1:it:it*nts^2]';
        for i = 1:it
            d2 = PP12*d(sl+(i-1),:);
            d2 = reshape(d2,nts,nts,rank1*rank2);
            d1 = conj(fft2(conj(d2)));
            fft3c(sl+(i-1),:) = reshape(d1,nts^2,rank1*rank2);
        end
        y2 = kron(VV12,ones(it,1)).*fft3c;
        y(ind,:) = y(ind,:) + sum(real(y2),2);
        
        %z7 = kron(kron(FF1,FF2),[vals3;zeros(nts-it,nts)])*c;
        %e7 = norm(z7(ind,:)-sum(y2,2))/norm(z7(ind,:))
        
        d = UUU.*repmat(c,1,rank1*rank2*rank3);
        fft3c = PPP*d;
        fft3c = conj(fft3(conj(reshape(fft3c,nts,nts,nts,rank1*rank2*rank3))));
        fft3c = reshape(fft3c,nts^3,rank1*rank2*rank3);
        y2 = VVV.*fft3c;
        y = y + sum(real(y2),2);
	    y = real(y);
        
        %z8 = kron(kron(FF1,FF2),FF3)*c;
        %e8 = norm(z8-sum(y2,2))/norm(z8)
    end

    
end
