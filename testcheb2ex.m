num=20;
da=-0.50;
db=-0.50;
tol=1e-12
str1='size';
str2='our_rank';
str3='nyu_rank';
str4='our_time';
str5='nyu_time';
str6='ratio_our/nyu';
str7='error_our';
str8='error_nyu';
str9='dir_time';
str10='cheb_rank';
str11='error_cheb';
fprintf('\n');
fprintf('start Chebyshev 2D transform test:');
fprintf('\n');
fprintf('%-6s%-11s%-11s%-15s%-15s%-15s%-15s%-14s%-10s\n',str1,str2,str3,str4,str5,str6,str7,str8,str9);

for m=6:12
    nts=2^m;
    if nts < 2^12
       it = 27;
    else
       it = 9;
    end
    
    nt=zeros(nts,1);
    nn=[nts,0]';
    [ts,jacobi1,jacobi2] = jacobiexample(nt,da,db);
    jacobi1=[zeros(nts,it) jacobi1];
    jacobi2=[zeros(nts,it) jacobi2];
    cheb2_our=kron(jacobi2,jacobi2);
    cheb2_nyu=kron(jacobi1,jacobi1);
    
    [ts1,ts2]=ndgrid(ts);
    ts=[ts1(:) ts2(:)];
    xs=mod(floor(ts*nts/2/pi),nts)+1;
    xs = sub2ind([nts nts],xs(:,1),xs(:,2));
    s=round(nts*ts);
    gamma=norm(nts*ts-s,inf);
    xi=log(log(10/tol)/gamma/7);
    lw=xi-log(xi)+log(xi)/xi+0.5*log(xi)^2/xi^2-log(xi)/xi^2;
    if m<14
       K=ceil(11*gamma*exp(lw));
    elseif m<16
       K=ceil(12*gamma*exp(lw));
    elseif m<18
       K=ceil(13*gamma*exp(lw));
    elseif m<21
       K=ceil(14*gamma*exp(lw));
    else
       K=ceil(15*gamma*exp(lw));
    end
    tR=K+2;
    mR=K;
    [U1,V1]=lowrank(cheb2_nyu,tol,tR,mR);
    [U2,V2]=lowrank(cheb2_our,tol,tR,mR);
    rank1=size(U1,2);
    rank2=size(U2,2);
    c=rand(nts^2,1);
    ncol = size(c,2);

    

    tic;
    for j=1:num
        d = reshape(repmat(conj(V2),1,ncol).*reshape(repmat(c,rank2,1),nts^2,rank2*ncol),nts,nts,rank2*ncol);
        fft2c = ifft2(d);
        fft2c = reshape(fft2c,nts^2,rank2*ncol);
        fft2c = fft2c(xs,:);
        result2 = nts^2*squeeze(sum(reshape(repmat(U2,1,ncol).*fft2c,nts^2,rank2,ncol),2));
    end
    timeour=toc/num;

    ex = exp(1i*nts/2*(ts(:,1)+ts(:,2)));
    U1=U1.*repmat(ex,1,rank1);
    tic;
    for j=1:num
        result1=zeros(nts^2,1);
        for i=1:rank1
            cj = nufft2dIInyumex(ts(:,1),ts(:,2),1,tol,reshape(conj(V1(:,i)).*c,nts,nts));
            result1 = result1 + U1(:,i).*cj;
        end
    end
    timenyu=toc/num;
    timeratio=timeour/timenyu;
            

    [k1 k2]=ndgrid(0:nts-1);
    k = [k1(:) k2(:)];
    F2=exp(1i*ts*k');
    tic;
    for i=1:num
    result3=(cheb2_nyu.*F2)*c;
    end
    timedir=toc/num;
    
    
  
    errornyu=norm(result1-result3)/norm(result3);
    errorour=norm(result2-result3)/norm(result3);
    fprintf('\n   %-5d %-9d  %-9d  %-1.6E   %-1.6E   %-1.6E   %-1.6E   %-1.6E  %-1.6E\n',m,rank2,rank1,timeour,timenyu,timeratio,errorour,errornyu,timedir);
%    gc=imagesc(real(jacobi1(:,it+1:end)));
%    saveas(gc,'image13.jpg');
%    gf=imagesc(real(jacobi1(:,it+1:end)*1i));
%    saveas(gf,'image13i.jpg');
%    bf=imagesc(abs(jacobi1(:,it+1:end)));
%    saveas(bf,'image13a.jpg');
end
