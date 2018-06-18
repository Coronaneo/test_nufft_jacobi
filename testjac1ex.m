format long
num=20;
da=0;
db=0;
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
fprintf('start Jacobi 1D transform test:');
fprintf('\n');
fprintf('da = %1.2f,db = %1.2f\n',da,db);
fprintf('%-6s%-11s%-11s%-11s%-15s%-15s%-15s%-15s%-15s%-14s%-10s\n',str1,str10,str2,str3,str4,str5,str6,str11,str7,str8,str9);
funnyu = @(rs,cs,n,da,db)funnyu1d(rs,cs,n,da,db);
funour = @(rs,cs,n,da,db)funour1d(rs,cs,n,da,db);
for m=7:30
    nts=2^m;
    if nts < 2^12
       it = 27;
    else
       it = 9;
    end
    
    nt=zeros(nts,1);
    c = randn(nts,1);
%    nn=[nts,0]';
%    [ts,jacobi1,jacobi2] = jacobiexample(nt,da,db);
%    jacobi1=[zeros(nts,it) jacobi1];
%    jacobi2=[zeros(nts,it) jacobi2];
%    cheb2_our=kron(jacobi2,jacobi2);
%    cheb2_nyu=kron(jacobi1,jacobi1);

%    nn = floor(nts^(1/2));
%    n1 = (randsample(nts-it,nn)+it-1)*1.000;
    
    n1 = randsample(nts,m);

    d = c(it+1:end);
    tic;
    for i=1:5
    [result3,ts]=directjac1(nt,d,da,db,n1);
    end
    timedir = nts/m*toc/5;

%    tic;
%    size(d)
%    d(1:5)i
%    size(result3)
%    for i=1:5
%    [result3,~]=directjac1(nt,d,da,db,n1);
%    end
%    size(result3)
%    result3(1:10)
%    ier
%    timedir=(nts-it)/nn*toc/5;


    xs=mod(floor(ts*nts/2/pi),nts)+1;
    s=round(nts*ts);
    gamma=norm(nts*ts-s,inf);
    xi=log(log(10/tol)/gamma/7);
    lw=xi-log(xi)+log(xi)/xi+0.5*log(xi)^2/xi^2-log(xi)/xi^2;
    if m<10
       K=ceil(10*gamma*exp(lw));
    elseif m<14
       K=ceil(11*gamma*exp(lw));
    elseif m<18
       K=ceil(12*gamma*exp(lw));
    elseif m<21
       K=ceil(13*gamma*exp(lw));
    elseif m<24
       K=ceil(14*gamma*exp(lw));
    elseif m<27
       K=ceil(15*gamma*exp(lw));
    else
       K=ceil(17*gamma*exp(lw));
    end
    tR=K+2;
    mR=K;
 
%    cs=[nts*it+it:nts*it+it+100]*1.00;
%    rs=[1:5]*1.00;
%    [M,ier]=extrcheb2(nt,rs,cs,1);  
%    com=norm(cheb2_our(1:5,nts*it+it:nts*it+it+100)-M)/norm(M)
   % M
%    ier

    [U1,V1]=lowrank(nts,funnyu,da,db,tol,tR,mR);
    [U2,V2]=lowrank(nts,funour,da,db,tol,tR,mR);
    rank1=size(U1,2);
    %V1=[zeros(nts*it,rank1);V1];
    rank2=size(U2,2);
    ncol = size(c,2);
%    U1(1:5)

    

    tic;
    for j=1:num
        d = repmat(conj(V2),1,ncol).*reshape(repmat(c,rank2,1),nts,rank2*ncol);
        fftc = ifft(d);
        fftc = fftc(xs,:);
        result2 = nts*squeeze(sum(reshape(repmat(U2,1,ncol).*fftc,nts,rank2,ncol),2));
    end
    timeour=toc/num;

    ex = exp(1i*nts/2*ts);
    U1=U1.*repmat(ex,1,rank1);
    tic;
    for j=1:num
        result1=zeros(nts,1);
        for i=1:rank1
            cj = nufft1dIInyumex(ts,1,tol,conj(V1(:,i)).*c);
            result1 = result1 + U1(:,i).*cj;
        end
    end
    timenyu=toc/num;
    timeratio=timeour/timenyu;

    [r,expvals,tss] = chebjacex(nt,da,db,tol);
    %r(1:it,:)=0;
    rank3 = size(r,2);
    xs=mod(floor(tss*nts/2/pi),nts)+1;
    b = repmat(r,1,ncol).*reshape(repmat(c,rank3,1),nts,rank3*ncol);     
    fftb = ifft(b);
    fftb = fftb(xs,:);
    result4 = nts*squeeze(sum(reshape(repmat(expvals,1,ncol).*fftb,nts,rank3,ncol),2));
    errorcheb = norm(result4(n1)-result3)/norm(result3);
    
          
    
    
%    error1=norm(result1-result2)/norm(result2)
    errornyu=norm(result1(n1)-result3)/norm(result3);
    errorour=norm(result2(n1)-result3)/norm(result3);
    fprintf('\n   %-5d %-9d  %-9d  %-9d  %-1.6E   %-1.6E   %-1.6E   %-1.6E   %-1.6E   %-1.6E  %-1.6E\n',m,rank3,rank2,rank1,timeour,timenyu,timeratio,errorcheb,errorour,errornyu,timedir);
%    gc=imagesc(real(jacobi1(:,it+1:end)));
%    saveas(gc,'image13.jpg');
%    gf=imagesc(real(jacobi1(:,it+1:end)*1i));
%    saveas(gf,'image13i.jpg');
%    bf=imagesc(abs(jacobi1(:,it+1:end)));
%    saveas(bf,'image13a.jpg');
end
