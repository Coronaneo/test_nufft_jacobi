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
fprintf('start Chebyshev 3D transform test:');
fprintf('\n');
fprintf('da = %1.2f,db = %1.2f\n',da,db);
fprintf('%-6s%-11s%-11s%-15s%-15s%-15s%-15s%-14s%-10s\n',str1,str2,str3,str4,str5,str6,str7,str8,str9);
funnyu = @(rs,cs,n,da,db)funnyu3d(rs,cs,n,da,db);
funour = @(rs,cs,n,da,db)funour3d(rs,cs,n,da,db);
for m=5:10
    nts=2^m;
    if nts < 2^12
       it = 27;
    else
       it = 9;
    end
    
    nt=zeros(nts,1);
    c = randn(nts^3,1);
%    nn=[nts,0]';
%    [ts,jacobi1,jacobi2] = jacobiexample(nt,da,db);
%    jacobi1=[zeros(nts,it) jacobi1];
%    kk=[0:nts-1];
%    ee=exp(1i*ts*kk);
    %size(ee)
%    jacobi1=jacobi1.*ee;
%    cheb_dir=kron(kron(jacobi1,jacobi1),jacobi1);
%    result4=cheb_dir*c;
%    norm(result4)
%    jacobi2=[zeros(nts,it) jacobi2];
%    cheb2_our=kron(jacobi2,jacobi2);
%    cheb2_nyu=kron(jacobi1,jacobi1);

%    nn = log(nts)/log(2);
%    n1 = (randsample(nts-it,nn)+it-1)*1.000;
%    n2 = (randsample(nts-it,nn)+it-1)*1.000;
%    n3 = (randsample(nts-it,nn)+it-1)*1.000;
    
    n1 = randsample(nts,m);
    n2 = randsample(nts,m);
    n3 = randsample(nts,m);
     
    d=zeros((nts-it)^3,1);
    for p=1:nts-it
        for q=1:nts-it
            d((p-1)*(nts-it)^2+(q-1)*(nts-it)+1:(p-1)*(nts-it)^2+q*(nts-it))=c(it*nts^2+(p-1)*nts^2+it*nts+(q-1)*nts+it+1:it*nts^2+(p-1)*nts^2+it*nts+q*nts);
        end
    end
%    norm(d)
    tic;
    for i=1:5
    [result3,ts,t]=directjac3(nt,d,da,db,n1,n2,n3);
    end
    timedir = nts^3/m^3*(toc-t)/5;
%    size(result3)
%    tic;
%    size(d)
%    d(1:5)
%    for i=1:2
%    [result3,~]=directjac3(nt,d,da,db,n1,n2,n3);
%    end
%    size(result3)
%    result3(1:10)
%    ier
%    timedir=((nts-it)/nn)^3*(toc/2);
%    norm(result3)
%    errordir=norm(result4-result3)/norm(result4)

    [ts1,ts2,ts3]=ndgrid(ts,ts,ts);
    ts=[ts1(:) ts2(:) ts3(:)];
    xs=mod(floor(ts*nts/2/pi),nts)+1;
    xs = sub2ind([nts nts nts],xs(:,1),xs(:,2),xs(:,3));
    s=round(nts*ts);
    gamma=norm(nts*ts-s,inf);
    xi=log(log(10/tol)/gamma/7);
    lw=xi-log(xi)+log(xi)/xi+0.5*log(xi)^2/xi^2-log(xi)/xi^2;
    if m<7
       K=ceil(30*gamma*exp(lw));
    elseif m<8
       K=ceil(40*gamma*exp(lw));
    elseif m<9
       K=ceil(50*gamma*exp(lw));
    elseif m<10
       K=ceil(60*gamma*exp(lw));
    else
       K=ceil(70*gamma*exp(lw));
    end
    tR=K+2;
    mR=K;
 
%    cs=[nts*it+it:nts*it+it+100]*1.00;
%    rs=[1:5]*1.00;
%    [M,ier]=extrcheb2(nt,rs,cs,1);  
%    com=norm(cheb2_our(1:5,nts*it+it:nts*it+it+100)-M)/norm(M)
   % M
%    ier

    [U1,V1]=lowrank(nts^3,funnyu,da,db,tol,tR,mR);
    [U2,V2]=lowrank(nts^3,funour,da,db,tol,3*tR,3*mR);
    rank1=size(U1,2);
    %V1=[zeros(nts*it,rank1);V1];
    rank2=size(U2,2);
    ncol = size(c,2);
%    U1(1:5)

    

    tic;
    for j=1:num
        d = reshape( repmat(conj(V2),1,ncol).*reshape(repmat(c,rank2,1), nts^3, rank2*ncol), nts,nts,nts,rank2*ncol);
        fft3c = ifft3(d);
        fft3c = reshape(fft3c,nts^3,rank2*ncol);
        fft3c = fft3c(xs,:);
        result2 = nts^3*squeeze(sum(reshape(repmat(U2,1,ncol).*fft3c,nts^3,rank2,ncol),2));
    end
    timeour=toc/num;
%norm(result2)
    ex = exp(1i*nts/2*(ts(:,1)+ts(:,2)+ts(:,3)));
    U1=U1.*repmat(ex,1,rank1);
    tic;
    for j=1:num
        result1=zeros(nts^3,1);
        for i=1:rank1
            cj = nufft3dIInyumex(ts(:,1),ts(:,2),ts(:,3),1,tol,reshape(conj(V1(:,i)).*c,nts,nts,nts));
            result1 = result1 + U1(:,i).*cj;
        end
    end
    timenyu=toc/num;
    timeratio=timeour/timenyu;
%norm(result1)
    
   
    n1=repmat((n1'-1)*nts^2,m^2,1);
    n2=remat(repmat((n2'-1)*nts,m,1),1,m);
    n3=repmat(n3,m^2,1);
    n=n1(:)+n2(:)+n3(:);       
    
    
%    error1=norm(result1-result2)/norm(result2)
    errornyu=norm(result1(n)-result3)/norm(result3);
    errorour=norm(result2(n)-result3)/norm(result3);
    fprintf('\n   %-5d %-9d  %-9d  %-1.6E   %-1.6E   %-1.6E   %-1.6E   %-1.6E  %-1.6E\n',m,rank2,rank1,timeour,timenyu,timeratio,errorour,errornyu,timedir);
%    gc=imagesc(real(jacobi1(:,it+1:end)));
%    saveas(gc,'image13.jpg');
%    gf=imagesc(real(jacobi1(:,it+1:end)*1i));
%    saveas(gf,'image13i.jpg');
%    bf=imagesc(abs(jacobi1(:,it+1:end)));
%    saveas(bf,'image13a.jpg');
end
