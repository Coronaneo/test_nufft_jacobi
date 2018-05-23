format long
num=10;
da=0.25;
db=0.30;
tol=1e-4
str1='size';
str2='our_rank';
str3='nyu_rank';
str4='our_time';
str5='nyu_time';
str6='ratio';
str7='error_our';
str8='error_nyu';
str9='dir_time';
fprintf('\n');
fprintf('%-6s%-11s%-11s%-15s%-15s%-14s%-15s%-14s%-10s\n',str1,str2,str3,str4,str5,str6,str7,str8,str9);

for m=9:20
    nts=2^m;
    nt=zeros(nts,1);
    nn=[nts,0]';
    [ts,jacobi1,jacobi2] = jacobiexample(nt,da,db);
    %s=round(nts*ts);
    %gamma=norm(nts*ts-s,inf);
    %xi=log(log(10/tol)/gamma/7);
    %lw=xi-log(xi)+log(xi)/xi+0.5*log(xi)^2/xi^2-log(xi)/xi^2;
    %K=ceil(5*gamma*exp(lw));
    tR=3*nts/4;
    mR=3*nts/4;
    jacobi1=[zeros(nts,27) jacobi1];
    jacobi2=[zeros(nts,27) jacobi2];
    [U1,V1]=lowrank(jacobi1,tol,tR,mR);
    [U2,V2]=lowrank(jacobi2,tol,tR,mR);
    rank1=size(U1,2);
    rank2=size(U2,2);
    c=rand(nts,1);
    ncol = size(c,2);


    xs=mod(round(ts*nts/2/pi),nts)+1;
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

    k=[0:nts-1];
    F2=exp(1i*ts*k);
    tic;
    for i=1:num
    result3=(jacobi1.*F2)*c;
    end
    timedir=toc/num;

    errornyu=norm(result1-result3)/norm(result3);
    errorour=norm(result2-result3)/norm(result3);
    fprintf('\n   %-5d %-9d  %-9d  %-1.6E   %-1.6E   %-1.6E  %-1.6E   %-1.6E  %-1.6E\n',m,rank2,rank1,timeour,timenyu,timeratio,errorour,errornyu,timedir);
end
