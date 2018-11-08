format long
num=10;
da=0.25;
db=0.25;
tol=1e-12
str1='size';
str2='our_rank';
str4='our_time';
str7='error_our';
str9='dir_time';
fprintf('\n');
fprintf('start 1D inverse uniform Jacobi polynomial transform test:');
fprintf('\n');
fprintf('da = %1.2f,db = %1.2f\n',da,db);
fprintf('%-6s%-11s%-15s%-15s%-15s\n',str1,str2,str7,str4,str9);
for m=7:13
    nts=2^m;
    if nts < 2^12
       it = 27;
    else
       it = 9;
    end
    
    nt=zeros(nts,1);
    c = randn(nts,1);

    [ts,wghts] = getts(nt,da,db);
    %ts = unique(rand(nts,1)*(pi-2/nts)+1/nts);
    nu = [it:nts-1]';
    n1 = randsample(nts-it,m);
    d = c./wghts;
    
    tic;
    
    [result3,t]=directinvjac1(nt,d,da,db,n1,ts,nu,wghts);
    result3 = real(result3);
%    norm(result3)    
    timedir = nts/m*(toc-t);




    %xs=mod(floor(ts*nts/2/pi),nts)+1;
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


    [fun,rank] = invJPT1D(nts,da,db,tR,mR,tol);

    tic;
    for j=1:num
        result2 = fun(c);
    end
    timeour=toc/num;

    errorour=norm(result2(n1+it)-result3)/norm(result3);
    fprintf('\n  %-5d %-9d  %-1.6E   %-1.6E   %-1.6E\n',m,rank,errorour,timeour,timedir);
    
end