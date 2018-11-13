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
fprintf('start 1D nonuniform Jacobi polynomial transform test:');
fprintf('\n');
fprintf('da = %1.2f,db = %1.2f\n',da,db);
fprintf('%-6s%-11s%-15s%-15s%-15s\n',str1,str2,str7,str4,str9);
%funour = @(rs,cs,n,da,db,ts,nu)funour1d(rs,cs,n,da,db,ts,nu);
vd = [7:16];
es = length(vd);
rank = zeros(es,1);
errorour = zeros(es,1);
timeour = zeros(es,1);
for ii=1:es
    m = vd(ii);
    nts=2^m;
    if nts < 2^12
       it = 9;
    else
       it = 27;
    end
    
    nt=zeros(nts,1);
    c = randn(nts,1);

    ts = unique(rand(nts,1)*(pi-2/nts)+1/nts);
    wghts = ones(nts,1);
    nu = [0:nts-1]';
    n1 = randsample(nts,m);

    d = c;
    tic;
    
    [result3,t]=directjac1(nt,d,da,db,n1,ts,nu,wghts);  
    %result3 = real(result3);
    timedir = nts/m*(toc-t);

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
 
    [fun,rank(ii)] = NJPT1D(nts,ts,da,db,tR,mR,tol);
    
    tic;
    for j=1:num
        result2 = fun(c);
    end
    timeour(ii)=toc/num;
       
    
    errorour(ii)=norm(result2(n1)-result3)/norm(result3);
    fprintf('\n  %-5d %-9d  %-1.6E   %-1.6E   %-1.6E\n',m,rank(ii),errorour(ii),timeour(ii),timedir);

end
   figure('visible','off');
    pic = figure;
    hold on;
    h(1) = plot(vd,vd-vd(1)+timeour(1),'--k','LineWidth',2);
    h(2) = plot(vd,2*vd-vd(1)*2-timeour(1),'--b','LineWidth',2);
    h(3) = plot(vd,log2(timeour),'-^r','LineWidth',2);
    legend('N log(N)','N^2','timeapp','Location','NorthWest');
    axis square;
    xlabel('log_2(N)'); ylabel('log_{2}(time)');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    saveas(pic,['testNJPT1D.eps'],'epsc');
