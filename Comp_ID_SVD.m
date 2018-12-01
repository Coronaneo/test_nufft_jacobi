format long
num=10;
da=0.25;
db=0.25;
tol=1e-8
str1='size';
str2='rank0';
str3='rank1';
str4='error';
str5='time_svdfac';
fprintf('\n');
fprintf('start compare:');
fprintf('\n');
fprintf('da = %1.2f,db = %1.2f\n',da,db);
fprintf('%-6s%-11s%-15s%-15s%-15s%-15s\n',str1,str2,str3,str4,str5);
%funour = @(rs,cs,n,da,db,ts,nu)funour1d(rs,cs,n,da,db,ts,nu);
vd = [7:16];
es = length(vd);
rank0 = zeros(es,1);
rank1 = zeros(es,1);
error01 = zeros(es,1);
timefac = zeros(es,1);
for ii=1:es
    m = vd(ii);
    nts=2^m;
    if nts < 2^12
       it = 10;
    else
       it = 28;
    end
    
    [cosvals1,sinvals1,r1] = jac_exp_extr(tol,1,nts,da,db);
    J1 = (cosvals1+1i*sinvals1)*r1;
    rank1(ii) = size(r1,1);
    size(cosvals1)   
    [cosvals0,sinvals0,r0] = jac_exp_extr(tol,0,nts,da,db);
    J0 = cosvals0+1i*sinvals0;
    size(cosvals0)
    norm(J0-J1)/norm(J1)
    tic
    for i = 1:num
        [U S V] = svdtrunc(J0,rank(ii),tol);
    end
    rank0(ii) = size(S,1);
    timefac(ii) = toc/num;
    
    nn = min(size(J1,1),size(J0,1));
    n1 = randsample(nn,m);
    error01(ii)=norm(J1-U*S*V')/norm(J1);
    fprintf('\n  %-5d %-9d  %-9d   %-1.6E   %-1.6E\n',m,rank0(ii),rank1(ii),error01(ii),timefac(ii));

end
    
