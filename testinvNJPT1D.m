format long
num=10;
da=0.25;
db=0.25;
tol=1e-12
str1='size';
str2='solved?';
str3='iter';
str4='error_our';
str5='time';
fprintf('\n');
fprintf('start 1D inverse nonuniform Jacobi polynomial transform test:');
fprintf('\n');
fprintf('da = %1.2f,db = %1.2f\n',da,db);
fprintf('%-6s%-11s%-15s%-15s%-15s\n',str1,str2,str3,str4,str5);
for m=7:13
    nts=2^m;
    ts = unique(rand(nts,1)*(pi-2/nts)+1/nts);
    b = randn(nts,1);
    [yfun,fun_pre] = invNJPT1D(nts,ts,da,db,tol);
    tic
    for i = 1:num
        [y,flag,relres,iter] = yfun(b);
    end
    time = toc/num;
    fprintf('\n %-5d %-9d %-1.6E %-1.6E %-1.6E\n',m,flag,iter,relres,time);
end