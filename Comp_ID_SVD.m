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
       it = 9;
    else
       it = 27;
    end
    
    [avals0,psivals0,r0,ts0,dnus0] = jac_exp_extr(tol,0,nts,da,db);
    J0 = avals0.*exp(1i*(repmat(dnus0.',size(ts0,1),1).*psivals0-ts0*dnus0.'));
    [cosvals1,sinvals1,r1,ts1,dnus1] = jac_exp_extr(tol,1,nts,da,db);
    ind = zeros(size(dnus1,1),1);
    for i = 1:size(dnus1,1)
	ind(i) = find(dnus0==dnus1(i));
    end
    %J1 = (avals0(:,ind).*exp(1i*(repmat(dnus1.',size(ts1,1),1).*psivals0(:,ind)-ts1*dnus1.')))*r1;
    J1 = J0(:,ind)*r1;
    rank1(ii) = size(r1,1);
    %size(avals0)
    norm(J0-J1)/norm(J1)
    tic
    for i = 1:num
        [U S V] = svdtrunc(J0,rank1(ii),tol);
    end
    rank0(ii) = size(S,1);
    timefac(ii) = toc/num;
    
    %nn = min(size(J1,1),size(J0,1));
    %n1 = randsample(nn,m);
    error01(ii)=norm(J1-U*S*V')/norm(J1);
    fprintf('\n  %-5d %-9d  %-9d   %-1.6E   %-1.6E\n',m,rank0(ii),rank1(ii),error01(ii),timefac(ii));

end
   k1 = 16;
   k2 = 24;
   chebdata1 = cos([k1-1:-1:0]*pi/(k1-1))';
   chebdata2 = cos([k2-1:-1:0]*pi/(k2-1))';
   cd0 = zeros(2,1000);
   nintscd = 0;
   if nts < 2^12
      nintscd        =  nintscd+1;
      cd0(1,nintscd) =  27.0;
      cd0(2,nintscd) =  81.0;
   else
      nintscd        =  nintscd+1;
      cd0(1,nintscd) =  9.0;
      cd0(2,nintscd) =  27.0;
   end


   while cd0(2,nintscd) < nts
      nintscd        = nintscd+1;
      cd0(1,nintscd) =  cd0(2,nintscd-1);
      cd0(2,nintscd) =  cd0(2,nintscd-1)*3.0;

   end

   cd0(2,nintscd) = min(nts,cd0(2,nintscd));

   cd = cd0(:,1:nintscd); 
   dnus = zeros(nintscd*k2,1);
   for intcd=1:nintscd
      c = cd(1,intcd);
      d = cd(2,intcd);

      for i=1:k2
          dnu                = chebdata2(i)*(d-c)/2 + (d+c)/2;
          idx                = i+(intcd-1)*k2;
          dnus(idx)  = dnu;
      end
   end

   dd         = 1/nts;
   dd         = min(0.01,dd);

   dd         = log(dd)/log(2);
   nints      = ceil(-dd)+1;
   nintsab    = 2*nints;

   nints0  = nintsab/2;
   dd      = 2.0;
   
   ab = zeros(2,nintsab);
   for int=1:nints0
      ab(1,int) = dd^(-nints0+int-1);
      ab(2,int) = dd^(-nints0+int);
   end 

   ab = pi/2*ab;

   for int=1:nints0
      ab(1,nintsab-int+1) = pi - ab(2,int);
      ab(2,nintsab-int+1) = pi - ab(1,int);
   end 
   idx = 0;
   ts = zeros(nintsab*k1,1);
   for intab=1:nintsab
       a = ab(1,intab);
       b = ab(2,intab);
       for i=1:k1
           t = chebdata1(i)*(b-a)/2 + (b+a)/2;
           idx = idx+1;
           ts(idx)    = t;
       end
   end

   da = 0.25;
   db = 0.25;
   nt = zeros(nts,1);
   t = getts(nt,da,db);
   w = [1/2 (-1).^[1:k1-2] 1/2*(-1)^(k1-1)]';
   [S,st,nu] = barcycheby(t,chebdata1,w,ab);
   y = interpjac1(nt,t,dnus0(1),da,db,-1);
   z = S*J0(:,1);
   norm(y-z)/norm(y)

   t1 = unique(rand(nts,1)*(pi-2/nts)+1/nts);
   p = @(x)exp(1i*x*pi+pi/2);
   S = barcycheby(t1,chebdata1,w,ab);
   y = S*p(ts0);
   z = p(t1);
   norm(y-z)/norm(z)

