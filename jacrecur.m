function [vals] = jacrecur(nts,ts,n,da,db)

%   Evaluate the functions (1) for dnu=0,1,...,n at a specified point t
%   using the well-known recurrence relations.


dnu = n;
vals = zeros(nts,n+1);
for ii =1:nts
t    = ts(ii);
x    = cos(t);

vals(ii,1) = sqrt(2^(-1-da-db)*(1+da+db))*sqrt(gamma(1+da+db)/(gamma(1+da)*gamma(1+db)));

if n > 0
vals(ii,2) = (sqrt(2^(-1-da-db)*(3+da+db))*(da-db+(2+da+db)*x)*sqrt(gamma(2+da+db)/(gamma(2+da)*gamma(2+db))))/2;
end

for i=3:n+1

dd1 = (sqrt((-1+da+db+2*i)*(1+da+db+2*i))*(4*(-1+da+db)*i*x+4*i^2*x+(da+db)*(da-db+(-2+da+db)*x)))...
    /(2*sqrt((i*(da+i)*(db+i))/(da+db+i))*(da+db+i)*(-2+da+db+2*i));

dd2 =(2^((da+db)/2)*(-1+da+i)*(-1+db+i)*(da+db+2*i))/(i*(da+db+i)*sqrt(-3+da+db+2*i)*(-2+da+db+2*i)...
    *sqrt((2^(da+db)*(-1+da+i)*(da+i)*(-1+db+i)*(db+i))/((-1+i)*i*(-1+da+db+i)*(da+db+i)*(1+da+db+2*i))));

vals(ii,i) = dd1*vals(ii,i-1) - dd2*vals(ii,i-2);


end

%  Scale by r(t) now
rval       = 2^((1+da+db)/2) * cos(t/2)^(db+0.5) * sin(t/2)^(da+0.5);
vals(ii,:) = vals(ii,:) * rval;


end

end