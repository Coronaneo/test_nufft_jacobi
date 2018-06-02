    function M = funnyu2d(rs,cs,n)
       rs=rs*1.000;
       cs=cs*1.000;
       ns=zeros(sqrt(n),1);
       M=extrcheb2(ns,rs,cs,0);
    end
