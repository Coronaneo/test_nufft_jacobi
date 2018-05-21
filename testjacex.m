nts=2^9
tR=60;
mR=60;
tol=1e-12
[ts jacobi1 jacobi2]=jacobiexample(nts);
jacobi1=[zeros(nts,27) jacobi1];
jacobi2=[zeros(nts,27) jacobi2];
[U1,V1]=lowrank(jacobi1,tol,tR,mR);
[U2,V2]=lowrank(jacobi2,tol,tR,mR);
rank1=size(U1,2)
rank2=size(U2,2)
c=rand(nts,1);
ncol = size(c,2);


xs=mod(round(ts*nts/2/pi),nts)+1;
d = repmat(conj(V2),1,ncol).*reshape(repmat(c,rank2,1),nts,rank2*ncol);
fftc = ifft(d);
fftc = fftc(xs,:);
result2 = nts*squeeze(sum(reshape(repmat(U2,1,ncol).*fftc,nts,rank2,ncol),2));

ex = exp(1i*nts/2*ts);
result1=zeros(nts,1)
do i=1:rank1
   cj = nufft1dIInyumex(ts,1,tol,conj(V1(:,i))*c);
   result1 = result1 + U1(:,i).*ex.*cj;
end

k=[0:nts-1];
F2=exp(1i*(ts.')*k);
result3=(jacobi1.*F2)*c;

error13=norm(result1-result3)/norm(result3)
error23=norm(result2-result3)/norm(result3)