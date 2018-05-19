format long
nts=2^9;
iflag=1
load jacobi1i.bin
load jacobi1r.bin
load jacobi2r.bin
load jacobi2i.bin
load ts.txt

jacobi1r=reshape(jacobi1r,[nts,nts-27]);
jacobi1i=reshape(jacobi1i,[nts,nts-27]);
jacobi2r=reshape(jacobi2r,[nts,nts-27]);
jacobi2i=reshape(jacobi2i,[nts,nts-27]);
jacobi1=jacobi1r+1i*jacobi1i;
jacobi2=jacobi2r+1i*jacobi2i;
jacobi1=[zeros(nts,27) jacobi1];
jacobi2=[zeros(nts,27) jacobi2];
size(jacobi1);
size(jacobi2);
tol=1e-12;
tR=60;
mR=60;
[U1,V1]=lowrank(jacobi1,tol,tR,mR);
[U2,V2]=lowrank(jacobi2,tol,tR,mR);
error1=norm(jacobi1-U1*V1')/norm(jacobi1)
error2=norm(jacobi2-U2*V2')/norm(jacobi2)
jacobi1(1:5,28)
%jacobi1(1:5,28)./jacobi2(1:5,28)
%error3=norm(jacobi1-U2*V2')/norm(jacobi1)
rank1=size(U1,2)
rank2=size(U2,2)
fout=fopen('rank.txt','w');
fprintf(fout,'%8i\n',rank1);
fprintf(fout,'%8i\n',rank2);
fclose(fout);

xs=mod(round(ts*nts/2/pi),nts)+1;
c=rand(nts,1);
fout=fopen('c.txt','w');
fprintf(fout,'%12.26f\n',c);
fclose(fout);
c(1:5)
V1(5:10,1)
ncol = size(c,2);
d=c;
%c=conj(V2).*repmat(c,1,rank2);
c = repmat(conj(V2),1,ncol).*reshape(repmat(c,rank2,1),nts,rank2*ncol);
if iflag < 0
   fftc = fft(c);
else
   fftc = ifft(c);
end
fftc = fftc(xs,:);
%fftc=sum(U2.*fftc,2);
fftc = nts*squeeze(sum(reshape(repmat(U2,1,ncol).*fftc,nts,rank2,ncol),2));
fout=fopen('result2.txt','w');
fprintf(fout,'%12.26f\n',fftc);
fclose(fout);

k=[0:nts-1];
ex=exp(1i*((xs-1)*2*pi/nts-ts).'*k);
jacobi3=jacobi2.*ex;
F2=exp(1i*(ts.')*k);
F=exp(1i*(2*pi/nts*(xs-1).')*k);
%jacobi3(1:5,27)
error4=norm(jacobi1-jacobi3)/norm(jacobi1)
error5=norm(((U1*V1').*F2)*d-((U2*V2').*F)*d)/norm((jacobi1.*F2)*d)
error6=norm(fftc-(jacobi1.*F2)*d)/norm((jacobi1.*F2)*d)
result3=(jacobi1.*F2)*d;
fout=fopen('result3.txt','w');
fprintf(fout,'%12.26f %12.26f %12.26f\n',real(result3));
fclose(fout);
fout=fopen('U1r.txt','w');
fprintf(fout,'%12.26f %12.26f %12.26f\n',real(U1));
 fclose(fout);
fout=fopen('U1i.txt','w');
fprintf(fout,'%12.26f %12.26f %12.26f\n',real(-1i*U1));
 fclose(fout);
fout=fopen('V1r.txt','w');
fprintf(fout,'%12.26f %12.26f %12.26f\n',real(V1));
 fclose(fout);
fout=fopen('V1i.txt','w');
fprintf(fout,'%12.26f %12.26f %12.26f\n',real(-1i*V1));
 fclose(fout);
