nts=2^9;
load jacobi1i.bin
load jacobi1r.bin
load jacobi2r.bin
load jacobi2i.bin
load xs.bin
jacobi1r=reshape(jacobi1r,[nts,nts-27]);
jacobi1i=reshape(jacobi1i,[nts,nts-27]);
jacobi2r=reshape(jacobi2r,[nts,nts-27]);
jacobi2i=reshape(jacobi2i,[nts,nts-27]);
jacobi1=jacobi1r+1i*jacobi1i;
jacobi2=jacobi2r+1i*jacobi2i;
size(jacobi1);
size(jacobi2);
tol=1e-12;
tR=60;
mR=60;
[U1,V1]=lowrank(jacobi1,tol,tR,mR);
[U2,V2]=lowrank(jacobi2,tol,tR,mR);
rank1=size(U1,2)
rank2=size(U2,2)
fout=fopen('rank.txt','w')£»
fprintf(fout,'%8i\n',rank1)
fprintf(fout,'%8i\n',rank2)
fclose£¨fout£©

xs=round(xs/2/pi*nts);
c=rand(nts,1);
fout=fopen('c.txt','w')£»
fprintf(fout,'%8.16f\n',c)
fclose£¨fout£©
ncol = size(c,2);
c = repmat(conj(V2),1,ncol).*reshape(repmat(c,rank2,1),nts,rank2*ncol);
if iflag < 0
   fftc = fft(c);
else
   fftc = ifft(c);
end
fftc = fftc(xs,:);
fftc = squeeze(sum(reshape(repmat(U2,1,ncol).*fftc,nts,rank2,ncol),2));
fout=fopen('result2.txt','w')£»
fprintf(fout,'%8.16f\n',fftc)
fclose£¨fout£©

fout=fopen('U1r.txt','w')£»
fprintf(fout,'%8.16f %8.16f %8.16f\n',U1.')
 fclose£¨fout£©
fout=fopen('U1i.txt','w')£»
fprintf(fout,'%8.16f %8.16f %8.16f\n',-1i*U1.')
 fclose£¨fout£©
fout=fopen('V1r.txt','w')£»
fprintf(fout,'%8.16f %8.16f %8.16f\n',V1.')
 fclose£¨fout£©
fout=fopen('V1i.txt','w')£»
fprintf(fout,'%8.16f %8.16f %8.16f\n',-1i*V1.')
 fclose£¨fout£©
