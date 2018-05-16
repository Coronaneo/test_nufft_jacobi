nts=2^9-26;
load jacobi1i.txt
load jacobi1r.txt
load jacobi2r.txt
load jacobi2i.txt
jacobi1r=reshape(jacobi1r,[nts,nts]);
jacobi1i=reshape(jacobi1i,[nts,nts]);
jacobi2r=reshape(jacobi2r,[nts,nts]);
jacobi2i=reshape(jacobi2i,[nts,nts]);
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
