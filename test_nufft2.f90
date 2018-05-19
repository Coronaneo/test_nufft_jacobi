program test_nufft2

implicit double precision (a-h,o-z)

integer nts,ier,iflag,i,xs(nts)
parameter (nts=2**9,iflag=1)
integer :: rank(2)
double precision eps,r13,r2,r23,r3
parameter (eps=1e-12)
double precision :: result3(nts),c(nts),ts(nts)
double precision,allocatable :: U1r(:,:),U1i(:,:),V1r(:,:),V1i(:,:)
complex*16,allocatable :: U1(:,:),V1(:,:)
double precision,allocatable :: U2r(:,:),U2i(:,:),V2r(:,:),V2i(:,:)
complex*16,allocatable :: U2(:,:),V2(:,:)
complex*16 :: result1(nts),ec(nts),cj(nts),result2(nts)

double complex in1, out1
dimension in1(nts), out1(nts)
integer*8 :: plan
integer FFTW_FORWARD,FFTW_MEASURE
parameter (FFTW_FORWARD=-1)
parameter (FFTW_MEASURE=0)


open(unit = 10,file = 'c.txt')
read(10,*) c
open(unit = 10,file = 'ts.txt')
read(10,*) ts
open(unit = 10,file = 'result3.txt')
read(10,*) result3
open(unit = 10,file = 'rank.txt')
read(10,*) rank
allocate(U1r(nts,rank(1)),U1i(nts,rank(1)))
allocate(V1r(nts,rank(1)),V1i(nts,rank(1)))
allocate(U1(nts,rank(1)),V1(nts,rank(1)))
allocate(U2r(rank(2),nts),U2i(rank(2),nts))
allocate(V2r(rank(2),nts),V2i(rank(2),nts))
allocate(U2(rank(2),nts),V2(rank(2),nts))
open(unit = 10,file = 'U1r.txt')
read(10,*) U1r
open(unit = 10,file = 'U1i.txt')
read(10,*) U1i
open(unit = 10,file = 'V1r.txt')
read(10,*) V1r
open(unit = 10,file = 'V1i.txt')
read(10,*) V1i
U1=U1r+dcmplx(0,1)*U1i
V1=V1r+dcmplx(0,1)*V1i
print *,'V1=',V1(5:10,1)
ec=exp(dcmplx(0,1)*nts/2*ts)

open(unit = 10,file = 'U2r.txt')
read(10,*) U2r
open(unit = 10,file = 'U2i.txt')
read(10,*) U2i
open(unit = 10,file = 'V2r.txt')
read(10,*) V2r
open(unit = 10,file = 'V2i.txt')
read(10,*) V2i
U2=U2r+dcmplx(0,1)*U2i
V2=V2r+dcmplx(0,1)*V2i

print *,'c =',c(1:5)
print *,'rank1 = ',rank(1),'rank2 = ',rank(2)
result1=0

call dfftw_plan_dft_1d(plan,nts,in1,out1,FFTW_BACKWARD,0)
xs=mod(floor(ts/2/pi*nts+0.5),nts)+1
call nufft1dIIapp(nts,plan,c,U2,V2,xs,nts,1,rank(2),result2)

do i=1,rank(1)
   call nufft1d2f90(nts,ts,cj,iflag,eps,nts,conjg(V1(:,i))*c,ier)
   result1=result1+U1(:,i)*ec*cj
end do
!result1=result1/nts
!re=result1(1:nts/2)
!result1(1:nts/2)=result1(nts/2+1:nts)
!result1(nts/2+1:nts)=re

r13=sqrt(sum((result3-real(result1))*(result3-real(result1))))
r2=sqrt(sum(result2*result2))
r3=sqrt(sum(result3*result3))
r23=sqrt(sum((result2-result3)*(result2-result3)))
print *,'relative error1 = ',r13/r3
print *,'relative error2 = ',r23/r3
print *,'result1=',real(result1(1:5))
print *,'result2=',result2(1:5)

end program
