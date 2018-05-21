program test_nufft2

implicit double precision (a-h,o-z)

integer nts,ier,iflag,i,j,num
parameter (nts=2**14,iflag=1,num=100)
integer :: rank(2),xs(nts)
double precision eps,r13,r2,r23,r3,pi
parameter (eps=1e-9)
double precision :: result3(nts),c(nts),ts(nts),arr(4),time1,time2
double precision,allocatable :: U1r(:,:),U1i(:,:),V1r(:,:),V1i(:,:)
complex*16,allocatable :: U1(:,:),V1(:,:)
double precision,allocatable :: U2r(:,:),U2i(:,:),V2r(:,:),V2i(:,:)
complex*16,allocatable :: U2(:,:),V2(:,:)
complex*16 :: result1(nts),ec(nts),cj(nts),result2(nts),cc(nts)

double complex in1, out1,in2,out2
dimension in1(nts), out1(nts),in2(12),out2(12)
integer*8 :: plan,plan2
integer FFTW_FORWARD,FFTW_MEASURE,FFTW_BACKWARD
!parameter (FFTW_FORWARD=-1)
parameter (FFTW_MEASURE=1)

character*8 date
character*10 time
character*5 zone
integer*4 values1(8),values2(8)
arr(1)=3600
arr(2)=60
arr(3)=1
arr(4)=0.001

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
allocate(U2r(nts,rank(2)),U2i(nts,rank(2)))
allocate(V2r(nts,rank(2)),V2i(nts,rank(2)))
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
!print *,'V1=',V1(5:10,1)
ec=exp(dcmplx(0,1)*nts/2*ts)

open(unit = 10,file = 'U2r.txt')
read(10,*) U2r
open(unit = 10,file = 'U2i.txt')
read(10,*) U2i
open(unit = 10,file = 'V2r.txt')
read(10,*) V2r
open(unit = 10,file = 'V2i.txt')
read(10,*) V2i
U2=transpose(U2r+dcmplx(0,1)*U2i)
V2=transpose(V2r+dcmplx(0,1)*V2i)
!print *,'V2=',V2(1,1:5)

call dfftw_plan_dft_1d(plan2,12,in2,out2,1,0)
do i=1,12
   in2(i)=dcmplx(i,1)
end do
call dfftw_execute_dft(plan2,in2,out2)
!print *,'out2 = ',out2
pi=acos(-1.0d0)
!print *,'c =',c(1:5)
print *,'rank1 = ',rank(1),'rank2 = ',rank(2)
result1=0
cc=c
call dfftw_plan_dft_1d(plan,nts,in1,out1,1,0)
xs=mod(floor(ts/2/pi*nts+0.5),nts)+1
!print *,'xs(1:5)=',xs(1:5)
call date_and_time(date,time,zone,values1)
do i=1,num
call nufft1dIIapp(nts,plan,cc,U2,V2,xs,nts,1,rank(2),result2)
end do
call date_and_time(date,time,zone,values2)
time1=sum((values2(5:8)-values1(5:8))*arr)
print *,' T_our         = ',time1/num
!result2=result2
call date_and_time(date,time,zone,values1)
do j=1,num
result1=0
do i=1,rank(1)
   call nufft1d2f90(nts,ts,cj,iflag,eps,nts,conjg(V1(:,i))*c,ier)
   result1=result1+U1(:,i)*ec*cj
end do
end do
call date_and_time(date,time,zone,values2)
time2=sum((values2(5:8)-values1(5:8))*arr)
print *,' T_nyu         = ',time2/num
print *,' T_our/T_nyu   = ',time1/time2
!result1=result1/nts
!re=result1(1:nts/2)
!result1(1:nts/2)=result1(nts/2+1:nts)
!result1(nts/2+1:nts)=re

r13=sqrt(sum((result3-real(result1))*(result3-real(result1))))
r2=sqrt(sum(result2*result2))
r3=sqrt(sum(result3*result3))
r23=sqrt(sum((real(result2)-result3)*(real(result2)-result3)))
print *,'relative error1 = ',r13/r3
print *,'relative error2 = ',r23/r3
!print *,'result1=',real(result1(1:5))
!print *,'result2=',real(result2(1:5))
!print *,'result3=',result3(1:5)
end program
