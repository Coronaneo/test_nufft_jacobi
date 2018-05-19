program test_nufft2

implicit double precision (a-h,o-z)

integer nts,ier,iflag,i
parameter (nts=2**9,iflag=1)
integer :: rank(2)
double precision eps,r13,r2,r23,r3
parameter (eps=1e-12)
double precision :: result2(nts),result3(nts),c(nts),ts(nts)
double precision,allocatable :: U1r(:,:),U1i(:,:),V1r(:,:),V1i(:,:)
complex*16,allocatable :: U1(:,:),V1(:,:),uu(:,:)
complex*16 :: result1(nts),ec(nts),cj(nts)

open(unit = 10,file = 'c.txt')
read(10,*) c
open(unit = 10,file = 'ts.txt')
read(10,*) ts
open(unit = 10,file = 'result2.txt')
read(10,*) result2
open(unit = 10,file = 'result3.txt')
read(10,*) result3
open(unit = 10,file = 'rank.txt')
read(10,*) rank
allocate(U1r(nts,rank(1)),U1i(nts,rank(1)),uu(nts/2,rank(1)))
allocate(V1r(nts,rank(1)),V1i(nts,rank(1)))
allocate(U1(nts,rank(1)),V1(nts,rank(1)))
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
!uu=U1(1:nts/2,:)
!U1(1:nts/2,:)=U1(nts/2+1:nts,:)
!U1(nts/2+1:nts,:)=uu
!re=ts(1:nts/2)
!ts(1:nts/2)=ts(nts/2+1:nts)
!ts(nts/2+1:nts)=re

print *,'c =',c(1:5)
print *,'rank1 = ',rank(1),'rank2 = ',rank(2)
result1=0

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
