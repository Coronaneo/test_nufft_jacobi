dprogram test_nufft
use utils
use chebyshev
implicit double precision (a-h,o-z)

type(chebexps_data)           :: chebdata
double precision, allocatable :: ab(:,:)
double precision, allocatable :: psivals(:),avals(:)
double precision, allocatable :: ts(:),avals0(:),psivals0(:),polvals(:),polvals0(:)
double precision, allocatable :: xx(:),xs(:)
complex*16,allocatable :: jacobi1(:,:),jacobi2(:,:)
integer nts

pi  = acos(-1.0d0)
nts = 2**15-16
allocate(jacobi1(nts,nts),jacobi2(nts,nts))
allocate(ts(nts),xs(nts),twhts(nts),avals0(nts),psivals0(nts),polvals(nts),polvals0(nts))

a = pi/2*(2**(-13))
b = pi - a				  

call jacobi_quad_mod(nts,da,db,ts,twhts)

do i=1,nts
xs(i) = floor(ts(i)/2/pi*nts+0.5)*2*pi/nts
end do
call quicksort(nts,ts)
call quicksort(nts,xs)
!call prin2("ts = ",ts)


do i = 27,nts

k  = i
call chebexps(k,chebdata)
if (i.ge.28) then
   deallocate(xx,psivals,avals,ab)
end if


! this must be integer because we use the recurrence relations to test accuracy
nu  = i

dnu = nu
da  = 0.25d0
db  = 0.25d0
p   = dnu + (da+db+1)/2
allocate(xx(0:nu))

dd      = nu
dd      = 1/(dd) * 2/pi
dd      = log(dd)/log(2.0d0)
nints   = ceiling(-dd)*2
nints   = max(nints,20)
allocate(ab(2,nints))

call jacobi_phase_disc(nints,ab)
!call prini("nints = ",nints)
!call prin2("before jacobi_phase, ab = ",ab)

allocate(psivals(k*nints),avals(k*nints))

call elapsed(t1)
call jacobi_phase(chebdata,dnu,da,db,nints,ab,avals,psivals)
call elapsed(t2)
!call prin2("time to construct phase = ",(t2-t1))






call elapsed(t1)
call jacobi_phase_eval(chebdata,dnu,da,db,nints,ab,avals,psivals,nts,ts,avals0,psivals0)
call elapsed(t2)
!call prin2("average eval time = ",(t2-t1)/nts)
!polvals0 = cos(psivals0)*avals0
!print *,size(avals0),size(psivals0)
jacobi1(:,i-26)=avals0*exp(dcmplx(0,1)*(psivals0-k*ts))
jacobi2(:,i-26)=avals0*exp(dcmplx(0,1)*(psivals0-k*xs))
end do

open(unit=10,file = "jacobi1r.txt")
write(10,*) real(jacobi1)
open(unit=10,file = "jacobi1i.txt")
write(10,*) real(-dcmplx(0,1)*jacobi1)

open(unit=10,file = "jacobi2r.txt")
write(10,*) real(jacobi2)
open(unit=10,file = "jacobi2i.txt")
write(10,*) real(-dcmplx(0,1)*jacobi2)


end program
