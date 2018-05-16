program test_jacobi_solve
use utils
use chebyshev
implicit double precision (a-h,o-z)

type(chebexps_data)           :: chebdata
double precision, allocatable :: ab(:,:)
double precision, allocatable :: psivals(:),avals(:)
double precision, allocatable :: ts(:),avals0(:),psivals0(:),polvals(:),polvals0(:)
double precision, allocatable :: xx(:)


k  = 20
call chebexps(k,chebdata)

pi  = acos(-1.0d0)

! this must be integer because we use the recurrence relations to test accuracy
nu  = 2**20

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
call prini("nints = ",nints)
call prin2("before jacobi_phase, ab = ",ab)

allocate(psivals(k*nints),avals(k*nints))

call elapsed(t1)
call jacobi_phase(chebdata,dnu,da,db,nints,ab,avals,psivals)
call elapsed(t2)
call prin2("time to construct phase = ",(t2-t1))


nts = 100
allocate(ts(nts),avals0(nts),psivals0(nts),polvals(nts),polvals0(nts))

a = ab(1,1)
b = ab(2,nints)

do i=1,nts
call random_number(dd)
ts(i) = a + (b-a)*dd
end do
ts(1) = 1.1d0/dnu
ts(2) = pi-1.1d0/dnu

call quicksort(nts,ts)
call prin2("ts = ",ts)

call elapsed(t1)
call jacobi_phase_eval(chebdata,dnu,da,db,nints,ab,avals,psivals,nts,ts,avals0,psivals0)
call elapsed(t2)
call prin2("average eval time = ",(t2-t1)/nts)
polvals0 = cos(psivals0)*avals0


dconst = sqrt( (2*dnu+da+db+1)/pi)
call elapsed(t1)
do i=1,nts
t  = ts(i)

if (dnu .gt. 2**14) then
if (t .gt. pi/2) then
call jacobi_asym2(dnu,da,db,pi-t,polvals(i),valq)
else
call jacobi_asym2(dnu,da,db,t,polvals(i),valq)
endif
polvals(i) = polvals(i) / dconst
else
call jacobi_recurrence(nu,da,db,t,xx)
polvals(i) = xx(nu)
endif


end do

call prin2("errs = ",abs(polvals-polvals0))

call elapsed(t2)
call prin2("asym expansion avg eval time = ",(t2-t1)/nts)
call prin2("max error = ",maxval(abs(polvals-polvals0)))

end program


