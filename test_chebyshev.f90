module test_chebyshev_subroutines

use utils
use chebyshev

contains


subroutine testfun(ts,vals)
implicit double precision (a-h,o-z)
double precision, intent(in)     :: ts(:)
double precision, intent(out)    :: vals(:)
vals = cos(ts) + sin(ts)**2 + ts**2 - 1
end subroutine


end module

program test_chebyshev

use utils
use chebyshev
use test_chebyshev_subroutines
use iso_c_binding

implicit double precision (a-h,o-z)

type(chebexps_data)           :: chebdata
double precision, allocatable :: ab(:,:),vals0(:,:),ts(:),vals(:)
double precision, allocatable :: z(:),w(:,:)

!
!  Build a piecewise discretization scheme
!

k = 16
call chebexps(k,chebdata)

! data xs / -0.1000000000000000D+01, &
!           -0.9781476007338056D+00, &
!           -0.9135454576426008D+00, &
!           -0.8090169943749473D+00, &
!           -0.6691306063588579D+00, &
!           -0.4999999999999998D+00, &
!           -0.3090169943749473D+00, &
!           -0.1045284632676533D+00, &
!            0.1045284632676537D+00, &
!            0.3090169943749475D+00, &
!            0.5000000000000001D+00, &
!            0.6691306063588582D+00, &
!            0.8090169943749475D+00, &
!            0.9135454576426009D+00, &
!            0.9781476007338057D+00, &
!            0.1000000000000000D+01  /

do i=1,16
write (*,"('      ',D44.36,', &')") chebdata%xs(i)
end do
stop

a = 0.0d0
b = 1.0d0

nints = 30
allocate(ab(2,nints))
do i=1,nints
ab(1,i) = (i-1.0d0)/(nints+0.0d0)
ab(2,i) = (i+0.0d0)/(nints+.0d0)
end do


!
!  Sample a test function on this grid
!

allocate(vals0(k,nints))

do int=1,nints
a = ab(1,int)
b = ab(2,int)
do i=1,k
x = chebdata%xs(i) * (b-a)/2 + (b+a)/2
vals0(i,int) = cos(133*x)
end do
end do

m = 100
allocate(ts(m),vals(m))
do i=1,m
ts(i) = (i-1.0d0)/(m-1.0d0)
end do

!
!  Call the interpolation routine
!

call elapsed(t1)
call chebpw_aint(chebdata,nints,ab,m,ts,vals0,vals)
call elapsed(t2)
call prin2("aint average time = ",(t2-t1)/m)

errmax = 0
do i=1,m
t = ts(i)
errabs = abs(vals(i)-cos(133*t))
errmax = max(errabs,errmaX)
end do
call prin2("aint errmax = ",errmax)

!
!
!

allocate(z(m),w(k,nints))

do i=1,m
call random_number(z(i))
end do

dsum1 = 0.0d0
do i=1,m
t = ts(i)
dsum1 = dsum1 + z(i)*cos(133*t)
end do

call elapsed(t1)
call chebpw_aintt(chebdata,nints,ab,m,ts,z,w)
call elapsed(t2)
call prin2("aintl average time = ",(t2-t1)/m)

dsum2 = 0.0d0
do int=1,nints
do i=1,k
t = ts(i)
dsum2 = dsum2 + w(i,int)*vals0(i,int)
end do
end do

call prin2("aintl error = ",dsum1-dsum2)
stop


! allocate(vals(k),coefs(k))

! call testfun(chebdata%xs,vals)
! coefs = matmul(chebdata%u, vals)

! nts = 100
! allocate(vals0(nts),vals1(nts),ts(nts))

! do i=1,nts
! call random_number(dd)
! ts(i) = -1 + dd*2
! end do


! a = -1.0d0
! b =  1.0d0

! call elapsed(t1)
! do i=1,nts
! call chebeval(a,b,k,chebdata%xs,vals,ts(i),vals1(i))
! end do
! call elapsed(t2)
! call prin2("average time for chebeval = ",(t2-t1)/nts)


! call elapsed(t1)
! do i=1,nts
! call chebeval2(a,b,k,coefs,ts(i),vals1(i))
! end do
! call elapsed(t2)
! call prin2("average time for chebeval2 = ",(t2-t1)/nts)


! call elapsed(t1)
! call testfun(ts,vals0)
! call elapsed(t2)
! call prin2("average time for testfun  = ",(t2-t1)/nts)

! call prin2("maximum error =",maxval(abs(vals1-vals0)))

end program
