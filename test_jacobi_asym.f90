program test_jacobi_asym
use utils
use chebyshev
implicit double precision (a-h,o-z)

double precision, allocatable :: ts(:),vals(:),vals0(:),xx(:)

pi  = acos(-1.0d0)





!
!  Test jacobi_asym_p
!

nts = 500
nu  = 101
dnu = nu
da  =-0.25d0
db  = 0.40d0
p   = dnu + (da+db+1)/2

allocate(ts(nts),vals(nts),vals0(nts),xx(0:nu))

do i=1,nts
call random_number(dd)
ts(i) = dd*2/p
end do

call elapsed(t1)
call jacobi_asym_p(dnu,da,db,nts,ts,vals)
call elapsed(t2)
call prin2("jacobi_asym_p average time = ",(t2-t1)/nts)


call elapsed(t1)
do i=1,nts
call jacobi_recurrence(nu,da,db,ts(i),xx)
vals0(i) = xx(nu)
end do
call elapsed(t2)
call prin2("recurrence average time = ",(t2-t1)/nts)
call prin2("max absolute error = ",maxval(abs(vals-vals0)))

deallocate(ts,vals,vals0,xx)

!
!  Test gammaratio1
!

a  = -0.25d0
x1 = 0.5d0
x2 = 1.0d0
x3 = 10.0d0
x4 = 49.0d0
x5 = 100.0d0
x6 = 149.0d0
x7 = 151.0d0
x8 = 500.0d0
x9 = 10000.0d0

call jacobi_gammaratio1(x1,a,val1)
call jacobi_gammaratio1(x2,a,val2)
call jacobi_gammaratio1(x3,a,val3)
call jacobi_gammaratio1(x4,a,val4)
call jacobi_gammaratio1(x5,a,val5)
call jacobi_gammaratio1(x6,a,val6)
call jacobi_gammaratio1(x7,a,val7)
call jacobi_gammaratio1(x8,a,val8)
call jacobi_gammaratio1(x9,a,val9)

val10  = 2.04553134422633734322220679712928332300701686492537835536334d0
val20  = 1.22541670246517764512909830336289052685123924810807061123012d0
val30  = 0.571424713101434844736495566266296060448231324515004087052830d0
val40  = 0.379177835639598205453552908743564571342287452400680471496010d0
val50  = 0.316723497898039427344880473514962279415749253192731879094123d0
val60  = 0.286523088822615169255926432456272775772957211430684240150443d0
val70  = 0.285565608142685009504078468007768172956825545134096802602960d0
val80  = 0.211540381785998200860427276724980713984423891854144871619748d0
val90  = 0.100001562551271019004853877472726826542911993522300821671556d0

print *,"jacobi_gammaratio1: "
print *,"      x1 = ",x1,"   ",(val1-val10)/val10
print *,"      x2 = ",x2,"   ",(val2-val20)/val20
print *,"      x3 = ",x3,"   ",(val3-val30)/val30
print *,"      x4 = ",x4,"   ",(val4-val40)/val40
print *,"      x5 = ",x5,"   ",(val5-val50)/val50
print *,"      x6 = ",x6,"   ",(val6-val60)/val60
print *,"      x7 = ",x7,"   ",(val7-val70)/val70
print *,"      x8 = ",x8,"   ",(val8-val80)/val80
print *,"      x9 = ",x9,"   ",(val9-val90)/val90

!
!  Test jacobi_bessel1
!

print *,"jacobi_bessel1: "

da  = 0.25d0
t   = 100.0d0

call jacobi_bessel1(da,t,val,der,der2)
val0  = 0.00636613805035886966269762012261578584014344724315586441506066d0
der0  = -0.0000636601872328833723843206720701625153487976929692112743034492d0
der20 = 1.27316795435782353970581792878196519524223398213775235434104d-6
print *,"   t = ",t," da = ",da
print *,"    ",val,val0,(val-val0)/val0
print *,"    ",der,der0,(der-der0)/der0
print *,"    ",der2,der20,(der2-der20)/der20

da  =-0.25d0
t   = 50.0d0

call jacobi_bessel1(da,t,val,der,der2)
val0  = 0.0127319182952151261511985831212112696498750080433176290206299d0
der0  = -0.000254619292300921414959871851995903738258529076697720669832133d0
der20 = 0.0000101836282723583706653047735391662562178967614470777200613716d0
print *,"   t = ",t," da = ",da
print *,"    ",val,val0,(val-val0)/val0
print *,"    ",der,der0,(der-der0)/der0
print *,"    ",der2,der20,(der2-der20)/der20

da  = 0.40d0
t   = 25.0d0

call jacobi_bessel1(da,t,val,der,der2)
val0  = 0.0254629619913649429010996918372719478670784981999212791226435d0
der0  = -0.00101837252942591658028465447877676844501483493915034414520885d0
der20 = 0.0000814523457977998243858825598507932181695815585996981189541041d0
print *,"   t = ",t," da = ",da
print *,"    ",val,val0,(val-val0)/val0
print *,"    ",der,der0,(der-der0)/der0
print *,"    ",der2,der20,(der2-der20)/der20


da  = 0.25d0
t   = 15.0d0

call jacobi_bessel1(da,t,val,der,der2)
val0  = 0.0424237603951985580917780706906959681369947594283358406447415d0
der0  = -0.00282592611342446429518611860553464270812361136212364359717273d0
der20 = 0.000376329528509670823360421801109311507489776319221235904570554d0
print *,"   t = ",t," da = ",da
print *,"    ",val,val0,(val-val0)/val0
print *,"    ",der,der0,(der-der0)/der0
print *,"    ",der2,der20,(der2-der20)/der20



!
!  Test jacobi_bessel2
!

print *,"jacobi_bessel2: "

da  = -0.25d0
t   = 25.0d0

call jacobi_bessel2(da,t,val,der,der2)
val0  = 0.000254457814901897588786177450713642504480912897856985915855330d0
der0  = -0.0000203414881956271618532442665482664446366713799672322709089305d0
der20 = 2.43856781932228419008153737779302365324160808625178655364978d-6
print *,"   t = ",t," da = ",da
print *,"    ",val,val0,(val-val0)/val0
print *,"    ",der,der0,(der-der0)/der0
print *,"    ",der2,der20,(der2-der20)/der20


da  = 0.25d0
t   = 15.0d0

call jacobi_bessel2(da,t,val,der,der2)
val0  = 0.00212002572996554144912269381427892082301171833820075247599872d0
der0  = -0.000282401044362130436701955178081006430143609192788607390732313d0
der20 = 0.0000564092313289393141011437971472190302321758889268282662842943d0
print *,"   t = ",t," da = ",da
print *,"    ",val,val0,(val-val0)/val0
print *,"    ",der,der0,(der-der0)/der0
print *,"    ",der2,der20,(der2-der20)/der20

!
!  Test the asymptotic expansion of P
!

print *,"jacobi_asym_pder: "

dnu = 100.0d0
da  = 0.25d0
db  = 0.0d0
t   = 1.0d-7

call jacobi_asym_pder(dnu,da,db,t,valp,derp)
valp0 = 0.0000207089580363564892806314796510142113227786980129630623834471d0
derp0 = 155.317185264286153306973046398222736224709840329893166106408d0
print *,"   dnu = ",dnu," da = ",da," db = ",db
print *,"   ",valp,valp0,(valp-valp0)/valp0
print *,"   ",derp,derp0,(derp-derp0)/derp0

dnu = 0.0d0
da  = 0.25d0
db  = 0.0d0
t   = 1.0d-7

call jacobi_asym_pder(dnu,da,db,t,valp,derp)
valp0 = 5.92655664405593664348183363124978703603348367425249756713283d-6
derp0 = 44.4491748304194137031766761854934290243240672861420649530163d0
print *,"   dnu = ",dnu," da = ",da," db = ",db
print *,"   ",valp,valp0,(valp-valp0)/valp0
print *,"   ",derp,derp0,(derp-derp0)/derp0

dnu = 0.1d0
da  =-0.25d0
db  = 0.15d0
t   = 1.0d-7

call jacobi_asym_pder(dnu,da,db,t,valp,derp)
valp0 = 0.0245668570381568038664692220346083842355189943549415191724403d0
derp0 = 61417.1425953913954947471011658803184524017194189407404976079d0
print *,"   dnu = ",dnu," da = ",da," db = ",db
print *,"   ",valp,valp0,(valp-valp0)/valp0
print *,"   ",derp,derp0,(derp-derp0)/derp0

dnu = 1000.25d0
da  =-0.25d0
db  = 0.15d0
t   = 1.0d-7

call jacobi_asym_pder(dnu,da,db,t,valp,derp)
valp0 = 0.00384554290897857127160654676450845847387227033896502108553830d0
derp0 = 9613.85701571783887805614046608865397531849811853133771275580d0
print *,"   dnu = ",dnu," da = ",da," db = ",db
print *,"   ",valp,valp0,(valp-valp0)/valp0
print *,"   ",derp,derp0,(derp-derp0)/derp0

dnu = 100000.25d0
da  =-0.25d0
db  = 0.15d0
t   = 1.0d-7

call jacobi_asym_pder(dnu,da,db,t,valp,derp)
valp0 = 0.001216237542894494544758170848577482543393d0
derp0 = 3039.78300541095088724813607834909049286535280645004331335792d0
print *,"   dnu = ",dnu," da = ",da," db = ",db
print *,"   ",valp,valp0,(valp-valp0)/valp0
print *,"   ",derp,derp0,(derp-derp0)/derp0



dnu = 2.0d0**20+.25d0
da  =-0.25d0
db  = 0.15d0
t   = 1.0d-7

call jacobi_asym_pder(dnu,da,db,t,valp,derp)
valp0 = 0.000673426516860561439549119918571668208856844519990784611026123d0
derp0 = 1634.09988327939569308220489611918378193706189809704371159661d0
print *,"   dnu = ",dnu," da = ",da," db = ",db
print *,"   ",valp,valp0,(valp-valp0)/valp0
print *,"   ",derp,derp0,(derp-derp0)/derp0


print *,"jacobi_asym_amp: "



dnu = 1000.d0
da  = 0.25d0
db  = 0.0d0

call elapsed(t1)
call jacobi_asym_amp(dnu,da,db,val,der,der2)
call elapsed(t2)
print *,"   dnu = ",dnu," da = ",da," db = ",db

val0  = 0.000999374743314905889856093661488700852715267864656318397991739d0
der0  = 3.08352339795809150402122004868273980943784495053728177062803d-9
der20 = -2.36266366067927152934702865786488334771648318332586234652635d-8
print *,"    ",val,val0,(val-val0)/val0
print *,"    ",der,der0,(der-der0)
print *,"    ",der2,der20,(der2-der20)

dnu = 100.d0
da  =-0.25d0
db  = 0.20d0

call elapsed(t1)
call jacobi_asym_amp(dnu,da,db,val,der,der2)
call elapsed(t2)
print *,"   dnu = ",dnu," da = ",da," db = ",db

val0  = 0.00995209118352737438760720086786124090039859808694432676743026d0
der0  = 3.04035313195607295098848273551341643102427296099324257039045d-6
der20 = -0.0000232532073065378258375441962532519680535626515540271903847359d0
print *,"    ",val,val0,(val-val0)/val0
print *,"    ",der,der0,der-der0
print *,"    ",der2,der20,(der2-der20)

dnu = 50.d0
da  =-0.25d0
db  = 0.20d0

call elapsed(t1)
call jacobi_asym_amp(dnu,da,db,val,der,der2)
call elapsed(t2)
print *,"   dnu = ",dnu," da = ",da," db = ",db

val0  = 0.0198068068857007600111257478501493718855875227158602953766196d0
der0  = 0.0000238336410817457254095252804309440307367001106888455224273120d0
der20 = -0.000181555668678042939593732109123015102644961810837941751604479d0
print *,"    ",val,val0,(val-val0)/val0
print *,"    ",der,der0,(der-der0)
print *,"    ",der2,der20,(der2-der20)

end program
