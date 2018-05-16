program test_jacobi_taylor

implicit double precision (a-h,o-z)


print *,"jacobi_taylor: "

dnu = 10.0d0
da  =-0.25d0
db  =-0.15d0
t   = 0.1d0

call jacobi_taylor_p(dnu,da,db,t,val,der,der2)
val0  = 0.256583326858048198871881934479075213d0
der0  = -1.66997344035118243165803880634212160d0
der20 = -32.0505033653418676698473232075796728d0
print *,"   dnu = ",dnu,"  da = ",da," db = ",db," t = ",t
print *,"    ",val,val0,(val-val0)/val0
print *,"    ",der,der0,(der-der0)/der0
print *,"    ",der2,der20,(der2-der20)/der20

dnu = 100.0d0
da  =-0.25d0
db  =-0.15d0
t   = 0.1d0

call jacobi_taylor_p(dnu,da,db,t,val,der,der2)
val0  = -0.0977497693238507529137258198749327592d0
der0  = 2.02169338473431363597301660272132816d0
der20 = 985.211386596241968307216946114166767d0
print *,"   dnu = ",dnu,"  da = ",da," db = ",db," t = ",t
print *,"    ",val,val0,(val-val0)/val0
print *,"    ",der,der0,(der-der0)/der0
print *,"    ",der2,der20,(der2-der20)/der20


dnu = 10.0d0
da  =-0.25d0
db  =-0.15d0
t   = 0.01d0

call jacobi_taylor_p(dnu,da,db,t,val,der,der2)
val0  = 0.213918763010064588679874422647004379d0
der0  = 5.19626112534518958709463504897678588d0
der20 = -423.807831642810366005718759839602619d0
print *,"   dnu = ",dnu,"  da = ",da," db = ",db," t = ",t
print *,"    ",val,val0,(val-val0)/val0
print *,"    ",der,der0,(der-der0)/der0
print *,"    ",der2,der20,(der2-der20)/der20

dnu = 50.0d0
da  =-0.25d0
db  =-0.15d0
t   = 0.01d0

call jacobi_taylor_p(dnu,da,db,t,val,der,der2)
val0  = 0.132462133478620672792711037600573095d0
der0  = 0.962488060337825102856569277835623720d0
der20 = -583.517223268710719916256768063740823d0
print *,"   dnu = ",dnu,"  da = ",da," db = ",db," t = ",t
print *,"    ",val,val0,(val-val0)/val0
print *,"    ",der,der0,(der-der0)/der0
print *,"    ",der2,der20,(der2-der20)/der20

dnu = 100.0d0
da  =-0.25d0
db  =-0.15d0
t   = 0.01d0

call jacobi_taylor_p(dnu,da,db,t,val,der,der2)
val0  = 0.0836217008263496146754034369345044098d0
der0  = -4.94158479419142854766686688974569961d0
der20 = -998.038588014301491815497938914199358d0
print *,"   dnu = ",dnu,"  da = ",da," db = ",db," t = ",t
print *,"    ",val,val0,(val-val0)/val0
print *,"    ",der,der0,(der-der0)/der0
print *,"    ",der2,der20,(der2-der20)/der20

dnu = 1000.0d0
da  =-0.25d0
db  =-0.15d0
t   = 0.01d0

call jacobi_taylor_p(dnu,da,db,t,val,der,der2)
val0  = -0.0311141414292331294328210093835700793d0
der0  = 5.54364549039345777737219946297963288d0
der20 = 31191.1539853651356376914353932013698d0
print *,"   dnu = ",dnu,"  da = ",da," db = ",db," t = ",t
print *,"    ",val,val0,(val-val0)/val0
print *,"    ",der,der0,(der-der0)/der0
print *,"    ",der2,der20,(der2-der20)/der20

dnu = 1000.0d0
da  =-0.25d0
db  =-0.15d0
t   = 0.001d0

call jacobi_taylor_p(dnu,da,db,t,val,der,der2)
val0  = 0.0265212782100443542561085936796288899d0
der0  = -15.5205279472275545631956717524305823d0
der20 = -31509.9349510618042431270095589975202d0
print *,"   dnu = ",dnu,"  da = ",da," db = ",db," t = ",t
print *,"    ",val,val0,(val-val0)/val0
print *,"    ",der,der0,(der-der0)/der0
print *,"    ",der2,der20,(der2-der20)/der20

end program
