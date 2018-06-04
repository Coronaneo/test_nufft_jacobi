include make.inc

OS = $(shell uname)

.SUFFIXES: .o .a .x  .f90 .F90 .f

ALLOBJ = dfft.o dfftpack.o next235.o 
ALLOBJ1 = utils.o amos.o chebyshev.o jacobi_asym.o jacobi_taylor.o jacobi_phase.o jacobi_quad.o
ALLOBJ2 = gspiv.o orthom.o idecomp.o jacobi_exp.o jacobi_transform.o
ALLOBJ3 = nufft1df90.o nufft2df90.o nufft3df90.o

all : ${ALLOBJ} ${ALLOBJ1} ${ALLOBJ2} ${ALLOBJ3} nufft1dIInyumex.mex jacobiexample.mex chebjacex.mex nufft2dIInyumex.mex directcheb2.mex extrcheb2.mex directcheb3.mex \
	extrcheb3.mex nufft3dIInyumex.mex directjac1.mex extrjac1.mex

nufft1dIInyumex.mex: nufft1dIInyumex.F90
	${MEX} ${FLAGS} nufft1dIInyumex.F90 $(ALLOBJ) nufft1df90.o

nufft2dIInyumex.mex: nufft2dIInyumex.F90
	${MEX} ${FLAGS} nufft2dIInyumex.F90 $(ALLOBJ) nufft2df90.o

nufft3dIInyumex.mex: nufft3dIInyumex.F90
	${MEX} ${FLAGS} nufft3dIInyumex.F90 $(ALLOBJ) nufft3df90.o

directjac1.mex:directjac1.F90
	${MEX} ${FLAGS} directjac1.F90 $(ALLOBJ1) $(LIBNAME)

directcheb2.mex:directcheb2.F90
	${MEX} ${FLAGS} directcheb2.F90 $(ALLOBJ1) $(LIBNAME)

directcheb3.mex:directcheb3.F90
	${MEX} ${FLAGS} directcheb3.F90 $(ALLOBJ1) $(LIBNAME)

extrjac1.mex:extrjac1.F90
	${MEX} ${FLAGS} extrjac1.F90 $(ALLOBJ1) $(LIBNAME)

extrcheb2.mex:extrcheb2.F90
	${MEX} ${FLAGS} extrcheb2.F90 $(ALLOBJ1) $(LIBNAME)

extrcheb3.mex:extrcheb3.F90
	${MEX} ${FLAGS} extrcheb3.F90 $(ALLOBJ1) $(LIBNAME)
	
jacobiexample.mex: jacobiexample.F90
	${MEX} ${FLAGS} jacobiexample.F90 $(ALLOBJ1) $(LIBNAME)

chebjacex.mex: chebjacex.F90
	${MEX} ${FLAGS} chebjacex.F90 $(ALLOBJ1) $(ALLOBJ2) $(LIBNAME) 
#libdfftpack.a:
#	cd dfftpack && $(MAKE) clean && $(MAKE) && cp libdfftpack.a ..

LINK_MACRO = $< nufft1dIInyumex.o jacobiexample.o chebjacex.o nufft2dIInyumex.o directcheb2.o extrcheb2.o directcheb3.o extrcheb3.o nufft3dIInyumex.o directjac1.o extrjac1.o -o $@

clean : 
	rm -f *.a
	rm -f *.o
	rm -f *.x
	rm -f *.mod
	rm -f *.mexa64

.f.o : 
	$(FORTRAN) $(OPTS) $(FPPFLAGS) -c $<  -o $@ $(LIBNAME)

.F90.o : 
	$(FORTRAN) $(OPTS) $(FPPFLAGS) -c $<  -o $@ $(LIBNAME)

.f90.o : 
	$(FORTRAN) $(OPTS) $(FPPFLAGS) -c $<  -o $@ $(LIBNAME)

.f.x : 
	$(FORTRAN) $(OPTS) $(FPPFLAGS) -pg $(LINK_MACRO)

.F90.x : 
	$(FORTRAN) $(OPTS) $(FPPFLAGS) -pg $(LINK_MACRO)

.f90.x : 
	$(FORTRAN) $(OPTS) $(LINK_MACRO)
