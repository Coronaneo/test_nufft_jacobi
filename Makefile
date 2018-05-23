include make.inc

OS = $(shell uname)

.SUFFIXES: .o .a .x  .f90 .F90 .f

ALLOBJ = dfft.o dfftpack.o dirft1d.o next235.o nufft1df90.o 
ALLOBJ1 = utils.o amos.o gspiv.o orthom.o idecomp.o chebyshev.o jacobi_asym.o jacobi_taylor.o jacobi_phase.o jacobi_quad.o jacobi_exp.o jacobi_transform.o


all : ${ALLOBJ} ${ALLOBJ1} nufft1dIInyumex.mex jacobiexample.mex

nufft1dIInyumex.mex: nufft1dIInyumex.F90
	${MEX} ${FLAGS} nufft1dIInyumex.F90 $(ALLOBJ)

jacobiexample.mex: jacobiexample.F90
	${MEX} ${FLAGS} jacobiexample.F90 $(ALLOBJ1) $(LIBNAME)

chebjacex.mex: chebjacex.F90
	${MEX} ${FLAGS} chebjacex.F90 $(ALLOBJ1) $(LIBNAME) 

LINK_MACRO = $< nufft1dIInyumex.o jacobiexample.o chebjacex.o -o $@

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
