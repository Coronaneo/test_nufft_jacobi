include make.inc

OS = $(shell uname)

.SUFFIXES: .o .a .x .f90 .F90 .f

ALLOBJ = dfft.o dfftpack.o dirft1d.o dirft2d.o dirft3d.o next235.o nufft1df90.o 
ALLOBJ1 = utils.o amos.o chebyshev.o jacobi_asym.o jacobi_taylor.o jacobi_phase.o jacobi_quad.o

all : ${ALLOBJ} ${ALLOBJ1} nufft1dIInyumex.mex jacobiexample.mex

nufft1dIInyumex.mex: nufft1dIInyumex.f90
	${MEX} ${FLAGS} nufft1dIInyumex.f90 $(ALLOBJ)

jacobiexample.mex: jacobiexample.f90
	${MEX} ${FLAGS} jacobiexample.f90 $(ALLOBJ1) $(LIBNAME)

LINK_MACRO = $< nufft1dIInyumex.o -o $@

clean : 
	rm -f *.a
	rm -f *.o
	rm -f *.x
	rm -f *.mod

.f.o : 
	$(FORTRAN) $(OPTS) $(FPPFLAGS) -c $<  -o $@

.F90.o : 
	$(FORTRAN) $(OPTS) $(FPPFLAGS) -c $<  -o $@

.f90.o : 
	$(FORTRAN) $(OPTS) $(FPPFLAGS) -c $<  -o $@

.f.x : 
	$(FORTRAN) $(OPTS) $(FPPFLAGS) -pg $(LINK_MACRO)

.F90.x : 
	$(FORTRAN) $(OPTS) $(FPPFLAGS) -pg $(LINK_MACRO)

.f90.x : 
	$(FORTRAN) $(OPTS) $(LINK_MACRO)
