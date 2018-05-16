#-fdefault-real-8
#-freal-8-real-10


# Choose the compiler and set compiler options

CPPCOMP  = g++
CPPOPTS  = -O3 

# FCOMP    = ifort
# FOPTS    = -O3
# LDOPT =  

FCOMP    = gfortran
FOPTS    = -w -Ofast -cpp 
LDOPT    = -lfftw3 -lblas

FC     = $(FCOMP)
FFLAGS = $(FOPTS)
export FC FFLAGS

# Set the list of programs to compile

PROGRAMS = test_chebyshev test_idecomp test_jacobi_asym test_jacobi_taylor \
           test_jacobi_phase test_jacobi_quad test_nufft test_jacobi_exp   \
           test_jacobi_transform

# Compile all of the test programs and the library

all	        : clean $(PROGRAMS) 



# List the dependencies for each module's test program

JACOBI_TRANSFORM_FILES = utils.o                       \
                         gspiv.o                       \
                         orthom.o                      \
                         amos.o                        \
	                 idecomp.o                     \
                         chebyshev.o                   \
                         jacobi_asym.o                 \
                         jacobi_taylor.o               \
                         jacobi_phase.o                \
                         jacobi_quad.o                 \
                         jacobi_exp.o                  \
                         jacobi_transform.o            \
                         libdfftpack.a


JACOBI_EXP_FILES       = utils.o                       \
                         gspiv.o                       \
                         orthom.o                      \
                         amos.o                        \
	                 idecomp.o                     \
                         chebyshev.o                   \
                         jacobi_asym.o                 \
                         jacobi_taylor.o               \
                         jacobi_phase.o                \
                         jacobi_quad.o                 \
                         jacobi_exp.o


JACOBI_QUAD_FILES      = utils.o                       \
                         amos.o                        \
                         chebyshev.o                   \
                         jacobi_asym.o                 \
                         jacobi_taylor.o               \
                         jacobi_phase.o                \
                         jacobi_quad.o


JACOBI_PHASE_FILES     = utils.o                       \
                         amos.o                        \
                         chebyshev.o                   \
                         jacobi_asym.o                 \
                         jacobi_taylor.o               \
                         jacobi_phase.o

JACOBI_TAYLOR_FILES    = utils.o                       \
                         chebyshev.o                   \
                         jacobi_asym.o                 \
                         jacobi_taylor.o

JACOBI_ASYM_FILES      = utils.o                       \
                         chebyshev.o                   \
                         jacobi_asym.o


NUFFT_FILES            = utils.o                       \
                         gspiv.o                       \
                         nrleastsq2.o                  \
                         chebyshev.o                   \
                         jacobi_asym.o                 \
                         jacobi_taylor.o               \
                         jacobi_phase.o                \
                         jacobi_quad.o                 \
                         nufft.o

IDECOMP_FILES          = utils.o                       \
                         orthom.o                      \
                         gspiv.o                       \
                         idecomp.o


CHEBYSHEV_FILES        = utils.o                       \
                         chebyshev.o





###################################################################################

libdfftpack.a   :
	cd dfftpack && $(MAKE) clean && $(MAKE) && cp libdfftpack.a ..

test_jacobi_transform.o : $(JACOBI_TRANSFORM_FILES) test_jacobi_transform.f90
test_jacobi_transform   : $(JACOBI_TRANSFORM_FILES) test_jacobi_transform.o

test_jacobi_exp.o       : $(JACOBI_EXP_FILES) test_jacobi_exp.f90
test_jacobi_exp         : $(JACOBI_EXP_FILES) test_jacobi_exp.o

test_jacobi_quad.o      : $(JACOBI_QUAD_FILES) test_jacobi_quad.f90
test_jacobi_quad        : $(JACOBI_QUAD_FILES) test_jacobi_quad.o

test_jacobi_phase.o     : $(JACOBI_PHASE_FILES) test_jacobi_phase.f90
test_jacobi_phase       : $(JACOBI_PHASE_FILES) test_jacobi_phase.o

test_jacobi_taylor.o    : $(JACOBI_TAYLOR_FILES) test_jacobi_taylor.f90
test_jacobi_taylor      : $(JACOBI_TAYLOR_FILES) test_jacobi_taylor.o

test_jacobi_asym.o      : $(JACOBI_ASYM_FILES) test_jacobi_asym.f90
test_jacobi_asym        : $(JACOBI_ASYM_FILES) test_jacobi_asym.o

test_nufft.o            : $(NUFFT_FILES) test_nufft.f90
test_nufft              : $(NUFFT_FILES) test_nufft.o

test_idecomp.o          : $(IDECOMP_FILES) test_idecomp.f90
test_idecomp            : $(IDECOMP_FILES) test_idecomp.o

test_chebyshev.o        : $(CHEBYSHEV_FILES) test_chebyshev.f90
test_chebyshev          : $(CHEBYSHEV_FILES) test_chebyshev.o


# Setup the general compilation rules

%		: %.o
	$(FCOMP) $(FOPTS) -o $@ $^ $(LDOPT)
	@echo  
	@echo 
	@echo "---------[ $@     ]--------------------------------------------"
	@echo 
	@./$@
	@echo 
	@echo "--------------------------------------------------------------------------"
	@echo 

%.o		: %.f90
	$(FCOMP) -c $(FOPTS)  $<

%.o		: %.f
	$(FCOMP) -c $(FOPTS)  $<

%.o		: %.cpp
	$(CPPCOMP) -c $(CPPOPTS)  $<

.PHONY: clean

clean:
	rm -f *.o *.mod *~ fort.* a.out *.a
	rm -f $(PROGRAMS)
	rm -f gn??? gn???.dat
