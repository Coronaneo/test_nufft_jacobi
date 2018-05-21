# ----------------------------machine dependent file

# define MEX, in mac mex may also refer to pdftex, a command conflict
MEX = /usr/local/MATLAB/R2017a/bin/mex

# Tested the flag on R2017a
FLAGS = 

FORTRAN   = gfortran
OPTS      = -O3
FPPFLAGS  = -DLAP -DSLAP -fPIC

LIBNAME = 
