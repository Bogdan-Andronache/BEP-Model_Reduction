TFEMPATH = /opt/tfem

#
# Mdefs.mk
#
# Some general macros

RM = rm -f
AR = ar vr
RANLIB = ranlib
OUTCLEAN = *.fig *.vtk *.out *.ps *.plt *.pos *.png *.mso

# Some tfem macros

EXTRATAGSSRC =
MODPATHINC = $(TFEMPATH)/modfiles

# the following lines are for building a library only

LIBCREATE = libdummy.a
MODPATHWRITE = .

# the following is for building a project tree including a library

ifdef PROJECTDIR
  LIBCREATE = $(PROJECTDIR)/lib/lib$(PROJECTNAME).a
  MODPATHWRITE = $(PROJECTDIR)/modfiles
  LIBSP = -L$(PROJECTDIR)/lib -l$(PROJECTNAME)
  FMODINCP = -I$(MODPATHWRITE)
endif

FMODINC = $(FMODINCP) -I$(MODPATHINC)

# Blas/Lapack library

LIBBLAS = -lopenblasp

# Metis

LIBMETIS = -lmetis

# UMFPACK

LIBUMFPACK =

# Compiler stuff

FMISC = -std=f2018 -fall-intrinsics -ffpe-trap=zero,overflow,invalid
FPPFLAGS = -DMETIS5 -DSKIT2
FC = gfortran
FMODWRITE = -J$(MODPATHWRITE)
FOPT = -O3 -march=haswell -funroll-loops -fstack-arrays
FSOLVERLIB = external_gfortran

FFLAGS = $(FCHECK) $(FOPT) $(FDEBUG) $(FMISC)
LIBS = $(LIBSP) -L$(TFEMPATH)/lib -ltfem -l$(FSOLVERLIB) $(LIBUMFPACK) $(LIBBLAS) $(LIBMETIS)
