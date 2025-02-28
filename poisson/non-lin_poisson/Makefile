include Mdefs.mk

.SUFFIXES:

.SUFFIXES: .o .f90 .F90 .f

.f90.o:
	$(FC) -c $(FFLAGS) $(FMODINC) $<

SRC = *.f90

OBJS = poisson_functions.o

MAINP = poisson_nonlin1

all: objects programs tags

objects: $(OBJS)

programs: $(MAINP)

tags: $(SRC)
	ctags --fortran-kinds=+i $(SRC) $(EXTRATAGSSRC)

poisson_nonlin1.o: poisson_functions.o
poisson_nonlin1: $(OBJS) poisson_nonlin1.o
	$(FC) $(FFLAGS) -o $@ poisson_nonlin1.o $(OBJS) $(LIBS)

clean:
	$(RM)  poisson_nonlin1 *.o *.mod *.bin tags

outclean:
	$(RM) $(OUTCLEAN)

allclean: clean outclean
