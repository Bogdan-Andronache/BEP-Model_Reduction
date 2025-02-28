include Mdefs.mk

.SUFFIXES:

.SUFFIXES: .o .f90 .F90 .f

.f90.o:
	$(FC) -c $(FFLAGS) $(FMODINC) $<


.F90.o:
	$(FC) -c $(FPPFLAGS) $(FFLAGS) $(FMODINC) $<

SRC = *.f90 *.F90

OBJS = functions.o poisson_functions.o

MAINP = poisson1 poisson12 poisson13 poisson1a poisson1b poisson1c poisson2 poisson20 poisson21 poisson22 poisson23 poisson27 poisson28 poisson29 poisson3 poisson30 poisson31 poisson32 poisson33 poisson34 poisson35 poisson36 poisson37 poisson4 poisson5 poisson6 poisson7 poisson8 poisson9 poisson14 poisson15 poisson16 poisson17 poisson18 poisson19 poisson24 poisson25 poisson26

all: objects programs tags

objects: $(OBJS)

programs: $(MAINP)

tags: $(SRC)
	ctags --fortran-kinds=+i $(SRC) $(EXTRATAGSSRC)

poisson1.o: poisson_functions.o
poisson1: $(OBJS) poisson1.o
	$(FC) $(FFLAGS) -o $@ poisson1.o $(OBJS) $(LIBS)

poisson12.o: poisson_functions.o
poisson12: $(OBJS) poisson12.o
	$(FC) $(FFLAGS) -o $@ poisson12.o $(OBJS) $(LIBS)

poisson13.o: poisson_functions.o
poisson13: $(OBJS) poisson13.o
	$(FC) $(FFLAGS) -o $@ poisson13.o $(OBJS) $(LIBS)

poisson1a.o: poisson_functions.o
poisson1a: $(OBJS) poisson1a.o
	$(FC) $(FFLAGS) -o $@ poisson1a.o $(OBJS) $(LIBS)

poisson1b.o: poisson_functions.o
poisson1b: $(OBJS) poisson1b.o
	$(FC) $(FFLAGS) -o $@ poisson1b.o $(OBJS) $(LIBS)

poisson1c.o: poisson_functions.o
poisson1c: $(OBJS) poisson1c.o
	$(FC) $(FFLAGS) -o $@ poisson1c.o $(OBJS) $(LIBS)

poisson2.o: poisson_functions.o
poisson2: $(OBJS) poisson2.o
	$(FC) $(FFLAGS) -o $@ poisson2.o $(OBJS) $(LIBS)

poisson20.o: poisson_functions.o
poisson20: $(OBJS) poisson20.o
	$(FC) $(FFLAGS) -o $@ poisson20.o $(OBJS) $(LIBS)

poisson21.o: functions.o poisson_functions.o
poisson21: $(OBJS) poisson21.o
	$(FC) $(FFLAGS) -o $@ poisson21.o $(OBJS) $(LIBS)

poisson22.o: functions.o poisson_functions.o
poisson22: $(OBJS) poisson22.o
	$(FC) $(FFLAGS) -o $@ poisson22.o $(OBJS) $(LIBS)

poisson23.o: poisson_functions.o
poisson23: $(OBJS) poisson23.o
	$(FC) $(FFLAGS) -o $@ poisson23.o $(OBJS) $(LIBS)

poisson27.o: poisson_functions.o
poisson27: $(OBJS) poisson27.o
	$(FC) $(FFLAGS) -o $@ poisson27.o $(OBJS) $(LIBS)

poisson28.o: poisson_functions.o
poisson28: $(OBJS) poisson28.o
	$(FC) $(FFLAGS) -o $@ poisson28.o $(OBJS) $(LIBS)

poisson29.o: poisson_functions.o
poisson29: $(OBJS) poisson29.o
	$(FC) $(FFLAGS) -o $@ poisson29.o $(OBJS) $(LIBS)

poisson3.o: functions.o poisson_functions.o
poisson3: $(OBJS) poisson3.o
	$(FC) $(FFLAGS) -o $@ poisson3.o $(OBJS) $(LIBS)

poisson30.o: poisson_functions.o
poisson30: $(OBJS) poisson30.o
	$(FC) $(FFLAGS) -o $@ poisson30.o $(OBJS) $(LIBS)

poisson31.o: poisson_functions.o
poisson31: $(OBJS) poisson31.o
	$(FC) $(FFLAGS) -o $@ poisson31.o $(OBJS) $(LIBS)

poisson32.o: poisson_functions.o
poisson32: $(OBJS) poisson32.o
	$(FC) $(FFLAGS) -o $@ poisson32.o $(OBJS) $(LIBS)

poisson33.o: poisson_functions.o
poisson33: $(OBJS) poisson33.o
	$(FC) $(FFLAGS) -o $@ poisson33.o $(OBJS) $(LIBS)

poisson34.o: poisson_functions.o
poisson34: $(OBJS) poisson34.o
	$(FC) $(FFLAGS) -o $@ poisson34.o $(OBJS) $(LIBS)

poisson35.o: functions.o poisson_functions.o
poisson35: $(OBJS) poisson35.o
	$(FC) $(FFLAGS) -o $@ poisson35.o $(OBJS) $(LIBS)

poisson36.o: functions.o poisson_functions.o
poisson36: $(OBJS) poisson36.o
	$(FC) $(FFLAGS) -o $@ poisson36.o $(OBJS) $(LIBS)

poisson37.o: functions.o poisson_functions.o
poisson37: $(OBJS) poisson37.o
	$(FC) $(FFLAGS) -o $@ poisson37.o $(OBJS) $(LIBS)

poisson4.o: functions.o poisson_functions.o
poisson4: $(OBJS) poisson4.o
	$(FC) $(FFLAGS) -o $@ poisson4.o $(OBJS) $(LIBS)

poisson5.o: poisson_functions.o
poisson5: $(OBJS) poisson5.o
	$(FC) $(FFLAGS) -o $@ poisson5.o $(OBJS) $(LIBS)

poisson6.o: poisson_functions.o
poisson6: $(OBJS) poisson6.o
	$(FC) $(FFLAGS) -o $@ poisson6.o $(OBJS) $(LIBS)

poisson7.o: poisson_functions.o
poisson7: $(OBJS) poisson7.o
	$(FC) $(FFLAGS) -o $@ poisson7.o $(OBJS) $(LIBS)

poisson8.o: poisson_functions.o
poisson8: $(OBJS) poisson8.o
	$(FC) $(FFLAGS) -o $@ poisson8.o $(OBJS) $(LIBS)

poisson9: $(OBJS) poisson9.o
	$(FC) $(FFLAGS) -o $@ poisson9.o $(OBJS) $(LIBS)

poisson14.o: poisson_functions.o
poisson14: $(OBJS) poisson14.o
	$(FC) $(FFLAGS) -o $@ poisson14.o $(OBJS) $(LIBS)

poisson15.o: poisson_functions.o
poisson15: $(OBJS) poisson15.o
	$(FC) $(FFLAGS) -o $@ poisson15.o $(OBJS) $(LIBS)

poisson16.o: poisson_functions.o
poisson16: $(OBJS) poisson16.o
	$(FC) $(FFLAGS) -o $@ poisson16.o $(OBJS) $(LIBS)

poisson17.o: poisson_functions.o
poisson17: $(OBJS) poisson17.o
	$(FC) $(FFLAGS) -o $@ poisson17.o $(OBJS) $(LIBS)

poisson18.o: poisson_functions.o
poisson18: $(OBJS) poisson18.o
	$(FC) $(FFLAGS) -o $@ poisson18.o $(OBJS) $(LIBS)

poisson19.o: poisson_functions.o
poisson19: $(OBJS) poisson19.o
	$(FC) $(FFLAGS) -o $@ poisson19.o $(OBJS) $(LIBS)

poisson24.o: poisson_functions.o
poisson24: $(OBJS) poisson24.o
	$(FC) $(FFLAGS) -o $@ poisson24.o $(OBJS) $(LIBS)

poisson25.o: poisson_functions.o
poisson25: $(OBJS) poisson25.o
	$(FC) $(FFLAGS) -o $@ poisson25.o $(OBJS) $(LIBS)

poisson26.o: poisson_functions.o
poisson26: $(OBJS) poisson26.o
	$(FC) $(FFLAGS) -o $@ poisson26.o $(OBJS) $(LIBS)

clean:
	$(RM)  poisson1 poisson12 poisson13 poisson1a poisson1b poisson1c poisson2 poisson20 poisson21 poisson22 poisson23 poisson27 poisson28 poisson29 poisson3 poisson30 poisson31 poisson32 poisson33 poisson34 poisson35 poisson36 poisson37 poisson4 poisson5 poisson6 poisson7 poisson8 poisson9 poisson14 poisson15 poisson16 poisson17 poisson18 poisson19 poisson24 poisson25 poisson26 *.o *.mod *.bin tags

outclean:
	$(RM) $(OUTCLEAN)

allclean: clean outclean
