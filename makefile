COMPILER = ifort
LD = ifort
FCFLAGS = -O3 -mkl

EXE = cpmc
OBJS = cpmc.o SYSTEM_mod.o ONE_BODY_mod.o INITIALIZATION_mod.o  WALK_mod.o
MODS = SYSTEM_mod.mod ONE_BODY_mod.mod INITIALIZATION_mod.mod WALK_mod.mod
 
$(EXE): $(OBJS)
	${LD} ${FCFLAGS} -o $(EXE) $(OBJS)

SYSTEM_mod.o SYSTEM_mod.mod: SYSTEM_mod.f90
	$(COMPILER) $(FCFLAGS) -c SYSTEM_mod.f90 $(LIBS)

ONE_BODY_mod.o ONE_BODY_mod.mod: SYSTEM_mod.mod ONE_BODY_mod.f90
	$(COMPILER) $(FCFLAGS) -c ONE_BODY_mod.f90 $(LIBS)

INITIALIZATION_mod.o INITIALIZATION_mod.mod: SYSTEM_mod.mod ONE_BODY_mod.mod \
        INITIALIZATION_mod.f90
	$(COMPILER) $(FCFLAGS) -c INITIALIZATION_mod.f90 $(LIBS)

WALK_mod.o WALK_mod.mod: SYSTEM_mod.mod ONE_BODY_mod.mod \
        INITIALIZATION_mod.mod WALK_mod.f90
	$(COMPILER) $(FCFLAGS) -c WALK_mod.f90 $(LIBS)

cpmc.o: $(MODS) cpmc.f90
	$(COMPILER) $(FCFLAGS) -c cpmc.f90 $(LIBS)

clean:
	rm -f core *.o *.mod

realclean: clean
	rm -f $(EXE)
