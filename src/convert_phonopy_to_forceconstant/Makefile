include Makefile.inc
CODE = convert_phonopy_to_forceconstant
PROG = ../../build/$(CODE)/$(CODE)
OBJECT_PATH=../../build/$(CODE)/

OBJS = $(OBJECT_PATH)main.o $(OBJECT_PATH)options.o $(OBJECT_PATH)io.o

LPATH = -L../../lib $(blaslapackLPATH) $(incLPATHhdf) $(incLPATHmpi)
IPATH = -I../../inc/libolle -I../../inc/libflap $(blaslapackIPATH) $(incIPATHhdf) $(incIPATHmpi)
LIBS = -lolle -lflap $(blaslapackLIBS) $(incLIBShdf) $(incLIBSmpi)

#OPT = -O0 -fbacktrace -fcheck=all -finit-real=nan -finit-derived -fmax-errors=10

F90 = $(FC) $(LPATH) $(IPATH) $(MODULE_FLAG) $(OBJECT_PATH) 
F90FLAGS = $(OPT) $(MODS) $(LIBS)

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.f90.o:
	$(F90) $(F90FLAGS) -c $< $(LIBS)

$(OBJECT_PATH)main.o: $(OBJECT_PATH)options.o $(OBJECT_PATH)io.o
	$(F90) $(F90FLAGS) -c main.f90 $(LIBS) -o $@
$(OBJECT_PATH)io.o:
	$(F90) $(F90FLAGS) -c io.f90 $(LIBS) -o $@
$(OBJECT_PATH)options.o:
	$(F90) $(F90FLAGS) -c options.f90 $(LIBS) -o $@
