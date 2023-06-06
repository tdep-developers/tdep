include Makefile.inc
CODE = extract_forceconstants
PROG = ../../build/$(CODE)/$(CODE)
OBJECT_PATH=../../build/$(CODE)/

OBJS = \
$(OBJECT_PATH)main.o\
$(OBJECT_PATH)options.o

LPATH = -L../../lib $(blaslapackLPATH) $(incLPATHhdf) $(incLPATHmpi)
IPATH = -I../../inc/libolle -I../../inc/libflap $(blaslapackIPATH) $(incLPATHhdf) $(incIPATHmpi)
LIBS = -lolle -lflap $(blaslapackLIBS) ${incIPATHhdf} $(incLIBShdf) $(incLIBSmpi)

#OPT = -O0 -fcheck=all -fbacktrace -finit-real=nan -finit-derived -fmax-errors=10
F90 = $(FC) $(LPATH) $(IPATH) $(MODULE_FLAG) $(OBJECT_PATH)
F90FLAGS = $(OPT) $(MODS) $(LIBS)

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) $(OBJECT_PATH)*.mod

.f90.o:
	$(F90) $(F90FLAGS) -c $< $(LIBS)

$(OBJECT_PATH)main.o: $(OBJECT_PATH)options.o
	$(F90) $(F90FLAGS) -c main.f90 $(LIBS) -o $@
$(OBJECT_PATH)options.o:
	$(F90) $(F90FLAGS) -c options.f90 $(LIBS) -o $@
