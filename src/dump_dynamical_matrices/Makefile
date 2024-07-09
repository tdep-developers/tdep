include Makefile.inc
CODE = dump_dynamical_matrices
PROG = ../../build/$(CODE)/$(CODE)
OBJECT_PATH=../../build/$(CODE)/

OBJS = $(OBJECT_PATH)main.o $(OBJECT_PATH)options.o

LPATH = -L../../lib $(blaslapackLPATH) $(incLPATHhdf) $(incLPATHmpi) $(incLPATHfft)
IPATH = -I../../inc/libolle -I../../inc/libflap $(blaslapackIPATH) $(incIPATHhdf) $(incIPATHmpi) $(incIPATHfft)
LIBS = -lolle -lflap $(blaslapackLIBS) $(incLIBShdf) $(incLIBSmpi) $(incLIBSfft)

#OPT = -O0 -fbacktrace -fcheck=all -Wall
F90 = $(FC) $(LPATH) $(IPATH) $(MODULE_FLAG) $(OBJECT_PATH)
F90FLAGS = $(OPT) $(MODS) $(LIBS)

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod $(OBJECT_PATH)*.mod


$(OBJECT_PATH)main.o: $(OBJECT_PATH)options.o 
	$(F90) $(F90FLAGS) -c main.f90 $(LIBS) -o $@
$(OBJECT_PATH)options.o:
	$(F90) $(F90FLAGS) -c options.f90 $(LIBS) -o $@

