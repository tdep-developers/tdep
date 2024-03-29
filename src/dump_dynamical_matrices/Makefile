include Makefile.inc
CODE = dump_dynamical_matrices
PROG = ../../build/$(CODE)/$(CODE)
OBJECT_PATH=../../build/$(CODE)/

OBJS = $(OBJECT_PATH)main.o $(OBJECT_PATH)options.o  $(OBJECT_PATH)ddb_io_out.o 


# some paths and stuff. Ok.
IPATH = -I../../inc/libolle -I../../inc/libflap $(blaslapackIPATH) $(incIPATHhdf) $(incIPATHmpi) $(incIPATHfft)

LPATH = -L../../lib $(blaslapackLPATH) $(incLPATHhdf) $(incLPATHmpi) $(incLPATHfft)
LIBS = -lolle -lflap $(blaslapackLIBS) $(incLIBShdf) $(incLIBSmpi) $(incLIBSfft)

# ok, I think I get this.
#OPT = -O0 -fbacktrace -fcheck=all -finit-real=nan -finit-derived -fmax-errors=10
F90 = $(FC) $(IPATH) $(MODULE_FLAG) $(OBJECT_PATH) #$(warnings_gcc)
F90FLAGS = $(OPT) $(MODS) $(LIBS)

# I guess this is the normal one
all: $(PROG)

$(PROG): $(OBJS)
	$(F90) -o $@ $(OBJS) $(LPATH) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod $(OBJECT_PATH)*.mod


$(OBJECT_PATH)main.o: $(OBJECT_PATH)options.o  $(OBJECT_PATH)ddb_io_out.o
	$(F90) $(F90FLAGS) -c main.f90  -o $@
$(OBJECT_PATH)options.o:
	$(F90) $(F90FLAGS) -c options.f90 -o $@
$(OBJECT_PATH)ddb_io_out.o:
	$(F90) $(F90FLAGS) -c ddb_io_out.F90 -o $@

