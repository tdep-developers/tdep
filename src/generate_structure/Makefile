include Makefile.inc
CODE = generate_structure
PROG = ../../build/$(CODE)/$(CODE)
OBJECT_PATH=../../build/$(CODE)/

# OBJS = $(OBJECT_PATH)main.o $(OBJECT_PATH)options.o $(OBJECT_PATH)magneticdisorder.o $(OBJECT_PATH)autocell.o 
OBJS = $(OBJECT_PATH)main.o $(OBJECT_PATH)options.o $(OBJECT_PATH)autocell.o 


#OPT = -O0 -fbacktrace -fcheck=all -finit-real=nan -finit-derived -fmax-errors=10
LPATH = -L../../lib $(blaslapackLPATH) $(incLPATHhdf) $(incLPATHmpi)
IPATH = -I../../inc/libolle -I../../inc/libflap $(blaslapackIPATH) $(incIPATHhdf) $(incIPATHmpi)
LIBS = ../../lib/libolle.a -lflap $(blaslapackLIBS) $(incLIBShdf) $(incLIBSmpi)

F90 = $(FC) $(LPATH) $(IPATH) $(MODULE_FLAG) $(OBJECT_PATH) # $(warnings_gcc)
F90FLAGS = $(OPT) $(MODS) $(LIBS)

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) $(OBJECT_PATH)*.mod

# $(OBJECT_PATH)main.o: $(OBJECT_PATH)options.o $(OBJECT_PATH)magneticdisorder.o $(OBJECT_PATH)autocell.o 
$(OBJECT_PATH)main.o: $(OBJECT_PATH)options.o  $(OBJECT_PATH)autocell.o 
	$(F90) $(F90FLAGS) -c main.f90 $(LIBS) -o $@
$(OBJECT_PATH)autocell.o:
	$(F90) $(F90FLAGS) -c autocell.f90 $(LIBS) -o $@
$(OBJECT_PATH)options.o:
	$(F90) $(F90FLAGS) -c options.f90 $(LIBS) -o $@
# $(OBJECT_PATH)magneticdisorder.o:
# 	$(F90) $(F90FLAGS) -c magneticdisorder.f90 $(LIBS) -o $@

