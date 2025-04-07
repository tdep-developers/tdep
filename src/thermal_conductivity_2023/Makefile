include Makefile.inc
CODE = thermal_conductivity_2023
PROG = ../../build/$(CODE)/$(CODE)
OBJECT_PATH=../../build/$(CODE)/

OBJS = \
$(OBJECT_PATH)main.o\
$(OBJECT_PATH)options.o\
$(OBJECT_PATH)scatteringstrengths.o\
$(OBJECT_PATH)pbe.o\
$(OBJECT_PATH)phononevents.o\
$(OBJECT_PATH)mfp.o

LPATH = -L../../lib $(blaslapackLPATH) $(incLPATHmpi) $(incLPATHhdf)
IPATH = -I../../inc/libolle -I../../inc/libflap $(blaslapackIPATH) $(incIPATHmpi) $(incIPATHhdf)
LIBS = ../../lib/libolle.a -lflap $(blaslapackLIBS) $(incLIBSmpi) $(incLIBShdf)

#OPT = -O0 -fbacktrace -fcheck=all -finit-real=nan -finit-derived
F90 = $(FC) $(LPATH) $(IPATH) $(MODULE_FLAG) $(OBJECT_PATH) #$(warnings_gcc)
F90FLAGS = $(OPT) $(MODS) $(LIBS) 

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) $(OBJECT_PATH)*.mod

$(OBJECT_PATH)main.o: \
$(OBJECT_PATH)options.o\
$(OBJECT_PATH)scatteringstrengths.o\
$(OBJECT_PATH)pbe.o\
$(OBJECT_PATH)phononevents.o\
$(OBJECT_PATH)mfp.o
	$(F90) $(OPT) $(F90FLAGS) -c main.f90 $(LIBS) -o $@
$(OBJECT_PATH)scatteringstrengths.o: $(OBJECT_PATH)phononevents.o
	$(F90) $(OPT) $(F90FLAGS) -c scatteringstrengths.f90 $(LIBS) -o $@
$(OBJECT_PATH)pbe.o: $(OBJECT_PATH)phononevents.o
	$(F90) $(OPT) $(F90FLAGS) -c pbe.f90 $(LIBS) -o $@
$(OBJECT_PATH)phononevents.o:
	$(F90) $(OPT) $(F90FLAGS) -c phononevents.f90 $(LIBS) -o $@
$(OBJECT_PATH)mfp.o: 
	$(F90) $(OPT) $(F90FLAGS) -c mfp.f90 $(LIBS) -o $@
$(OBJECT_PATH)options.o:
	$(F90) $(OPT) $(F90FLAGS) -c options.f90 $(LIBS) -o $@

