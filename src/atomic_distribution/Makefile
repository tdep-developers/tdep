include Makefile.inc
CODE = atomic_distribution
PROG = ../../build/$(CODE)/$(CODE)
OBJECT_PATH=../../build/$(CODE)/

OBJS = $(OBJECT_PATH)main.o $(OBJECT_PATH)options.o \
       $(OBJECT_PATH)pairmapping.o $(OBJECT_PATH)pair_distribution.o $(OBJECT_PATH)mean_square_displacement.o \
       $(OBJECT_PATH)diffraction.o

LPATH = -L../../lib $(blaslapackLPATH) $(incLPATHhdf) $(incLPATHfft) $(incLPATHmpi)
IPATH = -I../../inc/libolle -I../../inc/libflap $(blaslapackIPATH) $(incIPATHhdf) $(incIPATHfft) $(incIPATHmpi)
LIBS = ../../lib/libolle.a -lflap $(blaslapackLIBS) $(incLIBShdf) $(incLIBSfft) $(incLIBSmpi)

#OPT = -O0 -fcheck=all -fbacktrace -finit-real=nan -finit-derived -fmax-errors=10 #-Wall
F90 = $(FC) $(LPATH) $(IPATH) $(MODULE_FLAG) $(OBJECT_PATH) #$(warnings_gcc)
F90FLAGS = $(OPT) $(MODS) $(LIBS)

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) $(OBJECT_PATH)*.mod

$(OBJECT_PATH)main.o: main.f90 $(OBJECT_PATH)options.o $(OBJECT_PATH)pairmapping.o $(OBJECT_PATH)pair_distribution.o $(OBJECT_PATH)mean_square_displacement.o $(OBJECT_PATH)diffraction.o
	$(F90) $(F90FLAGS) -c main.f90 $(LIBS) -o $@
$(OBJECT_PATH)pair_distribution.o: pair_distribution.f90 $(OBJECT_PATH)pairmapping.o
	$(F90) $(F90FLAGS) -c pair_distribution.f90 $(LIBS) -o $@
$(OBJECT_PATH)mean_square_displacement.o: mean_square_displacement.f90 $(OBJECT_PATH)pairmapping.o
	$(F90) $(F90FLAGS) -c mean_square_displacement.f90 $(LIBS) -o $@
$(OBJECT_PATH)diffraction.o: diffraction.f90 $(OBJECT_PATH)pairmapping.o
	$(F90) $(F90FLAGS) -c diffraction.f90 $(LIBS) -o $@
$(OBJECT_PATH)pairmapping.o: pairmapping.f90
	$(F90) $(F90FLAGS) -c pairmapping.f90 $(LIBS) -o $@
$(OBJECT_PATH)options.o: options.f90
	$(F90) $(F90FLAGS) -c options.f90 $(LIBS) -o $@
