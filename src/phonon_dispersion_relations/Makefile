# I made this so I have this line under control.
include Makefile.inc
# ok, I understand this one.
CODE = phonon_dispersion_relations
# this one as well
PROG = ../../build/$(CODE)/$(CODE)
OBJECT_PATH=../../build/$(CODE)/

# list of stuff. Used for something at least.
OBJS = \
$(OBJECT_PATH)main.o\
$(OBJECT_PATH)options.o\
$(OBJECT_PATH)densityplots.o\
$(OBJECT_PATH)densityplots_stuntscattering.o\
$(OBJECT_PATH)unfold_phonons.o\
$(OBJECT_PATH)type_fast_interpolation.o\
$(OBJECT_PATH)velocitydos.o\
$(OBJECT_PATH)activity.o

# some paths and stuff. Ok.
LPATH = -L../../lib $(blaslapackLPATH) $(incLPATHhdf) $(incLPATHmpi) $(incLPATHfft)
IPATH = -I../../inc/libolle -I../../inc/libflap $(blaslapackIPATH) $(incIPATHhdf) $(incIPATHmpi) $(incIPATHfft)
LIBS = ../../lib/libolle.a -lflap $(blaslapackLIBS) $(incLIBShdf) $(incLIBSmpi) $(incLIBSfft)

# ok, I think I get this.
#OPT = -O0 -fbacktrace -fcheck=all -finit-real=nan -finit-derived -fmax-errors=10
F90 = $(FC) $(LPATH) $(IPATH) $(MODULE_FLAG) $(OBJECT_PATH) #$(warnings_gcc)
F90FLAGS = $(OPT) $(MODS) $(LIBS)

# I guess this is the normal one
all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

# ok, I understand this
clean:
	rm -f $(PROG) $(OBJS) *.mod

# This I understand, but it seems stupid.
$(OBJECT_PATH)main.o:\
    $(OBJECT_PATH)options.o\
    $(OBJECT_PATH)densityplots.o\
    $(OBJECT_PATH)unfold_phonons.o\
    $(OBJECT_PATH)type_fast_interpolation.o\
    $(OBJECT_PATH)velocitydos.o\
	$(OBJECT_PATH)activity.o
	$(F90) $(F90FLAGS) -c main.f90 $(LIBS) -o $@
$(OBJECT_PATH)densityplots.o: $(OBJECT_PATH)options.o
	$(F90) $(F90FLAGS) -c densityplots.f90 $(LIBS) -o $@
$(OBJECT_PATH)densityplots_stuntscattering.o: $(OBJECT_PATH)densityplots.o
	$(F90) $(F90FLAGS) -c densityplots_stuntscattering.f90 $(LIBS) -o $@
$(OBJECT_PATH)velocitydos.o:
	$(F90) $(F90FLAGS) -c velocitydos.f90 $(LIBS) -o $@
$(OBJECT_PATH)unfold_phonons.o: $(OBJECT_PATH)type_fast_interpolation.o
	$(F90) $(F90FLAGS) -c unfold_phonons.f90 $(LIBS) -o $@
$(OBJECT_PATH)type_fast_interpolation.o:
	$(F90) $(F90FLAGS) -c type_fast_interpolation.f90 $(LIBS) -o $@
$(OBJECT_PATH)activity.o:
	$(F90) $(F90FLAGS) -c activity.f90 $(LIBS) -o $@	
$(OBJECT_PATH)options.o:
	$(F90) $(F90FLAGS) -c options.f90 $(LIBS) -o $@
