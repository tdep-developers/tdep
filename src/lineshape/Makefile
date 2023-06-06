include Makefile.inc
CODE = lineshape
PROG = ../../build/$(CODE)/$(CODE)
OBJECT_PATH=../../build/$(CODE)/

OBJS = $(OBJECT_PATH)main.o\
 $(OBJECT_PATH)options.o\
 $(OBJECT_PATH)io.o\
 $(OBJECT_PATH)scatteringrates.o\
 $(OBJECT_PATH)phonondamping.o\
 $(OBJECT_PATH)phonondamping_aux.o\
 $(OBJECT_PATH)phonondamping_dos.o\
 $(OBJECT_PATH)phonondamping_grid.o\
 $(OBJECT_PATH)phonondamping_path.o\
 $(OBJECT_PATH)dielscatter.o\
 $(OBJECT_PATH)dielscatter_helper.o\
 $(OBJECT_PATH)lineshape_helper.o\
 $(OBJECT_PATH)lo_realspace_selfenergy.o\
 $(OBJECT_PATH)lo_thermal_transport.o

LPATH = -L../../lib $(blaslapackLPATH) $(incLPATHhdf) $(incLPATHmpi) $(incLPATHfft)
IPATH = -I../../inc/libolle -I../../inc/libflap $(blaslapackIPATH) ${incIPATHhdf} $(incIPATHmpi) $(incIPATHfft)
LIBS = -lolle -lflap $(blaslapackLIBS) ${incLIBShdf} $(incLIBSmpi) $(incLIBSfft)

#OPT = -Ofast
#OPT = -O0 -fbacktrace -fcheck=all -finit-real=nan -finit-derived -fmax-errors=10 -Wall
F90 = $(FC) $(LPATH) $(IPATH) $(MODULE_FLAG) $(OBJECT_PATH) #$(warnings_gcc)
F90FLAGS = $(OPT) $(MODS) $(LIBS)

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) $(OBJECT_PATH)*.mod

$(OBJECT_PATH)main.o:\
	$(OBJECT_PATH)options.o\
	$(OBJECT_PATH)io.o\
	$(OBJECT_PATH)scatteringrates.o\
	$(OBJECT_PATH)phonondamping.o\
	$(OBJECT_PATH)dielscatter.o\
	$(OBJECT_PATH)lineshape_helper.o\
	$(OBJECT_PATH)lo_realspace_selfenergy.o\
	$(OBJECT_PATH)lo_thermal_transport.o
	$(F90) $(F90FLAGS) -c main.f90 $(LIBS) -o $@
$(OBJECT_PATH)io.o: $(OBJECT_PATH)phonondamping.o $(OBJECT_PATH)dielscatter.o $(OBJECT_PATH)lineshape_helper.o
	$(F90) $(F90FLAGS) -c io.f90 $(LIBS) -o $@
$(OBJECT_PATH)phonondamping.o:\
	phonondamping.f90\
	$(OBJECT_PATH)scatteringrates.o\
	$(OBJECT_PATH)options.o\
	$(OBJECT_PATH)lineshape_helper.o\
	$(OBJECT_PATH)lo_realspace_selfenergy.o\
	$(OBJECT_PATH)lo_thermal_transport.o
	$(F90) $(F90FLAGS) -c phonondamping.f90 $(LIBS) -o $@
$(OBJECT_PATH)phonondamping_aux.o: $(OBJECT_PATH)phonondamping.o
	$(F90) $(F90FLAGS) -c phonondamping_aux.f90 $(LIBS) -o $@
$(OBJECT_PATH)phonondamping_dos.o: $(OBJECT_PATH)phonondamping.o
	$(F90) $(F90FLAGS) -c phonondamping_dos.f90 $(LIBS) -o $@
$(OBJECT_PATH)phonondamping_grid.o: $(OBJECT_PATH)phonondamping.o
	$(F90) $(F90FLAGS) -c phonondamping_grid.f90 $(LIBS) -o $@
$(OBJECT_PATH)phonondamping_path.o: $(OBJECT_PATH)phonondamping.o
	$(F90) $(F90FLAGS) -c phonondamping_path.f90 $(LIBS) -o $@	
$(OBJECT_PATH)scatteringrates.o: $(OBJECT_PATH)lineshape_helper.o
	$(F90) $(F90FLAGS) -c scatteringrates.f90 $(LIBS) -o $@
$(OBJECT_PATH)dielscatter.o:\
	dielscatter.f90\
	$(OBJECT_PATH)options.o\
	$(OBJECT_PATH)phonondamping.o\
	$(OBJECT_PATH)lineshape_helper.o
	$(F90) $(F90FLAGS) -c dielscatter.f90 $(LIBS) -o $@
$(OBJECT_PATH)dielscatter_helper.o:	$(OBJECT_PATH)dielscatter.o
	$(F90) $(F90FLAGS) -c dielscatter_helper.f90 $(LIBS) -o $@
$(OBJECT_PATH)lineshape_helper.o:\
	lineshape_helper.f90
	$(F90) $(F90FLAGS) -c lineshape_helper.f90 $(LIBS) -o $@
$(OBJECT_PATH)lo_realspace_selfenergy.o:
	$(F90) $(F90FLAGS) -c lo_realspace_selfenergy.f90 $(LIBS) -o $@
$(OBJECT_PATH)lo_thermal_transport.o: $(OBJECT_PATH)lineshape_helper.o
	$(F90) $(F90FLAGS) -c lo_thermal_transport.f90 $(LIBS) -o $@
$(OBJECT_PATH)options.o:
	$(F90) $(F90FLAGS) -c options.f90 $(LIBS) -o $@
