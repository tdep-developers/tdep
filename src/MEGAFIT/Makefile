include Makefile.inc
CODE = MEGAFIT
PROG = ../../build/$(CODE)/$(CODE)
OBJECT_PATH=../../build/$(CODE)/

OBJS = \
$(OBJECT_PATH)main.o\
$(OBJECT_PATH)options.o\
$(OBJECT_PATH)helperobjects.o\
$(OBJECT_PATH)gridenergy.o\
$(OBJECT_PATH)type_polynomial_interpolation.o\
$(OBJECT_PATH)type_gridsim.o\
$(OBJECT_PATH)diagnostics.o

LPATH = -L../../lib $(blaslapackLPATH) $(incLPATHhdf) $(incLPATHmpi)
IPATH = -I../../inc/libolle -I../../inc/libflap $(blaslapackIPATH) $(incLPATHhdf) $(incIPATHmpi)
LIBS = -lolle -lflap $(blaslapackLIBS) ${incIPATHhdf} $(incLIBShdf) $(incLIBSmpi)

#OPT = -O0 -fcheck=all -fbacktrace -fno-realloc-lhs -finit-real=nan -finit-derived 
#OPT = -O0 -fcheck=all -fbacktrace -fno-realloc-lhs -finit-real=nan -finit-derived -fmax-errors=10 -Werror #-g -pedantic -fbounds-check --warn-all #-ffpe-trap=zero,overflow,underflow,denormal
F90 = $(FC) $(LPATH) $(IPATH) $(MODULE_FLAG) $(OBJECT_PATH) # -DAGRESSIVE_SANITY # $(warnings_gcc)
F90FLAGS = $(OPT) $(MODS) $(LIBS)

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) $(OBJECT_PATH)*.mod

.f90.o:
	$(F90) $(F90FLAGS) -c $< $(LIBS)

$(OBJECT_PATH)main.o: $(OBJECT_PATH)options.o $(OBJECT_PATH)type_gridsim.o $(OBJECT_PATH)gridenergy.o $(OBJECT_PATH)diagnostics.o
	$(F90) $(F90FLAGS) -c main.f90 $(LIBS) -o $@
$(OBJECT_PATH)gridenergy.o: $(OBJECT_PATH)type_gridsim.o $(OBJECT_PATH)type_polynomial_interpolation.o
	$(F90) $(F90FLAGS) -c gridenergy.f90 $(LIBS) -o $@
$(OBJECT_PATH)diagnostics.o: $(OBJECT_PATH)type_gridsim.o $(OBJECT_PATH)type_polynomial_interpolation.o
	$(F90) $(F90FLAGS) -c diagnostics.f90 $(LIBS) -o $@
$(OBJECT_PATH)type_gridsim.o: $(OBJECT_PATH)type_polynomial_interpolation.o $(OBJECT_PATH)helperobjects.o
	$(F90) $(F90FLAGS) -c type_gridsim.f90 $(LIBS) -o $@
$(OBJECT_PATH)type_polynomial_interpolation.o:
	$(F90) $(F90FLAGS) -c type_polynomial_interpolation.f90 $(LIBS) -o $@
$(OBJECT_PATH)helperobjects.o:
	$(F90) $(F90FLAGS) -c helperobjects.f90 $(LIBS) -o $@	
$(OBJECT_PATH)options.o:
	$(F90) $(F90FLAGS) -c options.f90 $(LIBS) -o $@
