include Makefile.inc
CODE = refine_structure
PROG = ../../build/$(CODE)/$(CODE)
OBJECT_PATH=../../build/$(CODE)/

OBJS = $(OBJECT_PATH)main.o\
$(OBJECT_PATH)refine.o\
$(OBJECT_PATH)options.o\
$(OBJECT_PATH)lo_spacegroup.o\
$(OBJECT_PATH)lo_spacegroup_genoperations.o\
$(OBJECT_PATH)lo_spacegroup_charactertable.o\
$(OBJECT_PATH)lo_spacegroup_applyoperation.o\
$(OBJECT_PATH)lo_spacegroup_helpers.o

LPATH = -L../../lib $(blaslapackLPATH) $(incLPATHhdf) $(incLPATHmpi)
IPATH = -I../../inc/libolle -I../../inc/libflap $(blaslapackIPATH) $(incIPATHhdf) $(incIPATHmpi)
LIBS = ../../lib/libolle.a -lflap $(blaslapackLIBS) $(incLIBShdf) $(incLIBSmpi)

#OPT = -Ofast
F90 = $(FC) $(LPATH) $(IPATH) $(MODULE_FLAG) $(OBJECT_PATH)
F90FLAGS = $(OPT) $(MODS) $(LIBS)

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.f90.o:
	$(F90) $(F90FLAGS) -c $< $(LIBS)

$(OBJECT_PATH)main.o: $(OBJECT_PATH)options.o $(OBJECT_PATH)refine.o $(OBJECT_PATH)lo_spacegroup.o
	$(F90) $(F90FLAGS) -c main.f90 $(LIBS) -o $@
$(OBJECT_PATH)lo_spacegroup.o:
	$(F90) $(F90FLAGS) -c lo_spacegroup.f90 $(LIBS) -o $@
$(OBJECT_PATH)lo_spacegroup_genoperations.o: $(OBJECT_PATH)lo_spacegroup.o
	$(F90) $(F90FLAGS) -c lo_spacegroup_genoperations.f90 $(LIBS) -o $@
$(OBJECT_PATH)lo_spacegroup_applyoperation.o: $(OBJECT_PATH)lo_spacegroup.o
	$(F90) $(F90FLAGS) -c lo_spacegroup_applyoperation.f90 $(LIBS) -o $@
$(OBJECT_PATH)lo_spacegroup_charactertable.o: $(OBJECT_PATH)lo_spacegroup.o
	$(F90) $(F90FLAGS) -c lo_spacegroup_charactertable.f90 $(LIBS) -o $@
$(OBJECT_PATH)lo_spacegroup_helpers.o: $(OBJECT_PATH)lo_spacegroup.o
	$(F90) $(F90FLAGS) -c lo_spacegroup_helpers.f90 $(LIBS) -o $@
$(OBJECT_PATH)refine.o: $(OBJECT_PATH)lo_spacegroup.o
	$(F90) $(F90FLAGS) -c refine.f90 $(LIBS) -o $@
$(OBJECT_PATH)options.o:
	$(F90) $(F90FLAGS) -c options.f90 $(LIBS) -o $@
