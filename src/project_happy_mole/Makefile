include Makefile.inc
CODE = project_happy_mole
PROG = ../../build/$(CODE)/$(CODE)
OBJECT_PATH=../../build/$(CODE)/

OBJS = $(OBJECT_PATH)main.o\
 $(OBJECT_PATH)options.o\
 $(OBJECT_PATH)scatteringrates.o\
 $(OBJECT_PATH)lineshape_helper.o\
 $(OBJECT_PATH)lo_evaluate_phonon_self_energy.o\
 $(OBJECT_PATH)lo_evaluate_phonon_self_energy_generate.o\
 $(OBJECT_PATH)lo_evaluate_phonon_self_energy_isotope.o\
 $(OBJECT_PATH)lo_evaluate_phonon_self_energy_threephonon.o\
 $(OBJECT_PATH)lo_evaluate_phonon_self_energy_fourphonon.o\
 $(OBJECT_PATH)lo_tetrahedron_interpolation.o\
 $(OBJECT_PATH)lo_selfenergy_interpolation.o\
 $(OBJECT_PATH)lo_selfenergy_interpolation_evaluate.o\
 $(OBJECT_PATH)lo_selfenergy_interpolation_path.o\
 $(OBJECT_PATH)scf_helper.o\
 $(OBJECT_PATH)create_selfenergy_interpolation.o

LPATH = -L../../lib $(blaslapackLPATH) $(incLPATHhdf) $(incLPATHmpi) $(incLPATHfft)
IPATH = -I../../inc/libolle -I../../inc/libflap $(blaslapackIPATH) ${incIPATHhdf} $(incIPATHmpi) $(incIPATHfft)
LIBS = -lolle -lflap $(blaslapackLIBS) ${incLIBShdf} $(incLIBSmpi) $(incLIBSfft)

#OPT = -Ofast
OPT = -O0 -fbacktrace -fcheck=all -finit-real=nan -finit-derived -fmax-errors=10 -Wall -fcheck=no-recursion
F90 = $(FC) $(LPATH) $(IPATH) $(MODULE_FLAG) $(OBJECT_PATH) #$(warnings_gcc)
F90FLAGS = $(OPT) $(MODS) $(LIBS)

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) $(OBJECT_PATH)*.mod

$(OBJECT_PATH)main.o:\
	../../lib/libolle.a\
	main.f90\
	$(OBJECT_PATH)options.o\
	$(OBJECT_PATH)scatteringrates.o\
	$(OBJECT_PATH)lineshape_helper.o\
	$(OBJECT_PATH)scf_helper.o\
	$(OBJECT_PATH)create_selfenergy_interpolation.o\
	$(OBJECT_PATH)lo_evaluate_phonon_self_energy.o\
	$(OBJECT_PATH)lo_selfenergy_interpolation.o
	$(F90) $(F90FLAGS) -c main.f90 $(LIBS) -o $@

$(OBJECT_PATH)scatteringrates.o:\
	scatteringrates.f90\
	$(OBJECT_PATH)lineshape_helper.o
	$(F90) $(F90FLAGS) -c scatteringrates.f90 $(LIBS) -o $@

$(OBJECT_PATH)lineshape_helper.o:\
	lineshape_helper.f90
	$(F90) $(F90FLAGS) -c lineshape_helper.f90 $(LIBS) -o $@

$(OBJECT_PATH)options.o:\
	options.f90
	$(F90) $(F90FLAGS) -c options.f90 $(LIBS) -o $@

$(OBJECT_PATH)scf_helper.o:\
	scf_helper.f90\
	$(OBJECT_PATH)lineshape_helper.o
	$(F90) $(F90FLAGS) -c scf_helper.f90 $(LIBS) -o $@

$(OBJECT_PATH)create_selfenergy_interpolation.o:\
	create_selfenergy_interpolation.f90\
	$(OBJECT_PATH)lo_evaluate_phonon_self_energy.o\
	$(OBJECT_PATH)lineshape_helper.o
	$(F90) $(F90FLAGS) -c create_selfenergy_interpolation.f90 $(LIBS) -o $@

$(OBJECT_PATH)lo_tetrahedron_interpolation.o:\
	lo_tetrahedron_interpolation.f90
	$(F90) $(F90FLAGS) -c lo_tetrahedron_interpolation.f90 $(LIBS) -o $@

$(OBJECT_PATH)lo_selfenergy_interpolation.o:\
	lo_selfenergy_interpolation.f90\
	$(OBJECT_PATH)lo_tetrahedron_interpolation.o
	$(F90) $(F90FLAGS) -c lo_selfenergy_interpolation.f90 $(LIBS) -o $@
$(OBJECT_PATH)lo_selfenergy_interpolation_evaluate.o:\
	lo_selfenergy_interpolation_evaluate.f90\
	$(OBJECT_PATH)lo_selfenergy_interpolation.o
	$(F90) $(F90FLAGS) -c lo_selfenergy_interpolation_evaluate.f90 $(LIBS) -o $@
$(OBJECT_PATH)lo_selfenergy_interpolation_path.o:\
	lo_selfenergy_interpolation_path.f90\
	$(OBJECT_PATH)lo_selfenergy_interpolation.o
	$(F90) $(F90FLAGS) -c lo_selfenergy_interpolation_path.f90 $(LIBS) -o $@

$(OBJECT_PATH)lo_evaluate_phonon_self_energy.o:\
	lo_evaluate_phonon_self_energy.f90\
	$(OBJECT_PATH)scatteringrates.o\
	$(OBJECT_PATH)options.o\
	$(OBJECT_PATH)lo_selfenergy_interpolation.o\
	$(OBJECT_PATH)lineshape_helper.o
	$(F90) $(F90FLAGS) -c lo_evaluate_phonon_self_energy.f90 $(LIBS) -o $@
$(OBJECT_PATH)lo_evaluate_phonon_self_energy_generate.o:\
	lo_evaluate_phonon_self_energy_generate.f90\
	$(OBJECT_PATH)lo_evaluate_phonon_self_energy.o
	$(F90) $(F90FLAGS) -c lo_evaluate_phonon_self_energy_generate.f90 $(LIBS) -o $@
$(OBJECT_PATH)lo_evaluate_phonon_self_energy_gaussian.o:\
	lo_evaluate_phonon_self_energy_gaussian.f90\
	$(OBJECT_PATH)lo_evaluate_phonon_self_energy.o
	$(F90) $(F90FLAGS) -c lo_evaluate_phonon_self_energy_gaussian.f90 $(LIBS) -o $@
$(OBJECT_PATH)lo_evaluate_phonon_self_energy_isotope.o:\
	lo_evaluate_phonon_self_energy_gaussian.f90\
	$(OBJECT_PATH)lo_evaluate_phonon_self_energy.o
	$(F90) $(F90FLAGS) -c lo_evaluate_phonon_self_energy_isotope.f90 $(LIBS) -o $@
$(OBJECT_PATH)lo_evaluate_phonon_self_energy_threephonon.o:\
	lo_evaluate_phonon_self_energy_threephonon.f90\
	$(OBJECT_PATH)lo_evaluate_phonon_self_energy.o
	$(F90) $(F90FLAGS) -c lo_evaluate_phonon_self_energy_threephonon.f90 $(LIBS) -o $@
$(OBJECT_PATH)lo_evaluate_phonon_self_energy_fourphonon.o:\
	lo_evaluate_phonon_self_energy_fourphonon.f90\
	$(OBJECT_PATH)lo_evaluate_phonon_self_energy.o
	$(F90) $(F90FLAGS) -c lo_evaluate_phonon_self_energy_fourphonon.f90 $(LIBS) -o $@
