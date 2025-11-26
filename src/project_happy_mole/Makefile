include Makefile.inc
CODE = project_happy_mole
PROG = ../../build/$(CODE)/$(CODE)
OBJECT_PATH=../../build/$(CODE)/

OBJS = $(OBJECT_PATH)main.o\
 $(OBJECT_PATH)options.o\
 $(OBJECT_PATH)lo_evaluate_phonon_self_energy.o\
 $(OBJECT_PATH)lo_evaluate_phonon_self_energy_generate.o\
 $(OBJECT_PATH)lo_evaluate_phonon_self_energy_threephonon.o\
 $(OBJECT_PATH)lo_evaluate_phonon_self_energy_fourphonon.o\
 $(OBJECT_PATH)lo_evaluate_phonon_self_energy_helpers.o\
 $(OBJECT_PATH)lo_tetrahedron_interpolation.o\
 $(OBJECT_PATH)lo_selfenergy_interpolation.o\
 $(OBJECT_PATH)lo_selfenergy_interpolation_evaluate.o\
 $(OBJECT_PATH)lo_selfenergy_interpolation_path.o\
 $(OBJECT_PATH)lo_selfenergy_interpolation_grid.o\
 $(OBJECT_PATH)lo_thermal_transport.o\
 $(OBJECT_PATH)lo_collision_matrix.o\
 $(OBJECT_PATH)lo_spectralfunction_convolution.o\
 $(OBJECT_PATH)lo_spectralfunction_helpers.o\
 $(OBJECT_PATH)lo_distributed_phonon_dispersion_relations.o\
 $(OBJECT_PATH)create_selfenergy_interpolation.o

LPATH = -L../../lib $(blaslapackLPATH) $(incLPATHhdf) $(incLPATHmpi) $(incLPATHfft)
IPATH = -I../../inc/libolle -I../../inc/libflap $(blaslapackIPATH) ${incIPATHhdf} $(incIPATHmpi) $(incIPATHfft)
LIBS = -lolle -lflap $(blaslapackLIBS) ${incLIBShdf} $(incLIBSmpi) $(incLIBSfft)

OPT = -Ofast #-fbacktrace
#OPT = -O0 -fbacktrace -fcheck=all -finit-real=nan -finit-derived -fmax-errors=10 -Wall -fcheck=no-recursion
F90 = $(FC) $(LPATH) $(IPATH) $(MODULE_FLAG) $(OBJECT_PATH) #$(warnings_gcc)
F90FLAGS = $(OPT) $(MODS) $(LIBS)
F90FLAGS_LOW = -O0 $(MODS) $(LIBS)

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) $(OBJECT_PATH)*.mod

$(OBJECT_PATH)main.o:\
	../../lib/libolle.a\
	main.f90\
	$(OBJECT_PATH)options.o\
	$(OBJECT_PATH)create_selfenergy_interpolation.o\
	$(OBJECT_PATH)lo_evaluate_phonon_self_energy.o\
	$(OBJECT_PATH)lo_thermal_transport.o\
	$(OBJECT_PATH)lo_collision_matrix.o\
	$(OBJECT_PATH)lo_spectralfunction_convolution.o\
	$(OBJECT_PATH)lo_spectralfunction_helpers.o\
	$(OBJECT_PATH)lo_distributed_phonon_dispersion_relations.o\
	$(OBJECT_PATH)lo_selfenergy_interpolation.o
	$(F90) $(F90FLAGS) -c main.f90 $(LIBS) -o $@

$(OBJECT_PATH)lo_distributed_phonon_dispersion_relations.o:\
	lo_distributed_phonon_dispersion_relations.f90
	$(F90) $(F90FLAGS) -c lo_distributed_phonon_dispersion_relations.f90 $(LIBS) -o $@

$(OBJECT_PATH)options.o:\
	options.f90
	$(F90) $(F90FLAGS_LOW) -c options.f90 $(LIBS) -o $@

$(OBJECT_PATH)lo_thermal_transport.o:\
	lo_thermal_transport.f90\
	$(OBJECT_PATH)lo_spectralfunction_helpers.o
	$(F90) $(F90FLAGS) -c lo_thermal_transport.f90 $(LIBS) -o $@

$(OBJECT_PATH)create_selfenergy_interpolation.o:\
	create_selfenergy_interpolation.f90\
	$(OBJECT_PATH)lo_evaluate_phonon_self_energy.o\
	$(OBJECT_PATH)lo_distributed_phonon_dispersion_relations.o\
	$(OBJECT_PATH)lo_spectralfunction_helpers.o
	$(F90) $(F90FLAGS) -c create_selfenergy_interpolation.f90 $(LIBS) -o $@

$(OBJECT_PATH)lo_tetrahedron_interpolation.o:\
	lo_tetrahedron_interpolation.f90
	$(F90) $(F90FLAGS_LOW) -c lo_tetrahedron_interpolation.f90 $(LIBS) -o $@

$(OBJECT_PATH)lo_selfenergy_interpolation.o:\
	lo_selfenergy_interpolation.f90\
	$(OBJECT_PATH)lo_tetrahedron_interpolation.o\
	$(OBJECT_PATH)lo_spectralfunction_helpers.o\
	$(OBJECT_PATH)lo_thermal_transport.o
	$(F90) $(F90FLAGS) -c lo_selfenergy_interpolation.f90 $(LIBS) -o $@
$(OBJECT_PATH)lo_selfenergy_interpolation_evaluate.o:\
	lo_selfenergy_interpolation_evaluate.f90\
	$(OBJECT_PATH)lo_selfenergy_interpolation.o
	$(F90) $(F90FLAGS) -c lo_selfenergy_interpolation_evaluate.f90 $(LIBS) -o $@
$(OBJECT_PATH)lo_selfenergy_interpolation_path.o:\
	lo_selfenergy_interpolation_path.f90\
	$(OBJECT_PATH)lo_selfenergy_interpolation.o
	$(F90) $(F90FLAGS) -c lo_selfenergy_interpolation_path.f90 $(LIBS) -o $@
$(OBJECT_PATH)lo_selfenergy_interpolation_grid.o:\
	lo_selfenergy_interpolation_grid.f90\
	$(OBJECT_PATH)lo_selfenergy_interpolation.o
	$(F90) $(F90FLAGS) -c lo_selfenergy_interpolation_grid.f90 $(LIBS) -o $@

$(OBJECT_PATH)lo_evaluate_phonon_self_energy.o:\
	lo_evaluate_phonon_self_energy.f90\
	$(OBJECT_PATH)options.o\
	$(OBJECT_PATH)lo_selfenergy_interpolation.o\
	$(OBJECT_PATH)lo_spectralfunction_convolution.o\
	$(OBJECT_PATH)lo_spectralfunction_helpers.o\
	$(OBJECT_PATH)lo_distributed_phonon_dispersion_relations.o
	$(F90) $(F90FLAGS) -c lo_evaluate_phonon_self_energy.f90 $(LIBS) -o $@
$(OBJECT_PATH)lo_evaluate_phonon_self_energy_generate.o:\
	lo_evaluate_phonon_self_energy_generate.f90\
	$(OBJECT_PATH)lo_evaluate_phonon_self_energy.o
	$(F90) $(F90FLAGS) -c lo_evaluate_phonon_self_energy_generate.f90 $(LIBS) -o $@
$(OBJECT_PATH)lo_evaluate_phonon_self_energy_gaussian.o:\
	lo_evaluate_phonon_self_energy_gaussian.f90\
	$(OBJECT_PATH)lo_evaluate_phonon_self_energy.o
	$(F90) $(F90FLAGS) -c lo_evaluate_phonon_self_energy_gaussian.f90 $(LIBS) -o $@
$(OBJECT_PATH)lo_evaluate_phonon_self_energy_threephonon.o:\
	lo_evaluate_phonon_self_energy_threephonon.f90\
	$(OBJECT_PATH)lo_evaluate_phonon_self_energy.o
	$(F90) $(F90FLAGS) -c lo_evaluate_phonon_self_energy_threephonon.f90 $(LIBS) -o $@
$(OBJECT_PATH)lo_evaluate_phonon_self_energy_fourphonon.o:\
	lo_evaluate_phonon_self_energy_fourphonon.f90\
	$(OBJECT_PATH)lo_evaluate_phonon_self_energy.o
	$(F90) $(F90FLAGS) -c lo_evaluate_phonon_self_energy_fourphonon.f90 $(LIBS) -o $@
$(OBJECT_PATH)lo_evaluate_phonon_self_energy_helpers.o:\
	lo_evaluate_phonon_self_energy_helpers.f90\
	$(OBJECT_PATH)lo_evaluate_phonon_self_energy.o
	$(F90) $(F90FLAGS) -c lo_evaluate_phonon_self_energy_helpers.f90 $(LIBS) -o $@

$(OBJECT_PATH)lo_collision_matrix.o:\
	lo_collision_matrix.f90\
	$(OBJECT_PATH)lo_selfenergy_interpolation.o\
	$(OBJECT_PATH)lo_evaluate_phonon_self_energy.o\
	$(OBJECT_PATH)lo_spectralfunction_convolution.o\
	$(OBJECT_PATH)lo_spectralfunction_helpers.o
	$(F90) $(F90FLAGS) -c lo_collision_matrix.f90 $(LIBS) -o $@

$(OBJECT_PATH)lo_spectralfunction_convolution.o:\
	lo_spectralfunction_convolution.f90
	$(F90) $(F90FLAGS) -c lo_spectralfunction_convolution.f90 $(LIBS) -o $@

$(OBJECT_PATH)lo_spectralfunction_helpers.o:\
	lo_spectralfunction_helpers.f90
	$(F90) $(F90FLAGS) -c lo_spectralfunction_helpers.f90 $(LIBS) -o $@
