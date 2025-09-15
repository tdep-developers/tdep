#
# Makefile for Ollelib
#
# Computer dependent stuff
include Makefile.inc
# Where to put things
OBJECT_PATH=../../build/libolle/
MODULE_PATH=../../inc/libolle/
DEST_PATH=../../lib/

LPATH = $(blaslapackLPATH) $(incLPATHhdf) $(incLPATHfft) $(incLPATHmpi)
IPATH = $(blaslapackIPATH) $(incIPATHhdf) $(incIPATHfft) $(incIPATHmpi)
LIBS = $(blaslapackLIBS) $(incLIBShdf) $(incLIBSfft) $(incLIBSmpi)

#OPT = -O0 -fbacktrace -fcheck=all -finit-real=nan -finit-derived -fmax-errors=10 --warn-all # -ffpe-trap=zero,overflow,underflow
#OPT = -O0 -fbacktrace -fcheck=all -finit-real=nan -finit-derived -fmax-errors=10 --pedantic -Wall -Wextra -pedantic -Wcast-align -Wdisabled-optimization -Wmissing-include-dirs -Wshadow -Wunused -fdiagnostics-show-option -Wstrict-overflow=0 -Wno-maybe-uninitialized -Wno-compare-reals -Wno-character-truncation -Werror #-Wrealloc-lhs
#OPT = -O0 -check bounds -check uninit -check pointers -traceback -g -fpe0
#OPT = -traceback -fpe0 #-check bounds -check uninit -check pointers -traceback -g -fpe0


# For debugging/throwing warnings and all that.
#FLG0 = -O0 -fbacktrace -fcheck=all -finit-real=nan -finit-derived -fmax-errors=10 --pedantic -Wall -Wextra -pedantic -Wcast-align -Wdisabled-optimization -Wmissing-include-dirs -Wshadow -Wunused -fdiagnostics-show-option -Wstrict-overflow=0 -Wno-maybe-uninitialized -Werror #-Wrealloc-lhs
#FLG1 = -O0 -fbacktrace -fcheck=all -finit-real=nan -finit-derived -fmax-errors=10 --pedantic -Wall -Wextra -pedantic -Wcast-align -Wdisabled-optimization -Wmissing-include-dirs -Wshadow -Wunused -fdiagnostics-show-option -Wstrict-overflow=0 -Wno-maybe-uninitialized -Wno-compare-reals -Werror #-Wrealloc-lhs

FFLAGS = $(OPT) $(LPATH) $(LIBS) $(IPATH) $(MODULE_FLAG) $(MODULE_PATH) $(FLG0)
FFLAGS1 = $(OPT) $(LPATH) $(LIBS) $(IPATH) $(MODULE_FLAG) $(MODULE_PATH) $(FLG1)
FFLAGS_LOW = $(OPTLOW) $(LPATH) $(LIBS) $(IPATH) $(MODULE_FLAG) $(MODULE_PATH) $(FLG0)

OBJSf = \
$(OBJECT_PATH)konstanter.o \
$(OBJECT_PATH)gottochblandat.o \
$(OBJECT_PATH)gottochblandat_tensors.o \
$(OBJECT_PATH)gottochblandat_sorting.o \
$(OBJECT_PATH)gottochblandat_linalg.o \
$(OBJECT_PATH)gottochblandat_calculus.o \
$(OBJECT_PATH)gottochblandat_physics.o \
$(OBJECT_PATH)gottochblandat_boxes.o \
$(OBJECT_PATH)lo_memtracker.o \
$(OBJECT_PATH)lo_timetracker.o \
$(OBJECT_PATH)lo_sorting.o \
$(OBJECT_PATH)lo_verletboxes.o \
$(OBJECT_PATH)lo_voronoi.o \
$(OBJECT_PATH)lo_randomnumbers.o \
$(OBJECT_PATH)lo_brents_method.o \
$(OBJECT_PATH)lo_distribute_and_classify_things.o \
$(OBJECT_PATH)lo_electron_density_of_states.o \
$(OBJECT_PATH)lo_electron_dispersion_relations.o \
$(OBJECT_PATH)lo_electron_bandstructure_on_path.o \
$(OBJECT_PATH)lo_symmetry_of_interactions.o \
$(OBJECT_PATH)lo_symmetry_of_interactions_helpers.o \
$(OBJECT_PATH)lo_symmetry_of_interactions_tuplets.o \
$(OBJECT_PATH)lo_symmetry_of_interactions_nullspace.o \
$(OBJECT_PATH)lo_longrange_electrostatics.o \
$(OBJECT_PATH)lo_longrange_electrostatics_dynmat.o \
$(OBJECT_PATH)lo_dielectric_interaction.o \
$(OBJECT_PATH)lo_dielectric_interaction_matrixelements.o \
$(OBJECT_PATH)quadratures_stencils.o \
$(OBJECT_PATH)quadratic_programming.o \
$(OBJECT_PATH)quadratic_programming_slvr.o \
$(OBJECT_PATH)geometryfunctions.o \
$(OBJECT_PATH)geometryfunctions_slicingroutines.o \
$(OBJECT_PATH)geometryfunctions_objectprocedures.o \
$(OBJECT_PATH)geometryfunctions_generatingfunctions.o \
$(OBJECT_PATH)dump_data.o \
$(OBJECT_PATH)cgal_wrappers.o \
$(OBJECT_PATH)hdf5_wrappers.o \
$(OBJECT_PATH)fftw_wrappers.o \
$(OBJECT_PATH)mpi_wrappers.o \
$(OBJECT_PATH)glmnet_wrappers.o \
$(OBJECT_PATH)ifc_solvers.o \
$(OBJECT_PATH)ifc_solvers_putget.o \
$(OBJECT_PATH)ifc_solvers_prepsolver.o \
$(OBJECT_PATH)ifc_solvers_lsq.o \
$(OBJECT_PATH)ifc_solvers_diel.o \
$(OBJECT_PATH)ifc_solvers_trainpred.o \
$(OBJECT_PATH)ifc_solvers_diagnostics.o \
$(OBJECT_PATH)ifc_solvers_fastugly.o \
$(OBJECT_PATH)type_distancetable.o \
$(OBJECT_PATH)type_equation_of_state.o \
$(OBJECT_PATH)type_equation_of_state_birchmur.o \
$(OBJECT_PATH)type_equation_of_state_vinet.o \
$(OBJECT_PATH)type_equation_of_state_2Dbirchmur.o \
$(OBJECT_PATH)type_graph_distancetable.o \
$(OBJECT_PATH)type_voronoi.o \
$(OBJECT_PATH)type_voronoi_distancetable.o \
$(OBJECT_PATH)type_blas_lapack_wrappers.o \
$(OBJECT_PATH)type_crystalstructure.o \
$(OBJECT_PATH)type_crystalstructure_atomdata.o \
$(OBJECT_PATH)type_crystalstructure_io.o \
$(OBJECT_PATH)type_crystalstructure_symmetry.o \
$(OBJECT_PATH)type_crystalstructure_alloy.o \
$(OBJECT_PATH)type_embeddedatom.o \
$(OBJECT_PATH)type_forceconstant_alloy_secondorder.o \
$(OBJECT_PATH)type_forceconstant_firstorder.o \
$(OBJECT_PATH)type_forceconstant_firstorder_io.o \
$(OBJECT_PATH)type_forceconstant_secondorder.o \
$(OBJECT_PATH)type_forceconstant_secondorder_aux.o \
$(OBJECT_PATH)type_forceconstant_secondorder_io.o \
$(OBJECT_PATH)type_forceconstant_secondorder_loto.o \
$(OBJECT_PATH)type_forceconstant_secondorder_dynamicalmatrix.o \
$(OBJECT_PATH)type_forceconstant_thirdorder.o \
$(OBJECT_PATH)type_forceconstant_fourthorder.o \
$(OBJECT_PATH)type_jij_secondorder.o \
$(OBJECT_PATH)type_lennardjones.o \
$(OBJECT_PATH)type_linkedlist.o \
$(OBJECT_PATH)type_mdsim.o \
$(OBJECT_PATH)type_qpointmesh.o \
$(OBJECT_PATH)type_qpointmesh_wedgegeneration.o \
$(OBJECT_PATH)type_qpointmesh_gridgeneration.o \
$(OBJECT_PATH)type_qpointmesh_integrationweights.o \
$(OBJECT_PATH)type_qpointmesh_bandstructure.o \
$(OBJECT_PATH)type_qpointmesh_commensurate.o \
$(OBJECT_PATH)type_qpointmesh_io.o \
$(OBJECT_PATH)lo_phonon_bandstructure_on_path.o \
$(OBJECT_PATH)type_phonon_dispersions.o \
$(OBJECT_PATH)type_phonon_dispersions_generation.o \
$(OBJECT_PATH)type_phonon_dos.o \
$(OBJECT_PATH)type_phonon_dos_integrations.o \
$(OBJECT_PATH)type_symmetryoperation.o \
$(OBJECT_PATH)type_effective_hamiltonian.o \
$(OBJECT_PATH)type_forcemap.o \
$(OBJECT_PATH)type_forcemap_generate.o \
$(OBJECT_PATH)type_forcemap_constraints.o \
$(OBJECT_PATH)type_forcemap_coefficient_pair.o \
$(OBJECT_PATH)type_forcemap_coefficient_triplet.o \
$(OBJECT_PATH)type_forcemap_coefficient_quartet.o \
$(OBJECT_PATH)type_forcemap_coefficient_quartet_01.o \
$(OBJECT_PATH)type_forcemap_coefficient_quartet_02.o \
$(OBJECT_PATH)type_forcemap_coefficient_quartet_03.o \
$(OBJECT_PATH)type_forcemap_coefficient_quartet_04.o \
$(OBJECT_PATH)type_forcemap_coefficient_quartet_05.o \
$(OBJECT_PATH)type_forcemap_coefficient_quartet_06.o \
$(OBJECT_PATH)type_forcemap_coefficient_quartet_07.o \
$(OBJECT_PATH)type_forcemap_coefficient_quartet_08.o \
$(OBJECT_PATH)type_forcemap_coefficient_quartet_09.o \
$(OBJECT_PATH)type_forcemap_coefficient_quartet_10.o \
$(OBJECT_PATH)type_forcemap_coefficient_quartet_11.o \
$(OBJECT_PATH)type_forcemap_coefficient_quartet_12.o \
$(OBJECT_PATH)type_forcemap_coefficient_diel.o \
$(OBJECT_PATH)type_forcemap_returntensors.o \
$(OBJECT_PATH)type_forcemap_io.o \
$(OBJECT_PATH)type_lassosolvers.o\
$(OBJECT_PATH)lo_fftgrid_helper.o
# maybe CGAL, maybe not.
ifeq ($(USECGAL),yes)
OBJScg = \
$(OBJECT_PATH)cgal_cleanup.o \
$(OBJECT_PATH)cgal_chull2.o \
$(OBJECT_PATH)cgal_chull3.o \
$(OBJECT_PATH)cgal_chull3_intersection.o \
$(OBJECT_PATH)cgal_deltri2.o \
$(OBJECT_PATH)cgal_deltri3.o \
$(OBJECT_PATH)cgal_tesselate_polyhedron.o
endif

OBJS = $(OBJSf) $(OBJSa) $(OBJScg)

# rules to make the actual library
static: $(OBJS)
	-rm $(DEST_PATH)libolle.a
	ar rc $(DEST_PATH)libolle.a $(OBJS)

UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
    LIB_EXT = dylib
    OSFLAGS = -dynamiclib -install_name @rpath/libolle.$(LIB_EXT)
    LOAD_CMD = -all_load
else ifeq ($(UNAME_S),Linux)
    OSFLAGS = -shared
    LIB_EXT = so
    LOAD_CMD = -Wl,--whole-archive
    END_LOAD_CMD = -Wl,--no-whole-archive
endif

# Get rid of all the other flags on there
COMPILER = $(shell echo $(FC) | cut -d' ' -f1)

shared: $(OBJS)
	-rm $(DEST_PATH)libolle.$(LIB_EXT)
	$(COMPILER) $(OSFLAGS) -o $(DEST_PATH)libolle.$(LIB_EXT) $(LOAD_CMD) $(DEST_PATH)libolle.a $(END_LOAD_CMD) $(LDFLAGS) $(LPATH) $(LIBS)

# Stuff written in C
$(OBJECT_PATH)type_ipi_sockets.o:
	$(CC) -c type_ipi_sockets.c -o $@

ifeq ($(USECGAL),yes)
# c++ stuff, mainly cgal things
$(OBJECT_PATH)cgal_cleanup.o:
	$(CPP) -c cgal_cleanup.cc -o $@
$(OBJECT_PATH)cgal_chull2.o:
	$(CPP) -c cgal_chull2.cc -o $@
$(OBJECT_PATH)cgal_chull3.o:
	$(CPP) -c cgal_chull3.cc -o $@
$(OBJECT_PATH)cgal_chull3_intersection.o:
	$(CPP) -c cgal_chull3_intersection.cc -o $@
$(OBJECT_PATH)cgal_deltri2.o:
	$(CPP) -c cgal_deltri2.cc -o $@
$(OBJECT_PATH)cgal_deltri3.o:
	$(CPP) -c cgal_deltri3.cc -o $@
$(OBJECT_PATH)cgal_tesselate_polyhedron.o:
	$(CPP) -c cgal_tesselate_polyhedron.cc -o $@
endif

# The general stuff, does not depend on too complicated things
$(OBJECT_PATH)konstanter.o:
	$(FC) $(FFLAGS) -c konstanter.f90 -o $@

$(OBJECT_PATH)lo_memtracker.o:\
    $(OBJECT_PATH)konstanter.o
	$(FC) $(FFLAGS) -c lo_memtracker.f90 -o $@

$(OBJECT_PATH)lo_timetracker.o:\
    $(OBJECT_PATH)konstanter.o\
	$(OBJECT_PATH)mpi_wrappers.o\
	$(OBJECT_PATH)gottochblandat.o
	$(FC) $(FFLAGS) -c lo_timetracker.f90 -o $@

$(OBJECT_PATH)lo_sorting.o:\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)lo_memtracker.o
	$(FC) $(FFLAGS) -c lo_sorting.f90 -o $@

$(OBJECT_PATH)lo_voronoi.o:\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)lo_memtracker.o\
	$(OBJECT_PATH)gottochblandat.o\
	$(OBJECT_PATH)geometryfunctions.o\
	$(OBJECT_PATH)type_distancetable.o\
	$(OBJECT_PATH)mpi_wrappers.o
	$(FC) $(FFLAGS) -c lo_voronoi.f90 -o $@

$(OBJECT_PATH)gottochblandat.o:\
	gottochblandat.f90 \
    $(OBJECT_PATH)konstanter.o
	$(FC) $(FFLAGS) -c gottochblandat.f90 -o $@
$(OBJECT_PATH)gottochblandat_tensors.o:\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)konstanter.o
	$(FC) $(FFLAGS) -c gottochblandat_tensors.f90 -o $@
$(OBJECT_PATH)gottochblandat_sorting.o:\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)konstanter.o
	$(FC) $(FFLAGS) -c gottochblandat_sorting.f90 -o $@
$(OBJECT_PATH)gottochblandat_linalg.o:\
	gottochblandat_linalg.f90 \
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)konstanter.o
	$(FC) $(FFLAGS1) -c gottochblandat_linalg.f90 -o $@
$(OBJECT_PATH)gottochblandat_calculus.o:\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)konstanter.o
	$(FC) $(FFLAGS) -c gottochblandat_calculus.f90 -o $@
$(OBJECT_PATH)gottochblandat_physics.o:\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)konstanter.o
	$(FC) $(FFLAGS) -c gottochblandat_physics.f90 -o $@
$(OBJECT_PATH)gottochblandat_boxes.o:\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)konstanter.o
	$(FC) $(FFLAGS) -c gottochblandat_boxes.f90 -o $@

$(OBJECT_PATH)mpi_wrappers.o:\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o
	$(FC) $(FFLAGS) -c mpi_wrappers.f90 -o $@

$(OBJECT_PATH)lo_randomnumbers.o:\
    $(OBJECT_PATH)konstanter.o
	$(FC) $(FFLAGS) -c lo_randomnumbers.f90 -o $@

$(OBJECT_PATH)lo_brents_method.o:\
    $(OBJECT_PATH)konstanter.o
	$(FC) $(FFLAGS) -c lo_brents_method.f90 -o $@

$(OBJECT_PATH)lo_distribute_and_classify_things.o:\
    $(OBJECT_PATH)konstanter.o \
    $(OBJECT_PATH)mpi_wrappers.o \
    $(OBJECT_PATH)lo_memtracker.o \
    $(OBJECT_PATH)lo_sorting.o \
    $(OBJECT_PATH)gottochblandat.o
	$(FC) $(FFLAGS) -c lo_distribute_and_classify_things.f90 -o $@

$(OBJECT_PATH)lo_symmetry_of_interactions.o:\
    lo_symmetry_of_interactions.f90 \
    $(OBJECT_PATH)konstanter.o \
    $(OBJECT_PATH)mpi_wrappers.o \
    $(OBJECT_PATH)lo_memtracker.o \
	$(OBJECT_PATH)type_crystalstructure.o \
	$(OBJECT_PATH)type_distancetable.o \
	$(OBJECT_PATH)type_voronoi_distancetable.o \
    $(OBJECT_PATH)gottochblandat.o
	$(FC) $(FFLAGS) -c lo_symmetry_of_interactions.f90 -o $@
$(OBJECT_PATH)lo_symmetry_of_interactions_helpers.o:\
	lo_symmetry_of_interactions_helpers.f90 \
    $(OBJECT_PATH)lo_symmetry_of_interactions.o \
	$(OBJECT_PATH)type_blas_lapack_wrappers.o
	$(FC) $(FFLAGS) -c lo_symmetry_of_interactions_helpers.f90 -o $@
$(OBJECT_PATH)lo_symmetry_of_interactions_tuplets.o:\
    $(OBJECT_PATH)lo_symmetry_of_interactions.o \
	$(OBJECT_PATH)lo_sorting.o
	$(FC) $(FFLAGS) -c lo_symmetry_of_interactions_tuplets.f90 -o $@
$(OBJECT_PATH)lo_symmetry_of_interactions_nullspace.o:\
    lo_symmetry_of_interactions_nullspace.f90\
    $(OBJECT_PATH)lo_symmetry_of_interactions.o
	$(FC) $(FFLAGS) -c lo_symmetry_of_interactions_nullspace.f90 -o $@

$(OBJECT_PATH)lo_longrange_electrostatics.o:\
	lo_longrange_electrostatics.f90 \
    $(OBJECT_PATH)konstanter.o \
    $(OBJECT_PATH)mpi_wrappers.o \
    $(OBJECT_PATH)lo_memtracker.o \
	$(OBJECT_PATH)lo_brents_method.o \
	$(OBJECT_PATH)type_crystalstructure.o \
    $(OBJECT_PATH)gottochblandat.o
	$(FC) $(FFLAGS) -c lo_longrange_electrostatics.f90 -o $@
$(OBJECT_PATH)lo_longrange_electrostatics_dynmat.o:\
	lo_longrange_electrostatics_dynmat.f90 \
    $(OBJECT_PATH)lo_longrange_electrostatics.o
	$(FC) $(FFLAGS) -c lo_longrange_electrostatics_dynmat.f90 -o $@

$(OBJECT_PATH)lo_dielectric_interaction.o:\
    $(OBJECT_PATH)konstanter.o \
    $(OBJECT_PATH)mpi_wrappers.o \
	$(OBJECT_PATH)gottochblandat.o \
	$(OBJECT_PATH)type_crystalstructure.o
	$(FC) $(FFLAGS) -c lo_dielectric_interaction.f90 -o $@
$(OBJECT_PATH)lo_dielectric_interaction_matrixelements.o:\
	$(OBJECT_PATH)lo_dielectric_interaction.o
	$(FC) $(FFLAGS) -c lo_dielectric_interaction_matrixelements.f90 -o $@

$(OBJECT_PATH)quadratures_stencils.o:\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o
	$(FC) $(FFLAGS) -c quadratures_stencils.f90 -o $@

$(OBJECT_PATH)type_linkedlist.o:\
    $(OBJECT_PATH)konstanter.o
	$(FC) $(FFLAGS_LOW) -c type_linkedlist.f90 -o $@

$(OBJECT_PATH)geometryfunctions.o: \
	geometryfunctions.f90\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o
	$(FC) $(FFLAGS) -c geometryfunctions.f90 -o $@
$(OBJECT_PATH)geometryfunctions_slicingroutines.o: \
    $(OBJECT_PATH)geometryfunctions.o\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o
	$(FC) $(FFLAGS_LOW) -c geometryfunctions_slicingroutines.f90 -o $@
$(OBJECT_PATH)geometryfunctions_objectprocedures.o: \
    $(OBJECT_PATH)geometryfunctions.o\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o
	$(FC) $(FFLAGS) -c geometryfunctions_objectprocedures.f90 -o $@
$(OBJECT_PATH)geometryfunctions_generatingfunctions.o: \
    $(OBJECT_PATH)geometryfunctions.o\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o
	$(FC) $(FFLAGS) -c geometryfunctions_generatingfunctions.f90 -o $@

$(OBJECT_PATH)dump_data.o:\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o
	$(FC) $(FFLAGS) -c dump_data.f90 -o $@

$(OBJECT_PATH)hdf5_wrappers.o:\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)konstanter.o
	$(FC) $(FFLAGS) -O0 -c hdf5_wrappers.f90 -o $@

$(OBJECT_PATH)fftw_wrappers.o:\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o
	$(FC) $(FFLAGS) -c fftw_wrappers.f90 -o $@

$(OBJECT_PATH)lo_verletboxes.o:\
    $(OBJECT_PATH)konstanter.o \
    $(OBJECT_PATH)lo_memtracker.o \
    $(OBJECT_PATH)mpi_wrappers.o
	$(FC) $(FFLAGS) -c lo_verletboxes.f90 -o $@

$(OBJECT_PATH)type_blas_lapack_wrappers.o:\
    $(OBJECT_PATH)konstanter.o
	$(FC) $(FFLAGS) -c type_blas_lapack_wrappers.f90 -o $@

$(OBJECT_PATH)glmnet_wrappers.o:\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)mpi_wrappers.o\
    $(OBJECT_PATH)type_blas_lapack_wrappers.o
	$(FC) $(FFLAGS1) -c glmnet_wrappers.f90 -o $@

$(OBJECT_PATH)quadratic_programming.o:\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)mpi_wrappers.o\
    $(OBJECT_PATH)type_blas_lapack_wrappers.o
	$(FC) $(FFLAGS) -c quadratic_programming.f90 -o $@
$(OBJECT_PATH)quadratic_programming_slvr.o:\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)quadratic_programming.o\
    $(OBJECT_PATH)type_blas_lapack_wrappers.o
	$(FC) $(FFLAGS1) -c quadratic_programming_slvr.f90 -o $@

$(OBJECT_PATH)type_lassosolvers.o:\
    $(OBJECT_PATH)konstanter.o \
    $(OBJECT_PATH)gottochblandat.o \
    $(OBJECT_PATH)type_blas_lapack_wrappers.o
	$(FC) $(FFLAGS) -c type_lassosolvers.f90 -o $@
$(OBJECT_PATH)type_distancetable.o:\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)type_blas_lapack_wrappers.o\
    $(OBJECT_PATH)geometryfunctions.o\
    $(OBJECT_PATH)mpi_wrappers.o
	$(FC) $(FFLAGS) -c type_distancetable.f90 -o $@

$(OBJECT_PATH)type_equation_of_state.o:\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)quadratures_stencils.o\
    $(OBJECT_PATH)type_blas_lapack_wrappers.o
	$(FC) $(FFLAGS) -c type_equation_of_state.f90 -o $@
$(OBJECT_PATH)type_equation_of_state_birchmur.o:\
    $(OBJECT_PATH)type_equation_of_state.o\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)quadratures_stencils.o\
    $(OBJECT_PATH)type_blas_lapack_wrappers.o
	$(FC) $(FFLAGS) -c type_equation_of_state_birchmur.f90 -o $@
$(OBJECT_PATH)type_equation_of_state_vinet.o:\
    $(OBJECT_PATH)type_equation_of_state.o\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)quadratures_stencils.o\
    $(OBJECT_PATH)type_blas_lapack_wrappers.o
	$(FC) $(FFLAGS) -c type_equation_of_state_vinet.f90 -o $@
$(OBJECT_PATH)type_equation_of_state_2Dbirchmur.o:\
    $(OBJECT_PATH)type_equation_of_state.o\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)quadratures_stencils.o\
    $(OBJECT_PATH)type_blas_lapack_wrappers.o
	$(FC) $(FFLAGS) -c type_equation_of_state_2Dbirchmur.f90 -o $@

$(OBJECT_PATH)type_voronoi.o:\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)geometryfunctions.o\
    $(OBJECT_PATH)type_distancetable.o\
    $(OBJECT_PATH)type_linkedlist.o
	$(FC) $(FFLAGS_LOW) -c type_voronoi.f90 -o $@
$(OBJECT_PATH)type_graph_distancetable.o:\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)geometryfunctions.o\
    $(OBJECT_PATH)type_voronoi.o\
    $(OBJECT_PATH)type_distancetable.o
	$(FC) $(FFLAGS) -c type_graph_distancetable.f90 -o $@
$(OBJECT_PATH)type_voronoi_distancetable.o:\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)geometryfunctions.o\
    $(OBJECT_PATH)type_distancetable.o\
    $(OBJECT_PATH)type_voronoi.o
	$(FC) $(FFLAGS) -c type_voronoi_distancetable.f90 -o $@
$(OBJECT_PATH)type_symmetryoperation.o:\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)geometryfunctions.o
	$(FC) $(FFLAGS) -c type_symmetryoperation.f90 -o $@
#
# If I have CGAL installed
#
ifeq ($(USECGAL),yes)
$(OBJECT_PATH)cgal_wrappers.o:\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)cgal_cleanup.o\
    $(OBJECT_PATH)cgal_chull2.o\
    $(OBJECT_PATH)cgal_chull3.o\
    $(OBJECT_PATH)cgal_deltri2.o\
    $(OBJECT_PATH)cgal_deltri3.o\
    $(OBJECT_PATH)cgal_tesselate_polyhedron.o
	$(FC) $(FFLAGS) -c cgal_wrappers.f90 -o $@
else
$(OBJECT_PATH)cgal_wrappers.o:\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o
	$(FC) $(FFLAGS) -c cgal_wrappers.f90 -o $@
endif
#
# Sockets for I-PI
#
#$(OBJECT_PATH)type_ipi.o:\
#    $(OBJECT_PATH)type_ipi_sockets.o
#	$(FC) $(FFLAGS) -c type_ipi.f90 -o $@
#
# A bit more advanced
#
$(OBJECT_PATH)type_crystalstructure.o:\
	type_crystalstructure.f90\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)geometryfunctions.o\
    $(OBJECT_PATH)type_voronoi.o\
    $(OBJECT_PATH)type_distancetable.o\
    $(OBJECT_PATH)type_symmetryoperation.o
	$(FC) $(FFLAGS) -c type_crystalstructure.f90 -o $@
$(OBJECT_PATH)type_crystalstructure_atomdata.o:\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)type_crystalstructure.o
	$(FC) $(FFLAGS) -c type_crystalstructure_atomdata.f90 -o $@
$(OBJECT_PATH)type_crystalstructure_io.o:\
	type_crystalstructure_io.f90\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)hdf5_wrappers.o\
    $(OBJECT_PATH)type_crystalstructure.o
	$(FC) $(FFLAGS) -c type_crystalstructure_io.f90 -o $@
$(OBJECT_PATH)type_crystalstructure_symmetry.o:\
	type_crystalstructure_symmetry.f90\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)type_crystalstructure.o
	$(FC) $(FFLAGS) -c type_crystalstructure_symmetry.f90 -o $@
$(OBJECT_PATH)type_crystalstructure_alloy.o:\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)type_crystalstructure.o
	$(FC) $(FFLAGS) -c type_crystalstructure_alloy.f90 -o $@

$(OBJECT_PATH)type_forceconstant_firstorder.o:\
	type_forceconstant_firstorder.f90\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)type_crystalstructure.o
	$(FC) $(FFLAGS) -c type_forceconstant_firstorder.f90 -o $@
$(OBJECT_PATH)type_forceconstant_firstorder_io.o:\
	type_forceconstant_firstorder_io.f90\
    $(OBJECT_PATH)type_forceconstant_firstorder.o
	$(FC) $(FFLAGS) -c type_forceconstant_firstorder_io.f90 -o $@

$(OBJECT_PATH)type_forceconstant_secondorder.o:\
	type_forceconstant_secondorder.f90\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)mpi_wrappers.o\
	$(OBJECT_PATH)lo_memtracker.o\
	$(OBJECT_PATH)lo_longrange_electrostatics.o\
    $(OBJECT_PATH)type_crystalstructure.o\
    $(OBJECT_PATH)type_qpointmesh.o
	$(FC) $(FFLAGS) -c type_forceconstant_secondorder.f90 -o $@
$(OBJECT_PATH)type_forceconstant_secondorder_dynamicalmatrix.o:\
	type_forceconstant_secondorder_dynamicalmatrix.f90\
    $(OBJECT_PATH)type_forceconstant_secondorder.o\
	$(OBJECT_PATH)lo_sorting.o\
    $(OBJECT_PATH)type_blas_lapack_wrappers.o
	$(FC) $(FFLAGS) -c type_forceconstant_secondorder_dynamicalmatrix.f90 -o $@
$(OBJECT_PATH)type_forceconstant_secondorder_aux.o:\
	type_forceconstant_secondorder_aux.f90\
    $(OBJECT_PATH)type_forceconstant_secondorder.o\
    $(OBJECT_PATH)type_distancetable.o\
    $(OBJECT_PATH)lo_randomnumbers.o\
	$(OBJECT_PATH)lo_sorting.o\
    $(OBJECT_PATH)type_blas_lapack_wrappers.o
	$(FC) $(FFLAGS) -c type_forceconstant_secondorder_aux.f90 -o $@
$(OBJECT_PATH)type_forceconstant_secondorder_io.o:\
	type_forceconstant_secondorder_io.f90\
    $(OBJECT_PATH)type_forceconstant_secondorder.o
	$(FC) $(FFLAGS) -c type_forceconstant_secondorder_io.f90 -o $@
$(OBJECT_PATH)type_forceconstant_secondorder_loto.o:\
	type_forceconstant_secondorder_loto.f90\
    $(OBJECT_PATH)type_forceconstant_secondorder.o
	$(FC) $(FFLAGS) -c type_forceconstant_secondorder_loto.f90 -o $@

$(OBJECT_PATH)type_forceconstant_thirdorder.o:\
    type_forceconstant_thirdorder.f90\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)type_crystalstructure.o\
    $(OBJECT_PATH)type_qpointmesh.o\
    $(OBJECT_PATH)type_distancetable.o
	$(FC) $(FFLAGS) -c type_forceconstant_thirdorder.f90 -o $@
$(OBJECT_PATH)type_forceconstant_fourthorder.o:\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)type_crystalstructure.o\
    $(OBJECT_PATH)type_distancetable.o
	$(FC) $(FFLAGS) -c type_forceconstant_fourthorder.f90 -o $@
$(OBJECT_PATH)type_jij_secondorder.o:\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)type_distancetable.o\
    $(OBJECT_PATH)type_crystalstructure.o
	$(FC) $(FFLAGS) -c type_jij_secondorder.f90 -o $@

$(OBJECT_PATH)type_mdsim.o:\
    type_mdsim.f90\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)lo_randomnumbers.o\
    $(OBJECT_PATH)hdf5_wrappers.o\
    $(OBJECT_PATH)mpi_wrappers.o\
    $(OBJECT_PATH)type_crystalstructure.o
	$(FC) $(FFLAGS) -c type_mdsim.f90 -o $@

$(OBJECT_PATH)type_qpointmesh.o:\
	type_qpointmesh.f90\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)mpi_wrappers.o\
    $(OBJECT_PATH)lo_memtracker.o\
    $(OBJECT_PATH)gottochblandat.o\
	$(OBJECT_PATH)geometryfunctions.o\
    $(OBJECT_PATH)type_crystalstructure.o\
    $(OBJECT_PATH)type_symmetryoperation.o
	$(FC) $(FFLAGS) -c type_qpointmesh.f90 -o $@
$(OBJECT_PATH)type_qpointmesh_wedgegeneration.o:\
    $(OBJECT_PATH)cgal_wrappers.o\
    $(OBJECT_PATH)lo_sorting.o\
    $(OBJECT_PATH)lo_verletboxes.o\
    $(OBJECT_PATH)type_qpointmesh.o
	$(FC) $(FFLAGS) -c type_qpointmesh_wedgegeneration.f90 -o $@
$(OBJECT_PATH)type_qpointmesh_gridgeneration.o:\
    type_qpointmesh_gridgeneration.f90\
    $(OBJECT_PATH)type_qpointmesh.o
	$(FC) $(FFLAGS) -c type_qpointmesh_gridgeneration.f90 -o $@
$(OBJECT_PATH)type_qpointmesh_integrationweights.o:\
    $(OBJECT_PATH)lo_voronoi.o\
	$(OBJECT_PATH)lo_sorting.o\
	$(OBJECT_PATH)lo_verletboxes.o\
	$(OBJECT_PATH)type_distancetable.o\
    $(OBJECT_PATH)type_qpointmesh.o
	$(FC) $(FFLAGS) -c type_qpointmesh_integrationweights.f90 -o $@
$(OBJECT_PATH)type_qpointmesh_bandstructure.o:\
    type_qpointmesh_bandstructure.f90\
	$(OBJECT_PATH)lo_sorting.o\
    $(OBJECT_PATH)type_qpointmesh.o
	$(FC) $(FFLAGS) -c type_qpointmesh_bandstructure.f90 -o $@
$(OBJECT_PATH)type_qpointmesh_commensurate.o:\
    $(OBJECT_PATH)type_qpointmesh.o
	$(FC) $(FFLAGS) -c type_qpointmesh_commensurate.f90 -o $@
$(OBJECT_PATH)type_qpointmesh_io.o:\
	type_qpointmesh_io.f90\
    $(OBJECT_PATH)hdf5_wrappers.o\
    $(OBJECT_PATH)type_qpointmesh.o
	$(FC) $(FFLAGS_LOW) -c type_qpointmesh_io.f90 -o $@

$(OBJECT_PATH)type_lennardjones.o:\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)type_crystalstructure.o\
    $(OBJECT_PATH)gottochblandat.o
	$(FC) $(FFLAGS) -c type_lennardjones.f90 -o $@
$(OBJECT_PATH)type_embeddedatom.o:\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)type_crystalstructure.o
	$(FC) $(FFLAGS) -c type_embeddedatom.f90 -o $@

$(OBJECT_PATH)type_phonon_dispersions.o:\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)mpi_wrappers.o\
	$(OBJECT_PATH)lo_memtracker.o\
	$(OBJECT_PATH)lo_sorting.o\
	$(OBJECT_PATH)type_blas_lapack_wrappers.o\
    $(OBJECT_PATH)type_crystalstructure.o\
    $(OBJECT_PATH)type_forceconstant_secondorder.o\
    $(OBJECT_PATH)type_forceconstant_thirdorder.o\
    $(OBJECT_PATH)type_qpointmesh.o\
	$(OBJECT_PATH)type_symmetryoperation.o\
    $(OBJECT_PATH)hdf5_wrappers.o
	$(FC) $(FFLAGS) -c type_phonon_dispersions.f90 -o $@
$(OBJECT_PATH)type_phonon_dispersions_generation.o:\
	type_phonon_dispersions_generation.f90\
    $(OBJECT_PATH)type_phonon_dispersions.o
	$(FC) $(FFLAGS) -c type_phonon_dispersions_generation.f90 -o $@

$(OBJECT_PATH)type_phonon_dos.o:\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)hdf5_wrappers.o\
    $(OBJECT_PATH)mpi_wrappers.o\
    $(OBJECT_PATH)type_crystalstructure.o\
    $(OBJECT_PATH)type_phonon_dispersions.o\
    $(OBJECT_PATH)type_qpointmesh.o\
    $(OBJECT_PATH)type_forceconstant_secondorder.o
	$(FC) $(FFLAGS) -c type_phonon_dos.f90 -o $@
$(OBJECT_PATH)type_phonon_dos_integrations.o:\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)mpi_wrappers.o\
    $(OBJECT_PATH)type_crystalstructure.o\
    $(OBJECT_PATH)type_phonon_dispersions.o\
    $(OBJECT_PATH)type_qpointmesh.o\
    $(OBJECT_PATH)type_phonon_dos.o
	$(FC) $(FFLAGS) -c type_phonon_dos_integrations.f90 -o $@

$(OBJECT_PATH)lo_phonon_bandstructure_on_path.o:\
    $(OBJECT_PATH)konstanter.o\
	$(OBJECT_PATH)gottochblandat.o\
	$(OBJECT_PATH)mpi_wrappers.o\
	$(OBJECT_PATH)lo_memtracker.o\
    $(OBJECT_PATH)hdf5_wrappers.o\
    $(OBJECT_PATH)type_crystalstructure.o\
    $(OBJECT_PATH)type_phonon_dispersions.o\
	$(OBJECT_PATH)type_forceconstant_secondorder.o\
	$(OBJECT_PATH)type_forceconstant_thirdorder.o\
	$(OBJECT_PATH)type_qpointmesh.o
	$(FC) $(FFLAGS) -c lo_phonon_bandstructure_on_path.f90 -o $@

$(OBJECT_PATH)type_forceconstant_alloy_secondorder.o:\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)type_crystalstructure.o\
    $(OBJECT_PATH)type_forceconstant_secondorder.o
	$(FC) $(FFLAGS) -c type_forceconstant_alloy_secondorder.f90 -o $@
# Electron stuff
$(OBJECT_PATH)type_effective_hamiltonian.o:\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)mpi_wrappers.o\
    $(OBJECT_PATH)hdf5_wrappers.o\
    $(OBJECT_PATH)type_crystalstructure.o\
    $(OBJECT_PATH)type_blas_lapack_wrappers.o\
    $(OBJECT_PATH)type_qpointmesh.o\
    $(OBJECT_PATH)type_symmetryoperation.o
	$(FC) $(FFLAGS) -c type_effective_hamiltonian.f90 -o $@
$(OBJECT_PATH)lo_electron_dispersion_relations.o:\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)lo_memtracker.o\
    $(OBJECT_PATH)lo_brents_method.o\
    $(OBJECT_PATH)mpi_wrappers.o\
	$(OBJECT_PATH)hdf5_wrappers.o\
    $(OBJECT_PATH)type_crystalstructure.o\
    $(OBJECT_PATH)type_effective_hamiltonian.o\
    $(OBJECT_PATH)type_qpointmesh.o
	$(FC) $(FFLAGS) -c lo_electron_dispersion_relations.f90 -o $@
$(OBJECT_PATH)lo_electron_density_of_states.o:\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)lo_memtracker.o\
    $(OBJECT_PATH)mpi_wrappers.o\
	$(OBJECT_PATH)hdf5_wrappers.o\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)type_crystalstructure.o\
    $(OBJECT_PATH)lo_electron_dispersion_relations.o\
    $(OBJECT_PATH)type_qpointmesh.o
	$(FC) $(FFLAGS) -c lo_electron_density_of_states.f90 -o $@
$(OBJECT_PATH)lo_electron_bandstructure_on_path.o:\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)lo_memtracker.o\
    $(OBJECT_PATH)mpi_wrappers.o\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)type_crystalstructure.o\
    $(OBJECT_PATH)lo_electron_dispersion_relations.o\
    $(OBJECT_PATH)type_effective_hamiltonian.o\
    $(OBJECT_PATH)type_qpointmesh.o
	$(FC) $(FFLAGS) -c lo_electron_bandstructure_on_path.f90 -o $@

$(OBJECT_PATH)type_forcemap.o:\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o\
	$(OBJECT_PATH)mpi_wrappers.o\
	$(OBJECT_PATH)lo_memtracker.o\
    $(OBJECT_PATH)type_crystalstructure.o\
    $(OBJECT_PATH)type_blas_lapack_wrappers.o\
    $(OBJECT_PATH)type_forceconstant_firstorder.o\
    $(OBJECT_PATH)type_forceconstant_secondorder.o\
    $(OBJECT_PATH)type_forceconstant_thirdorder.o\
    $(OBJECT_PATH)type_forceconstant_fourthorder.o\
	$(OBJECT_PATH)lo_dielectric_interaction.o\
    $(OBJECT_PATH)type_jij_secondorder.o\
	$(OBJECT_PATH)lo_symmetry_of_interactions.o
	$(FC) $(FFLAGS) -c type_forcemap.f90 -o $@
$(OBJECT_PATH)type_forcemap_generate.o:\
    $(OBJECT_PATH)type_forcemap.o
	$(FC) $(FFLAGS) -c type_forcemap_generate.f90 -o $@
$(OBJECT_PATH)type_forcemap_constraints.o:\
    $(OBJECT_PATH)type_forcemap.o
	$(FC) $(FFLAGS) -c type_forcemap_constraints.f90 -o $@
$(OBJECT_PATH)type_forcemap_coefficient_pair.o:\
    $(OBJECT_PATH)type_forcemap.o
	$(FC) $(FFLAGS) -c type_forcemap_coefficient_pair.f90 -o $@
$(OBJECT_PATH)type_forcemap_coefficient_triplet.o:\
    $(OBJECT_PATH)type_forcemap.o
	$(FC) $(FFLAGS) -c type_forcemap_coefficient_triplet.f90 -o $@
$(OBJECT_PATH)type_forcemap_coefficient_quartet.o:\
    $(OBJECT_PATH)type_forcemap.o
	$(FC) $(FFLAGS) -c type_forcemap_coefficient_quartet.f90 -o $@
$(OBJECT_PATH)type_forcemap_coefficient_diel.o:\
    $(OBJECT_PATH)type_forcemap.o
	$(FC) $(FFLAGS) -c type_forcemap_coefficient_diel.f90 -o $@
$(OBJECT_PATH)type_forcemap_returntensors.o:\
    $(OBJECT_PATH)type_forcemap.o
	$(FC) $(FFLAGS) -c type_forcemap_returntensors.f90 -o $@
$(OBJECT_PATH)type_forcemap_io.o:\
    $(OBJECT_PATH)hdf5_wrappers.o\
    $(OBJECT_PATH)type_forcemap.o
	$(FC) $(FFLAGS) -c type_forcemap_io.f90 -o $@
$(OBJECT_PATH)type_forcemap_coefficient_quartet_01.o:\
    $(OBJECT_PATH)type_forcemap.o
	$(FC) $(FFLAGS) -c type_forcemap_coefficient_quartet_01.f90 -o $@
$(OBJECT_PATH)type_forcemap_coefficient_quartet_02.o:\
    $(OBJECT_PATH)type_forcemap.o
	$(FC) $(FFLAGS) -c type_forcemap_coefficient_quartet_02.f90 -o $@
$(OBJECT_PATH)type_forcemap_coefficient_quartet_03.o:\
    $(OBJECT_PATH)type_forcemap.o
	$(FC) $(FFLAGS) -c type_forcemap_coefficient_quartet_03.f90 -o $@
$(OBJECT_PATH)type_forcemap_coefficient_quartet_04.o:\
    $(OBJECT_PATH)type_forcemap.o
	$(FC) $(FFLAGS) -c type_forcemap_coefficient_quartet_04.f90 -o $@
$(OBJECT_PATH)type_forcemap_coefficient_quartet_05.o:\
    $(OBJECT_PATH)type_forcemap.o
	$(FC) $(FFLAGS) -c type_forcemap_coefficient_quartet_05.f90 -o $@
$(OBJECT_PATH)type_forcemap_coefficient_quartet_06.o:\
    $(OBJECT_PATH)type_forcemap.o
	$(FC) $(FFLAGS) -c type_forcemap_coefficient_quartet_06.f90 -o $@
$(OBJECT_PATH)type_forcemap_coefficient_quartet_07.o:\
    $(OBJECT_PATH)type_forcemap.o
	$(FC) $(FFLAGS) -c type_forcemap_coefficient_quartet_07.f90 -o $@
$(OBJECT_PATH)type_forcemap_coefficient_quartet_08.o:\
    $(OBJECT_PATH)type_forcemap.o
	$(FC) $(FFLAGS) -c type_forcemap_coefficient_quartet_08.f90 -o $@
$(OBJECT_PATH)type_forcemap_coefficient_quartet_09.o:\
    $(OBJECT_PATH)type_forcemap.o
	$(FC) $(FFLAGS) -c type_forcemap_coefficient_quartet_09.f90 -o $@
$(OBJECT_PATH)type_forcemap_coefficient_quartet_10.o:\
    $(OBJECT_PATH)type_forcemap.o
	$(FC) $(FFLAGS) -c type_forcemap_coefficient_quartet_10.f90 -o $@
$(OBJECT_PATH)type_forcemap_coefficient_quartet_11.o:\
    $(OBJECT_PATH)type_forcemap.o
	$(FC) $(FFLAGS) -c type_forcemap_coefficient_quartet_11.f90 -o $@
$(OBJECT_PATH)type_forcemap_coefficient_quartet_12.o:\
    $(OBJECT_PATH)type_forcemap.o
	$(FC) $(FFLAGS) -c type_forcemap_coefficient_quartet_12.f90 -o $@

$(OBJECT_PATH)ifc_solvers.o:\
	ifc_solvers.f90 \
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)geometryfunctions.o\
    $(OBJECT_PATH)mpi_wrappers.o\
    $(OBJECT_PATH)hdf5_wrappers.o\
    $(OBJECT_PATH)glmnet_wrappers.o\
    $(OBJECT_PATH)quadratic_programming.o\
    $(OBJECT_PATH)type_crystalstructure.o\
    $(OBJECT_PATH)type_symmetryoperation.o\
    $(OBJECT_PATH)type_qpointmesh.o\
    $(OBJECT_PATH)type_forceconstant_firstorder.o\
    $(OBJECT_PATH)type_forceconstant_secondorder.o\
    $(OBJECT_PATH)type_forceconstant_thirdorder.o\
    $(OBJECT_PATH)type_forceconstant_fourthorder.o\
    $(OBJECT_PATH)type_forcemap.o\
    $(OBJECT_PATH)type_mdsim.o\
    $(OBJECT_PATH)type_distancetable.o\
    $(OBJECT_PATH)type_blas_lapack_wrappers.o
	$(FC) $(FFLAGS) -c ifc_solvers.f90 -o $@
$(OBJECT_PATH)ifc_solvers_trainpred.o:\
	ifc_solvers_trainpred.f90 \
    $(OBJECT_PATH)ifc_solvers.o
	$(FC) $(FFLAGS) -c ifc_solvers_trainpred.f90 -o $@
$(OBJECT_PATH)ifc_solvers_putget.o:\
	ifc_solvers_putget.f90 \
    $(OBJECT_PATH)ifc_solvers.o
	$(FC) $(FFLAGS) -c ifc_solvers_putget.f90 -o $@
$(OBJECT_PATH)ifc_solvers_prepsolver.o:\
	ifc_solvers_prepsolver.f90 \
    $(OBJECT_PATH)ifc_solvers.o
	$(FC) $(FFLAGS) -c ifc_solvers_prepsolver.f90 -o $@
$(OBJECT_PATH)ifc_solvers_lsq.o:\
    $(OBJECT_PATH)ifc_solvers.o
	$(FC) $(FFLAGS) -c ifc_solvers_lsq.f90 -o $@
$(OBJECT_PATH)ifc_solvers_diel.o:\
	ifc_solvers_diel.f90\
    $(OBJECT_PATH)ifc_solvers.o\
	$(OBJECT_PATH)type_blas_lapack_wrappers.o\
	$(OBJECT_PATH)lo_longrange_electrostatics.o
	$(FC) $(FFLAGS) -c ifc_solvers_diel.f90 -o $@
$(OBJECT_PATH)ifc_solvers_diagnostics.o:\
	ifc_solvers_diagnostics.f90 \
    $(OBJECT_PATH)ifc_solvers.o
	$(FC) $(FFLAGS) -c ifc_solvers_diagnostics.f90 -o $@
$(OBJECT_PATH)ifc_solvers_fastugly.o:\
	ifc_solvers_fastugly.f90 \
    $(OBJECT_PATH)ifc_solvers.o
	$(FC) $(FFLAGS) -c ifc_solvers_fastugly.f90 -o $@

$(OBJECT_PATH)lo_fftgrid_helper.o:\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)lo_randomnumbers.o
	$(FC) $(FFLAGS) -c lo_fftgrid_helper.f90 -o $@

clean:
	rm -rf $(OBJECT_PATH)*.o $(MODULE_PATH)*mod
