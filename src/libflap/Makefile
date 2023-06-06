#
# Olles makefile for FLAP.
#
# Computer dependent stuff
include Makefile.inc

PROGNAME=libflap

OBJECT_PATH=../../build/libflap/
MODULE_PATH=../../inc/libflap/
DEST_PATH=../../lib/

#LPATH = $(blaslapackLPATH)
#IPATH = $(blaslapackIPATH)
#LIBS = $(blaslapackLIBS)

FFLAGS = $(MODULE_FLAG) $(MODULE_PATH)
#FFLAGS = $(OPT) $(MODULE_FLAG) $(MODULE_PATH)

# Stuff that needs to be compiled
OBJS = \
$(OBJECT_PATH)penf.o \
$(OBJECT_PATH)flap.o \
$(OBJECT_PATH)flap_command_line_argument_t.o \
$(OBJECT_PATH)flap_command_line_arguments_group_t.o \
$(OBJECT_PATH)flap_command_line_interface_t.o \
$(OBJECT_PATH)flap_object_t.o \
$(OBJECT_PATH)flap_utils_m.o

# rules for the library
$(DEST_PATH)libflap.a: $(OBJS)
	-rm $(DEST_PATH)libflap.a
	ar -rcs $(DEST_PATH)libflap.a $(OBJS); ranlib $(DEST_PATH)libflap.a

# I did all the rules manually. Because I hate myself less than I hate make.
$(OBJECT_PATH)penf.o:
	$(FC) $(FFLAGS) -c penf.F90 -o $@ 
$(OBJECT_PATH)flap_object_t.o: $(OBJECT_PATH)penf.o
	$(FC) $(FFLAGS) -c flap_object_t.f90 -o $@ 
$(OBJECT_PATH)flap_utils_m.o: $(OBJECT_PATH)penf.o
	$(FC) $(FFLAGS) -c flap_utils_m.f90 -o $@
$(OBJECT_PATH)flap_command_line_argument_t.o: $(OBJECT_PATH)penf.o $(OBJECT_PATH)flap_object_t.o $(OBJECT_PATH)flap_utils_m.o
	$(FC) $(FFLAGS) -c flap_command_line_argument_t.F90 -o $@
$(OBJECT_PATH)flap_command_line_arguments_group_t.o:\
	$(OBJECT_PATH)flap_command_line_argument_t.o\
	$(OBJECT_PATH)penf.o
	$(FC) $(FFLAGS) -c flap_command_line_arguments_group_t.f90 -o $@
$(OBJECT_PATH)flap_command_line_interface_t.o: $(OBJECT_PATH)penf.o $(OBJECT_PATH)flap_command_line_argument_t.o $(OBJECT_PATH)flap_command_line_arguments_group_t.o
	$(FC) $(FFLAGS) -c flap_command_line_interface_t.F90 -o $@
$(OBJECT_PATH)flap.o: $(OBJECT_PATH)flap_command_line_argument_t.o $(OBJECT_PATH)flap_command_line_arguments_group_t.o $(OBJECT_PATH)flap_command_line_interface_t.o
	$(FC) $(FFLAGS) -c flap.f90 -o $@
#
clean:
	rm -rf $(OBJECT_PATH)*.o $(MODULE_PATH)*.mod
