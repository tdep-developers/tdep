#!/bin/bash
#
# My feeble attempt at a build script. I really should learn cmake or something, but it's
# so boring. The only way this can get improved is if I blackmail someone into fixing it.
#
# it's probably a good idea to abort this script right away if something goes wrong.
set -e

#should everything be cleaned? I keep this in place since I might want more options in the future.
CLEAN=YES
MANPAGE=YES
NTHREADS_MAKE=1
for i in "$@"
do
case $i in
    --noclean)
    CLEAN=NO
    shift # past argument with no value
    ;;
    --clean)
    CLEAN=YES
    shift # past argument with no value
    ;;
    --nomanpage)
    MANPAGE=NO
    shift # past argument with no value
    ;;
    --nthreads_make)
    shift
    NTHREADS_MAKE=$1
# would be cleaner if we checked this was an integer
    shift
    ;;
    *)
            # unknown option
    ;;
esac
done
echo "clean everything? ${CLEAN}"
echo "number of threads: ${NTHREADS_MAKE}"

# Grab wich branch we are on, and which revision
gitbranch=`git rev-parse --abbrev-ref HEAD`
gitrevision=`git rev-parse HEAD`

# Make sure there are directories available. I need them to put stuff in them.
important_directories="
bin
inc
lib
doc
man
man/man1
build"

echo "... checking that I succeded in making the directories I need."
for i in $important_directories
do
    [ -d "${i}" ] || mkdir ${i}
    [ -d "${i}" ] && echo "found ${i}, good"
done

# read the machine-specific settings
source important_settings
echo "parsed the important settings"

# ok, that's a decent start. Start by building the main library, that's the tricky part
cd src/libolle

# make the git information accessible so that it can be added during compilation
cat>gitinformation<<EOF
#define gitbranch '${gitbranch}'
#define gitrevision '${gitrevision}'
EOF

# will I use cgal?
if [ ${USECGAL} == "yes" ]
then
    echo "Ok, will try to include the CGAL stuff"
    PRECOMPILER_FLAGS="${PRECOMPILER_FLAGS} -Dusecgal"
fi
# paste a file with all the Machine-specific info
cat>Makefile.inc<<EOF
# Specify fortran compiler
FC = ${FORTRAN_COMPILER} ${FCFLAGS} ${FCFLAGS_EXTRA} ${PRECOMPILER_FLAGS}
F77C = ${FORTRAN_COMPILER}

# And the flag that tells the compiler where to put modules
MODULE_FLAG = ${MODULE_FLAG}

# Specify optimization
OPT = ${OPTIMIZATION_LEVEL}
OPTLOW = ${OPTIMIZATION_SENSITIVE}

# BLAS/LAPACK paths
blaslapackLPATH = ${PATH_TO_BLASLAPACK_LIB}
blaslapackIPATH = ${PATH_TO_BLASLAPACK_INC}
blaslapackLIBS = ${BLASLAPACK_LIBS}
# MPI
incLPATHmpi = ${PATH_TO_MPI_LIB}
incIPATHmpi = ${PATH_TO_MPI_INC}
incLIBSmpi = ${MPI_LIBS}
# fftw
incLPATHfft = ${PATH_TO_FFTW_LIB}
incIPATHfft = ${PATH_TO_FFTW_INC}
incLIBSfft = ${FFTW_LIBS}
# HDF5
incLPATHhdf = ${PATH_TO_HDF5_LIB}
incIPATHhdf = ${PATH_TO_HDF5_INC}
incLIBShdf = ${HDF5_LIBS}
# c-compiler
CC = ${C_COMPILER} ${C_FLAGS}

# paths to things.
incbasedir = ../../inc/
bindir = ../../bin
libdir = ../../lib
EOF
# add the extra CGAL stuff
if [ ${USECGAL} == "yes" ]
then
cat>>Makefile.inc<<EOF
# CGAL stuff
USECGAL = yes
incLPATHcgal=${PATH_TO_CGAL_LIB}
incIPATHcgal=${PATH_TO_CGAL_INC}
CPP = ${CPP_COMPILER} ${CPP_FLAGS} \$(incLPATHcgal) \$(incIPATHcgal)
LDFLAGS =\$(incLPATHcgal) \$(incIPATHcgal) ${CGALLINKLINE}
EOF
fi

# a place to put all the .o and .mod files in, to not make the source messy
[ -d ../../build/libolle ] || mkdir ../../build/libolle
# and somewhere for the .mod files
[ -d ../../inc/libolle ] || mkdir ../../inc/libolle
# now, in theory, everything should be in place to compile
echo " "
echo "building libolle"
if [ ${CLEAN} = "YES" ]
then
    make clean && make -j ${NTHREADS_MAKE}
else
    make -j ${NTHREADS_MAKE}
fi
# good to know where this library ends up
libolledir=`pwd`
# now it should be ok! great success, hopefully!
cd ../../

# SECRET BLACK BELT TRICK: link some stuff to the inc thing
cd inc/libolle
    ln -sf ../../src/libolle/precompilerdefinitions .
cd ../../

# build the command line parser library, FLAP, since that is used by all the codes
cd src/libflap
    # make some space
    [ -d ../../build/libflap ] || mkdir ../../build/libflap
    [ -d ../../inc/libflap ] || mkdir ../../inc/libflap
    echo " "
    echo "building FLAP"
    # again, a neat place to put the object files
    [ -d ../../build/libflap ] || mkdir ../../build/libflap
    cp ../libolle/Makefile.inc .
    # and build it.
    if [ ${CLEAN} = "YES" ]
    then
        make clean && make -j ${NTHREADS_MAKE}
    else
        make -j ${NTHREADS_MAKE}
    fi
cd ../../

# now I should build a lot of codes. If we made it this far, it should be easy.
listofcodes="
dump_dynamical_matrices
phonon_dispersion_relations
crystal_structure_info
generate_structure
canonical_configuration
lineshape
samples_from_md
extract_forceconstants
atomic_distribution
pack_simulation
refine_structure
thermal_conductivity
modecoupling_transport
anharmonic_free_energy
phasespace_surface
"

#some things that are not quite ok yet
notdone="
generate_kpoints
"

# go through them and compile
for code in ${listofcodes}
do
    echo " "
    echo "Building ${code}"
    # go to the correct directory
    cd src/${code} || exit 0
        # make sure there is space in the build thingy
        [ -d ../../build/${code} ] || mkdir ../../build/${code}
	# make sure we have the system-specific things
	cp ../libolle/Makefile.inc .
        # build it
        if [ ${CLEAN} = "YES" ]
        then
            make clean && make -j ${NTHREADS_MAKE}
        else
            make -j ${NTHREADS_MAKE}
        fi
        if [ ${MANPAGE} = "YES" ]
        then
            # if I got a binary, run it with --manpage to create the manpage and move it to the right place
            [ -f ../../build/${code}/${code} ] && ../../build/${code}/${code} --manpage
            [ -f ${code} ] && ./${code} --manpage
            [ -f "${code}.1" ] && mv ${code}.1 ../../man/man1 # manpage
            [ -f "${code}.md" ] && mv ${code}.md ../../man/ # same thing in markdown format.
        fi
    cd ../../
    # link it to bin?
    [ -f build/${code}/${code} ] && ln -sf ../build/${code}/${code} bin/${code}
done

basedir=`pwd`

echo " "
echo "Printing bashrc_tdep, append these lines to your .bashrc for stuff to work nicely"

cat>bashrc_tdep<<EOF
MANPATH=\$MANPATH:${basedir}/man
PATH=\$PATH:${basedir}/bin
export MANPATH
export PATH
alias gnuplot='gnuplot -persist'
EOF

cat bashrc_tdep

echo " "
echo "Everything should be set up and good to go!"
