project('tdep',
        'fortran', 'c', 'cpp',
        meson_version: '>=1.1',
        default_options: ['buildtype=custom', 'debug=true', 'optimization=3'])
add_project_arguments('-cpp', language: 'fortran')
add_project_arguments('-ffree-line-length-none', language: 'fortran')
# add_project_arguments('-std=f2008', language: 'fortran')
add_project_arguments('-fallow-argument-mismatch', language: 'fortran')

# global dependencies:

# mpi
dep_mpi = dependency('mpi', language: 'fortran')

# linear algebra
dep_linalg_found = false
# 1) try accelerate framework for macos (if wanted)
if build_machine.system() == 'darwin' and get_option('accfrmwrk') == true
  dep_linalg = dependency('appleframeworks', modules : 'accelerate', required: false)
  dep_linalg_found = dep_linalg.found()
endif

# 2) try blas/openblas + lapack otherwise
if not dep_linalg_found
  dep_blas = dependency('blas', required: false)
  if not dep_blas.found()
    message('Did not find the blas library. Trying with openblas')
    dep_blas = dependency('openblas', required: true) # required as it is the last try
  endif
  dep_lapack = dependency('lapack', required: true)
  dep_linalg = [dep_blas, dep_lapack]
endif
# # add some common path to pkg-config search path
# env = environ()
# env.prepend('PKG_CONFIG_PATH', '/usr/lib/x86_64-linux-gnu')

# message(get_option('linalg_flavor'))

# fft
dep_fftw = dependency('fftw3')

# hdf5
dep_hdf5 = dependency('hdf5', language: 'fortran', version: '>1.10.7')

# dep_all = [dep_mpi, dep_blas, dep_lapack, dep_fftw, dep_hdf5]
dep_all = [dep_mpi, dep_linalg, dep_fftw, dep_hdf5]

subdir('src/libolle')
subdir('src/libflap')

subdir('src/anharmonic_free_energy')
subdir('src/atomic_distribution')
subdir('src/canonical_configuration')
subdir('src/crystal_structure_info')
subdir('src/dump_dynamical_matrices')
subdir('src/extract_forceconstants')
subdir('src/generate_structure')
subdir('src/lineshape')
subdir('src/pack_simulation')
subdir('src/phasespace_surface')
subdir('src/phonon_dispersion_relations')
subdir('src/refine_structure')
subdir('src/samples_from_md')
subdir('src/thermal_conductivity')
subdir('src/thermal_conductivity_2023')
