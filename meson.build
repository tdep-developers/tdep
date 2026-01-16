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
dep_blas = dependency('blas', required: false)
if not dep_blas.found()
  message('Did not find the blas library. Trying with openblas')
  dep_blas = dependency('openblas', required: true) # required as it is the last try
endif
dep_lapack = dependency('lapack')
# # add some common path to pkg-config search path
# env = environ()
# env.prepend('PKG_CONFIG_PATH', '/usr/lib/x86_64-linux-gnu')

# message(get_option('linalg_flavor'))

# fft
dep_fftw = dependency('fftw3')

# hdf5 + netcdf
dep_hdf5 = dependency('hdf5', language: 'fortran', version: '>1.10.7')
dep_netcdf = dependency('netcdf', language: 'fortran')

dep_all = [dep_mpi, dep_blas, dep_lapack, dep_fftw, dep_hdf5, dep_netcdf]

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
