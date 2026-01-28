exec_name='phonon_dispersion_relations'
executable(exec_name,
           'activity.f90',
           'densityplots.f90',
           'densityplots_stuntscattering.f90',
           'main.f90',
           'options.f90',
           'type_fast_interpolation.f90',
           'unfold_phonons.f90',
           'velocitydos.f90',
           link_with: [libolle, libflap],
           include_directories: ['../libolle', '../libflap'],
           dependencies: dep_all,
           install: true,
           install_dir: meson.project_build_root() / 'bin')

install_symlink(exec_name, pointing_to: meson.project_build_root() / 'src' / exec_name / exec_name, install_dir: meson.project_source_root() / 'bin')
