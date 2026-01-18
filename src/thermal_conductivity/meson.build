exec_name='thermal_conductivity'
executable(exec_name,
           'cumulative_kappa.f90',
           'kappa.f90',
           'main.f90',
           'options.f90',
           'scattering.f90',
           link_with: [libolle, libflap],
           include_directories: ['../libolle', '../libflap'],
           dependencies: dep_all,
           install: true,
           install_dir: meson.project_build_root() / 'bin')

install_symlink(exec_name, pointing_to: meson.project_build_root() / 'src' / exec_name / exec_name, install_dir: meson.project_source_root() / 'bin')
