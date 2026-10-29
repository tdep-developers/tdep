exec_name='anharmonic_free_energy'
executable(exec_name,
           'epot.f90',
           'fourthorder.f90',
           'main.f90',
           'options.f90',
           'thermodynamic_helpers.f90',
           'thirdorder.f90',
           link_with: [libolle, libflap],
           include_directories: ['../libolle', '../libflap'],
           dependencies: dep_all)

install_symlink(exec_name, pointing_to: meson.project_build_root() / 'src' / exec_name / exec_name, install_dir: meson.project_source_root() / 'bin')
