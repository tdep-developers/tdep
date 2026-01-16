exec_name='thermal_conductivity_2023'
executable(exec_name,
           'main.f90',
           'mfp.f90',
           'options.f90',
           'pbe.f90',
           'phononevents.f90',
           'scatteringstrengths.f90',
           link_with: [libolle, libflap],
           include_directories: ['../libolle', '../libflap'],
           dependencies: dep_all,
           install: true,
           install_dir: meson.project_build_root() / 'bin')

install_symlink(exec_name, pointing_to: meson.project_build_root() / 'src' / exec_name / exec_name, install_dir: meson.project_source_root() / 'bin')
