executable('thermal_conductivity_2023',
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

install_symlink('thermal_conductivity_2023', pointing_to: meson.project_build_root() / 'bin/thermal_conductivity_2023', install_dir: meson.project_source_root() / 'bin')
