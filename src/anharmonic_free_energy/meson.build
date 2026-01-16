executable('anharmonic_free_energy',
           'energy.f90',
           'epot.f90',
           'main.f90',
           'options.f90',
           link_with: [libolle, libflap],
           include_directories: ['../libolle', '../libflap'],
           dependencies: dep_all,
           install: true,
           install_dir: meson.project_build_root() / 'bin')

install_symlink('anharmonic_free_energy', pointing_to: meson.project_build_root() / 'bin/anharmonic_free_energy', install_dir: meson.project_source_root() / 'bin')
