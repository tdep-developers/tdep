executable('dump_dynamical_matrices',
           'main.f90',
           'options.f90',
           link_with: [libolle, libflap],
           include_directories: ['../libolle', '../libflap'],
           dependencies: dep_all,
           install: true,
           install_dir: meson.project_build_root() / 'bin')

install_symlink('dump_dynamical_matrices', pointing_to: meson.project_build_root() / 'bin/dump_dynamical_matrices', install_dir: meson.project_source_root() / 'bin')
