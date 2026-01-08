executable('canonical_configuration',
           'main.f90',
           'options.f90',
           'semirandom.f90',
           link_with: [libolle, libflap],
           include_directories: ['../libolle', '../libflap'],
           dependencies: dep_all,
           install: true,
           install_dir: meson.project_build_root() / 'bin')

install_symlink('canonical_configuration', pointing_to: meson.project_build_root() / 'bin/canonical_configuration', install_dir: meson.project_source_root() / 'bin')
