executable('extract_forceconstants',
           'main.f90',
           'options.f90',
           link_with: [libolle, libflap],
           include_directories: ['../libolle', '../libflap'],
           dependencies: dep_all,
           install: true,
           install_dir: meson.project_build_root() / 'bin')

install_symlink('extract_forceconstants', pointing_to: meson.project_build_root() / 'bin/extract_forceconstants', install_dir: meson.project_source_root() / 'bin')
