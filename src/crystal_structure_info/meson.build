executable('crystal_structure_info',
           'main.f90',
           'options.f90',
           link_with: [libolle, libflap],
           include_directories: ['../libolle', '../libflap'],
           dependencies: dep_all,
           install: true,
           install_dir: meson.project_build_root() / 'bin')

install_symlink('crystal_structure_info', pointing_to: meson.project_build_root() / 'bin/crystal_structure_info', install_dir: meson.project_source_root() / 'bin')
