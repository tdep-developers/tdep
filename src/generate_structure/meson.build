executable('generate_structure',
           'autocell.f90',
           'main.f90',
           'options.f90',
           link_with: [libolle, libflap],
           include_directories: ['../libolle', '../libflap'],
           dependencies: dep_all,
           install: true,
           install_dir: meson.project_build_root() / 'bin')

install_symlink('generate_structure', pointing_to: meson.project_build_root() / 'bin/generate_structure', install_dir: meson.project_source_root() / 'bin')
