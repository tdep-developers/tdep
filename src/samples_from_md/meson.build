executable('samples_from_md',
           'main.f90',
           'options.f90',
           link_with: [libolle, libflap],
           include_directories: ['../libolle', '../libflap'],
           dependencies: dep_all,
           install: true,
           install_dir: meson.project_build_root() / 'bin')

install_symlink('samples_from_md', pointing_to: meson.project_build_root() / 'bin/samples_from_md', install_dir: meson.project_source_root() / 'bin')
