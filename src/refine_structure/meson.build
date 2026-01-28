exec_name='refine_structure'
executable(exec_name,
           'lo_spacegroup_applyoperation.f90',
           'lo_spacegroup_charactertable.f90',
           'lo_spacegroup.f90',
           'lo_spacegroup_genoperations.f90',
           'lo_spacegroup_helpers.f90',
           'main.f90',
           'options.f90',
           'refine.f90',
           link_with: [libolle, libflap],
           include_directories: ['../libolle', '../libflap'],
           dependencies: dep_all,
           install: true,
           install_dir: meson.project_build_root() / 'bin')

install_symlink(exec_name, pointing_to: meson.project_build_root() / 'src' / exec_name/ exec_name, install_dir: meson.project_source_root() / 'bin')
