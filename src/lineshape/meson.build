exec_name='lineshape'
executable(exec_name,
           'dielscatter.f90',
           'dielscatter_helper.f90',
           # 'dielscatter_matrixelement.f90',
           'io.f90',
           'lineshape_helper.f90',
           'lo_realspace_selfenergy.f90',
           'lo_thermal_transport.f90',
           'main.f90',
           'options.f90',
           'phonondamping_aux.f90',
           'phonondamping_dos.f90',
           'phonondamping.f90',
           'phonondamping_grid.f90',
           'phonondamping_path.f90',
           'scatteringrates.f90',
           link_with: [libolle, libflap],
           include_directories: ['../libolle', '../libflap'],
           dependencies: dep_all,
           install: true,
           install_dir: meson.project_build_root() / 'bin')

install_symlink(exec_name, pointing_to: meson.project_build_root() / 'src' / exec_name / exec_name, install_dir: meson.project_source_root() / 'bin')
