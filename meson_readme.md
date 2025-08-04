Run first: meson_setup.py
Then:
meson setup build
cd build
meson compile

If some dependencies are not found, please make sure that they are in your PKG_CONFIG_PATH. For example, put something of the form in your .bashrc/.bash_profile :
export PKG_CONFIG_PATH="/path/to/your/netcdf/:${PKG_CONFIG_PATH}"
Meson will first try to find dependencies via pkg-config. If it does not find them, it will try to use CMake (if installed).

Once the configuration step is done, everything should go smoothly. The binaries will be in build/bin/executable_name.
