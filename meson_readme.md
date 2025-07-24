Run first: meson_setup.py
Then:
meson setup build
cd build
meson compile

If some dependencies are not found, pleasee make sure that they are in yourr PKG_CONFIG_PATH. For example, put something of the form in your .bashrc/.bash_profile :
export PKG_CONFIG_PATH="/path/to/your/netcdf/:${PKG_CONFIG_PATH}"

Once the configuration step is done, everything should go smoothly. The binaries will be in build/src/code_name/executable_name.
