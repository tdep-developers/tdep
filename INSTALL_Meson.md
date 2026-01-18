Installation with the Meson build system
===

Alternativaly to the `build_things.sh` script, there is also the possibility to use [Meson](https://mesonbuild.com/). It is a build automation tool, and  it supports incremental builds. The dependencies should be installed in standard locations (e.g. `/usr/local/`) or specified in the `PKG_CONFIG_PATH`.

First setup the git version for the code:
```setup_git_version.sh```
Then you can run the configuration step:
```meson setup build```
And finally compile the code:
```
cd build
meson compile
meson install
```

If some dependencies are not found, please make sure that they are in your PKG_CONFIG_PATH. For example, put a line of the form in your `.bashrc` / `.bash_profile` :
```export PKG_CONFIG_PATH="/path/to/your/hdf5/:${PKG_CONFIG_PATH}"```
Depending on the method used to install the required libraries, they may not be automatically put inside the search path (Homebrew is known to not always do it).
You can make sure that `pkg-config` is able to find the dependencies by running: `pkg-config --libs hdf5`
Meson will first try to find dependencies via `pkg-config`. If it does not find them, it will try to use CMake (if installed/loaded).

Once the configuration step is done, everything should go smoothly. The binaries will be in build/bin/executable_name.
