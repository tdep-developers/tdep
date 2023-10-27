
### Short description

Small utility to pack a simulation to hdf5.

### Command line options:




Optional switches:

* `--stride value`, `-s value`  
    default value 1  
    Pack every N configuration instead of all.

* `--nrand value`  
    default value -1  
    Pack N random configurations instead of all.

* `--temperature value`  
    default value -1  
    Override the simulation temperature.

* `--variable_cell`, `-npt`  
    default value .false.  
    Make sure to store variable cell information.

* `--magnetic_moments`, `-mag`  
    default value .false.  
    Make sure to store projected magnetic moments.

* `--dielectric`  
    default value .false.  
    Make sure to store dielectric constant and Born charges.

* `--molecular_dynamics`, `-md`  
    default value .false.  
    Make sure to specify that this is real molecular dynamics.

* `--notidy`  
    default value .false.  
    Per default the simulation is cleaned, drift removed and so on. This switch skips that.

* `--help`, `-h`  
    Print this help message

* `--version`, `-v`  
    Print version
### Examples

`pack_simulation` 

## Long summary

This utility takes the plain-text input files and packs them to hdf5. The options allow you to pack a subset and specify what information to include. The resulting output file can take the place of the regular input files. This is by no means a necessary tool, the main utility is reduced file size and increased io speed.

### Input files

* [infile.ucposcar](../files.md#infile.ucposcar)
* [infile.ssposcar](../files.md#infile.ucposcar)
* [infile.positions](../files.md#infile.positions)
* [infile.forces](../files.md#infile.forces)
* [infile.stat](../files.md#infile.stat)
* [infile.meta](../files.md#infile.meta)

### Output files

#### `outfile.sim.hdf5`

This is merely an archive that contains all the information in the input files but packed to a single file. For the casual use it has little benefit except using considerably less space.
