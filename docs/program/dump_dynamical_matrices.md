
### Short description

Small utility to extract the dynamical matrices from TDEP, in formats readable by other codes. For the moment the ABINIT DDB format is supported.

### Command line options:




Optional switches:

* `--qpoint_grid value#1 value#2 value#3`, `-qg value#1 value#2 value#3`  
    default value 26 26 26
    Density of q-point mesh for Brillouin zone integrations.

* `--meshtype value, for value in: 1,2,3,4`
    default value 1
    Type of q-point mesh. 1 Is a Monkhorst-Pack mesh, 2 an FFT mesh and 3 my fancy wedge-based mesh with approximately the same density the grid-based meshes. 4 build the commensurate mesh of an approximately cubic supercell.

* `--readqpoints`  
    default value .false.  
    Instead of generating a q-mesh, read it in fractional coordinates from a file called "infile.dynmatqpoints"

* `--help`, `-h`  
    Print this help message

* `--version`, `-v`  
    Print version

### Examples

`dump_dynamical_matrices` 

## Long summary

This utility takes the infile.forceconstant hdf5 file, reads it in and converts it to reciprocal space dynamical matrices, then dumps these on a set of q points, either a regular grid or a set of user selected points. The goal is to take renormalized high T TDEP phonon frequencies and eigenvectors and re-use them in other contexts.


### Input files

* [infile.forceconstant](../files.md#infile.forceconstant)

### Output files

#### `outfile.many_dynamical_matrices_DDB`


