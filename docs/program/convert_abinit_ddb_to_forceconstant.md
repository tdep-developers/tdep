
### Short description

Generates "outfile.ucposcar" and "outfile.forceconstant" from abinit DDB files.

### Command line options:




Optional switches:

* `--qpoint_grid value#1 value#2 value#3`, `-qg value#1 value#2 value#3`  
    default value -1 -1 -1  
    Phonon q-point mesh used in Abinit.

* `--files value#1 [value#2...]`, `-f value#1 [value#2...]`  
    DDB filename.

* `--truncate`  
    default value .false.  
    Truncate the realspace IFCs.

* `--cutoff value`, `-rc value`  
    default value -1  
    Force a cutoff on the realspace ifcs.

* `--help`, `-h`  
    Print this help message

* `--version`, `-v`  
    Print version
author: Olle Hellman
display: none
graph: none
propname: convert abinit
propnamelink: <a href="../program/convert_abinit_ddb_to_forceconstant.html">convert abinit</a>
{!man/convert_abinit_ddb_to_forceconstant.md!}

### Longer summary

This is a utility to convert Abinit DDB files to the TDEP forceconstant format. A short disclaimer should be that it is not trivial to figure out what information is written to the abinit DDB files by parsing the DDB file alone. The conversion script has been tested on hundreds of Abinit calculations, but your mileage might vary.

The script will try to figure out the q-mesh used in the DFPT calculations automatically, but in case you used a shifted mesh, you might want to specify it manually, since the information about the q-mesh is not printed to the DDB file.

@todo Add io
