author: Olle Hellman
display: none
graph: none
propname: pack simulation
propnamelink: <a href="../program/pack_simulation.html">pack simulation</a>
{!man/pack_simulation.md!}

## Long summary

This utility takes the plain-text input files and packs them to hdf5

### Input files

* [infile.ucposcar](../page/files.html#infile.ucposcar)
* [infile.ssposcar](../page/files.html#infile.ucposcar)
* [infile.positions](../page/files.html#infile.positions)
* [infile.forces](../page/files.html#infile.forces)
* [infile.stat](../page/files.html#infile.stat)
* [infile.meta](../page/files.html#infile.meta)

### Output files

#### `outfile.sim.hdf5`

This is merely an archive that contains all the information in the input files but packed to a single file. For the casual use it has little benefit except using considerably less space.
