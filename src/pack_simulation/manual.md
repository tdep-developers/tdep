
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
