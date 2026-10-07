
### Short description

Compute the symmetry-respecting average structure from an MD simulation with fixed cell (NVT).

### Command line options:




Optional switches:

* `--stride value`, `-s value`  
    default value 1  
    Use every N configuration instead of all.

* `--help`, `-h`  
    Print this help message

* `--version`, `-v`  
    Print version
### Examples

`average_structure` 

### Longer summary

This code computes the average structure from a molecular dynamics simulation such that the result has the symmetry of the reference structure. The internal degrees of freedom, i.e. the atomic displacements that are allowed by the space group of `infile.ucposcar`, are determined first. The net translation of the cell is removed by requiring that the displacements sum to zero.

For each timestep, the displacements of all atoms in the supercell from the reference positions are projected onto these internal degrees of freedom. The projections are averaged over the trajectory and the averaged structure is constructed from them. Displacements that break the symmetry are projected out, and the averaged structure keeps the space group of the reference. The lattice is not changed, so only simulations with fixed cell (NVT) are supported.

If the structure has no internal degrees of freedom, i.e. all positions are fixed by symmetry, the input structure is written unchanged.

### Input files

* [infile.ucposcar](../files.md#infile.ucposcar), determines the symmetry
* [infile.ssposcar](../files.md#infile.ssposcar)
* [infile.meta](../files.md#infile.meta)
* [infile.stat](../files.md#infile.stat)
* [infile.positions](../files.md#infile.positions)
* [infile.forces](../files.md#infile.forces)
* `infile.refposcar` (optional), the reference positions the displacements are measured from. Same format, number of atoms, species and atom order as `infile.ucposcar`. Defaults to the positions in `infile.ucposcar`.

### Output files

* `outfile.ucposcar`, the averaged unit cell
* `outfile.ssposcar`, the averaged supercell, consistent with `infile.ssposcar`
* `outfile.internal`, the time in fs (first column) and the value of each internal degree of freedom (remaining columns) for every timestep
