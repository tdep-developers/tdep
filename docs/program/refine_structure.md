
### Short description

Small utility to ensure that the input structure satisfies all symmetries to high precision.

### Command line options:




Optional switches:

* `--tolerance_lattice value`, `-tl value`  
    default value 1E-5  
    Tolerance for the lattice.

* `--tolerance_internal value`, `-ti value`  
    default value 1E-5  
    Tolerance for internal degrees of freedom.

* `--unitcell value`, `-uc value`  
    default value infile.ucposcar  
    Filename for the unitcell to refine.

* `--supercell value`, `-ss value`  
    default value none  
    Filename for supercell to refine.

* `--prototype value`, `-pf value`  
    default value none  
    Prototype unitcell where you know the symmetry is correct.

* `--help`, `-h`  
    Print this help message

* `--version`, `-v`  
    Print version
### Examples

`refine_structure` 

### Longer summary

This is a utility to ensure that the symmetries of the crystal structure is satisfied to as high precision as possible. A vast majority of issues users experience, can be solved with this tool.

In essence, it takes an input crystal structure that might be specified like this

```
hcp Fe
 2.4
   0.50000  -0.866025  0.00000
   0.50000   0.866025  0.00000
   0.00000   0.000000  1.63299
Fe
   2
Direct
   0.333333   0.66666   0.25000
   0.666666   0.33333   0.75000
```

And returns

```
hcp Fe
       2.399999720250
    0.50000000000000    -0.86602540378444     0.00000000000000
    0.50000000000000     0.86602540378444     0.00000000000000
    0.00000000000000     0.00000000000000     1.63299057103649
 Fe
 2
Direct coordinates
  0.33333333333333   0.66666666666667   0.25000000000000  site 1 species 1: Fe
  0.66666666666667   0.33333333333333   0.75000000000000  site 2 species 1: Fe
```

!!! note
    Explain that we can use a prototype structure to determine the spacegroup, as in pick spacegroup from a unit cell I know is ok, and impose that on the current cell.


