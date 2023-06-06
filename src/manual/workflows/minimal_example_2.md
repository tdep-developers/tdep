title: Minimal example II
author: Olle Hellman

## Create the input files

In this example we want to produce input files similar to those provided in [the first minimal example](minimal_example_1.html) by setting up the required files, running VASP and postprocessing. But since this is a tutorial, you do not need to actually run any DFT calculations, the output files are provided in `examples/example_2_ScF3`.

### Preparing input

Start by making an empty working directory. The first step is to prepare the unit cell. In this particular case, use this and and name it `infile.ucposcar`

```
ScF3
   4.011000
  1.000000000000000   0.000000000000000   0.000000000000000
  0.000000000000000   1.000000000000000   0.000000000000000
  0.000000000000000   0.000000000000000   1.000000000000000
 Sc F
   1   3
Direct
  0.000000000000000   0.000000000000000   0.000000000000000
  0.500000000000000   0.000000000000000   0.000000000000000
  0.000000000000000   0.000000000000000   0.500000000000000
  0.000000000000000   0.500000000000000   0.000000000000000
```

The format is specified [here](../files.html#infile.ucposcar). To run molecular dynamics, you need a larger simulation cell. I have provided [tools](../../program/generate_structure.html) for that:

```
generate_structure -d 4 4 4
```
again, change the output to input

```
mv outfile.ssposcar infile.ssposcar
```
This is the preferred way of generating supercells. Any script should do, but some of them that exist out there do not use double precision for some odd reason. The [symmetry analysis](../../program/extract_forceconstants.html) is rather sensitive, and assumes that the unit- and supercell describe exactly the same lattice, to machine precision. Almost the same is not good enough.

This can not be stressed enough, in my experience this is the most common source of problems, and it is easily avoided. When in doubt, use at least 12 decimal places. A common problem is that .cif-files sometimes write too few digits for e.g. \(\sqrt(3)/2\) for my codes to correctly differentiate between crystal lattices. For example, this is really bad:

```
hcp Fe
 2.4
   0.50000  -0.866025  0.00000
   0.50000   0.866025  0.00000
   0.00000   0.000000  1.63299
Fe
   2
Direct
   0.333333333333333   0.666666666666667   0.250000000000000
   0.666666666666667   0.333333333333333   0.750000000000000
```

whereas this is good:

```
hcp Fe
 2.4
   0.500000000000000  -0.866025403784439  0.000000000000000
   0.500000000000000   0.866025403784439  0.000000000000000
   0.000000000000000   0.000000000000000  1.632993161855452
Fe
   2
Direct
   0.333333333333333   0.666666666666667   0.250000000000000
   0.666666666666667   0.333333333333333   0.750000000000000
```

It is a deliberate choice to set tolerances this tight: there is absolutely no reason why you should not be defining structures to high precision. Loose tolerances can cause subtle errors that are far more annoying to track down than entering those extra digits.

#### Use the primitive unit cell

The symmetry analysis is clever enough to ensure that you get the same forceconstants regardless of the cell you use, but you should use the smallest possible unit cell. There are [tools](../../program/remap_forceconstant.html) available to rearrange forceconstants to larger unit cells, but not from a larger cell to a smaller. So by using the smallest possible cell, you can always change your mind later.

### Running a VASP simulation

With the structure defined, we can run some molecular dynamics. In this pariticular example we do it cooking-show style: the output is already prepared. Nonetheless, this is what the input looks like (INCAR only, there is nothing special in the other input files):

```
ENCUT = 600
ISMEAR = 0
ISIF = 2
IBRION = 0
NSW = 10000
EDIFF = 1E-5
IALGO = 48

POTIM = 2
TEBEG = 150
SMASS = 0
```

There are some things that might seem nonstandard. The energy cutoff is very high. I start convergence testing at about double the default, `ENCUT = 2*ENMAX`, where the default is read from the POTCAR files. You can not use tetrahedron integration, you must use some smearing. Fermi smearing makes the most sense, but requires a really tight k-point grid. It differs from what you might normally do in that we require high-quality forces, and not much else. Speaking of that, the tag `ISIF = 2` is crucial, since without it forces will not get printed to `OUTCAR`.

At this stage, you would run the simulations. Starting molecular dynamics from ideal lattice positions usually require a long equilibration time that has to be discarded. This time can be minimized by using a better seed. A program to do this is provided, see [canonical configuration](../../program/canonical_configuration.html).

### Parsing VASP output

Pretend that you started VASP and let it run for a while, and copy `examples/example_2_ScF3/OUTCAR` to the working directory and run

```
process_outcar_5.3.py OUTCAR
```

and all the input files should be created. From this point, repeat the steps in [the first example](minimal_example_1.html). Note that these are not converged calculations by any means, the provided example is just 80 time steps using a single k-point for the electronic structure calculations (github is really not the place to put large output files, if I can figure out a place to host large files, I will put production quality files here).

### TLDR

The short version if you are impatient:

* get a unitcell
* make a supercell
* run VASP
* run `process_outcar_5.3.py`

I could have started with this, but I wanted you to understand a bit what you should and should not do.

### Things to understand:

Some things to consider:

*	Read the [documentation](../1_utilities.html) for `process_outcar_5.3.py`, figure out how to stitch a continuation job together.
* 	Maybe you ran VASP with some settings that makes the OUTCAR parsing fail. Make sure you understand how to construct the [input files](../files.html) manually.

