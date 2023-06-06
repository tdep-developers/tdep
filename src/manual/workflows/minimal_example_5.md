title: Minimal example V
author: Olle Hellman

### Thermal conductivity

@todo Introduction, say what papers to read and so on.

@todo Specify what needs to exist

@todo Switch to Si, makes more sense.

Here we will reuse the data in `example_1_fcc_Al`. Copy the files to a new, clean folder, and repeat the steps from [the first tutorial](minimal_example_1.html) to generate second and third order force constants.

Note that this example is for fcc Al. I chose that since with one atom per unit cell, the thermal conductivity calculations are really fast. Unfortunately, since thermal conductivity in Al is dominated by the electronic component, none of the values can be compared with anything in a meaningful way -- but all the procedures in this tutorial are valid for any compound.

You should end up with (among other things)

* [infile.ucposcar](../page/files.html#infile.ucposcar)
* [infile.forceconstant](extract_forceconstants.html#infile.forceconstant)
* [infile.forceconstant_thirdorder](extract_forceconstants.html#infile.forceconstant_thirdorder)

This set of files is the minimum required to calculate thermal conductivity.

```
mpirun thermal_conductivity -qg 5 5 5 --temperature 300
```

and it should appear, rather quickly. Make yourself familiar with the options for [thermal conductivity](../../program/thermal_conductivity.html), there are quite a few. A series of temperatures can be calculated with

```
mpirun thermal_conductivity -qg 5 5 5 --temperature_range 10 1000 50 --logtempaxis
```

@todo Convergence with respect to q-point density. Point out that you have to plot with respect to 1/q, not q.

@todo Think about different integration schemes. Per default, the scattering rates are integrated using the tetrahedron method.[^Lehmann1972]

@todo Make them test tetrahedron vs gaussian.

[^Lehmann1972]: [Lehmann, G., & Taut, M. (1972). On the Numerical Calculation of the Density of States and Related Properties. Physica Status Solidi (B), 54(2), 469–477.]( http://doi.org/10.1002/pssb.2220540211)
