
This page details the format for all input files.

### <a name="infile.ucposcar"></a> infile.ucposcar

An example of a crystal structure (this is exactly the VASP 5+ file format):

	Bi2Te3
	10.314046162
	 0.243360208934  0.000000000000  0.969935981757
	-0.121680104467  0.210756123207  0.969935981757
	-0.121680104467 -0.210756123207  0.969935981757
	Bi Te
	2 3
	direct coordinates
	0.599898812406 0.599898812406 0.599898812406
	0.400101187594 0.400101187594 0.400101187594
	0.791308614612 0.791308614612 0.791308614612
	0.208691385388 0.208691385388 0.208691385388
	0.000000000000 0.000000000000 0.000000000000

The first line is a comment, the second line a global scaling factor \(a\), and the three following lines the lattice vectors \(\mathbf{a}_1\),\(\mathbf{a}_2\),\(\mathbf{a}_3\). These will be multiplied by \(a\). Next is the specfication of what elements there are, and below that how many of each. The line direct coordinates specify that the coordinates of the atoms are in fractional coordinates. Each atom has it's own line, with two first specifying the Bi positions and the last three Te positions.

### <a name="infile.ssposcar"></a> infile.ssposcar

This file hold the ideal positions for the atoms in the supercell. The file format is identical to the unitcell, only that the cell is larger. There is no need for the unit and supercell to be perfectly commensurate: for example with the trigonal unitcell in the example above, I used a supercell based on the hexagonal unit cell. The only requirement is that the two cells describe the same lattice. And by the same, I do not mean almost the same, I mean the same to machine precision. You can not build the supercell from ideal positions and use relaxed positions for the unitcell, for example.

### <a name="infile.positions"></a> infile.positions

This file hold the positions. It looks something like this:

```
   8.2561407159200030E-005   2.2226748311341836E-003  0.99884779259537781
  0.99849619391764632        1.3299790395140015E-003  0.20126465487180664
  0.99942091881904371       0.99921744312226823       0.39925524034421311
  0.99726097226994259        9.3000027090734956E-004  0.60117595839769800
  0.99998106248734664       0.99877980293422230       0.80215876222138049
   2.7059648117732179E-004  0.20042853805668126        9.0607694878967109E-004
   4.7165392641841262E-004  0.20049530215652531       0.20005541126842488
  0.99978563110732910       0.20120110797551810       0.39998540127778304
   5.9697200870672466E-005  0.20002590513205151       0.60020279388531417
   1.0832057274455832E-003  0.19752092901626805       0.80125606638730629
   ...
```

and goes on for many many lines. The format is as follows:

| row ||||
|-----|-------------|---|---|
| 1   | \( r_1^x \) | \( r_1^y \) | \( r_1^z \) |
| \(\vdots\)   | \(\vdots\) | \(\vdots\) | \(\vdots\) |
| N   | \( r_N^x \) | \( r_N^y \) | \( r_N^z \) |
| N+1 | \( r_1^x \) | \( r_1^y \) | \( r_1^z \) |
| \(\vdots\)   | \(\vdots\) | \(\vdots\) | \(\vdots\) |
| 2N   | \( r_N^x \) | \( r_N^y \) | \( r_N^z \) |
| \(\vdots\)   | \(\vdots\) | \(\vdots\) | \(\vdots\) |

That is the positions in fractional coordinates for each atom in the simulation cell, in the same order as in [infile.ssposcar](#infile.ssposcar). Once all the positions are specified, it starts over with the next configuration, for a total of number of atoms times number of configurations lines.

### <a name="infile.forces"></a> infile.forces

The format for this file is nearly identical to [infile.positions](#infile.positions), only that

| row ||||
|-----|-------------|---|---|
| 1   | \( f_1^x \) | \( f_1^y \) | \( f_1^z \) |
| \(\vdots\)   | \(\vdots\) | \(\vdots\) | \(\vdots\) |
| N   | \( f_N^x \) | \( f_N^y \) | \( f_N^z \) |
| N+1 | \( f_1^x \) | \( f_1^y \) | \( f_1^z \) |
| \(\vdots\)   | \(\vdots\) | \(\vdots\) | \(\vdots\) |
| 2N   | \( f_N^x \) | \( f_N^y \) | \( f_N^z \) |
| \(\vdots\)   | \(\vdots\) | \(\vdots\) | \(\vdots\) |

each line is the force on each atom in Cartesian coordinates, in eV/Å. Again, a total of number of atoms times number of configurations lines.

### <a name="infile.stat"></a> infile.stat

This file contains all the energy and stresses in the following format:

|row|   |   |   |   |   |   |   |   |   |   |   |   |   |   |
|----|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | i | \(t\) | \(E_t\) | \(E_p\) | \(E_k\) | \(T\) | \(P\) | \(\sigma_{xx}\) | \(\sigma_{yy}\) |\(\sigma_{zz}\) | \(\sigma_{xz}\) | \(\sigma_{yz}\) | \(\sigma_{xy}\) |
| 2 | i | \(t\) | \(E_t\) | \(E_p\) | \(E_k\) | \(T\) | \(P\) | \(\sigma_{xx}\) | \(\sigma_{yy}\) |\(\sigma_{zz}\) | \(\sigma_{xz}\) | \(\sigma_{yz}\) | \(\sigma_{xy}\) |
| 3 | i | \(t\) | \(E_t\) | \(E_p\) | \(E_k\) | \(T\) | \(P\) | \(\sigma_{xx}\) | \(\sigma_{yy}\) |\(\sigma_{zz}\) | \(\sigma_{xz}\) | \(\sigma_{yz}\) | \(\sigma_{xy}\) |
| \(\vdots\) | \(\vdots\) | \(\vdots\) | \(\vdots\) | \(\vdots\) | \(\vdots\) | \(\vdots\) | \(\vdots\) | \(\vdots\) | \(\vdots\) | \(\vdots\) | \(\vdots\) | \(\vdots\) |\(\vdots\) |

The columns are, in order configuration (integer), time (in fs), total energy, potential energy, kinetic energy, temperature, pressure and stress tensor. Energies are in eV/supercell, temperature in K, pressures and stress tensor in GPa. There is one line per configuration in the simulation.

One line for every configuration in the simulation. The energies are in eV/supercell, temperature in K, pressure and stress in GPa.

!!! note
	The information in this file is only crucial when calculating the free energy. For other applications this can safely be filled with mock data, in case the method you use to calculate forces does not output e.g. stress tensors.

!!! note
	In case you use stochastic sampling, the time column is of no importance and can be set to any number. However, the column still needs to exist in the file.

### <a name="infile.meta"></a> infile.meta

Some information about the MD run, an example:

```
240    # N atoms
2000   # N timesteps
1.0    # timestep in fs
300    # temperature in K
```

The first line is the number of atoms in the simulation, the second line the number of timesteps. Then the timestep and the temperature.

!!! note
	The timestep and temperature can safely be set to 0 in case you use stochastic sampling.

### <a name="infile.lotosplitting"></a> infile.lotosplitting

This is an example input file for the long-range electrostatic corrections, also for Bi2Te3.

	 5.0  0.0  0.0 # three lines for dielectric tensor
	 0.0  5.0  0.0
	 0.0  0.0  7.5
	 8.0  0.0  0.0 # Born effective charge for atom 1, three lines
	 0.0  8.0  0.0
	 0.0  0.0  2.5
	 8.0  0.0  0.0 # Born effective charge for atom 2, three lines
	 0.0  8.0  0.0
	 0.0  0.0  2.5
	-4.5  0.0  0.0 # Born effective charge for atom 3, three lines
	 0.0 -4.5  0.0
	 0.0  0.0 -0.5
	-4.5  0.0  0.0 # Born effective charge for atom 4, three lines
	 0.0 -4.5  0.0
	 0.0  0.0 -0.5
	-7.0  0.0  0.0 # Born effective charge for atom 5, three lines
	 0.0 -7.0  0.0
	 0.0  0.0 -5.1

The first three lines is the dielectric tensor, after that three lines per atom in [infile.ucposcar](#infile.ucposcar) with the Born effective charges. The dielectric tensor is unitless, and the Born effective charges are in electron charges.

It is important to note that in general, the Born effective charges are not necessarily symmetric:

$$
Z_i^{\alpha\beta} = \frac{\partial^2 U}{\partial \epsilon_i^{\alpha} \partial u_i^{\beta}}
$$

That is derivative with respect to electric field and position. This means it matters how you enter them. In this input file, the convention is as follows:

$$
\begin{pmatrix}
	\partial_{\epsilon x}\partial_{ux} & \partial_{\epsilon x}\partial_{uy} & \partial_{\epsilon x}\partial_{uz} \\
	\partial_{\epsilon y}\partial_{ux} & \partial_{\epsilon y}\partial_{uy} & \partial_{\epsilon y}\partial_{uz} \\
	\partial_{\epsilon z}\partial_{ux} & \partial_{\epsilon z}\partial_{uy} & \partial_{\epsilon z}\partial_{uz}
\end{pmatrix}
$$
### <a name="infile.qpoints_dispersion"></a> infile.qpoints_dispersion

Many programs output properties along a path in the BZ. Per default, this is generated procedurally. In case the default is not satisfactory, you have the possibility to specify a path.

```
FCC                         ! Bravais lattice type
  100                       ! Number of points on each path
    4                       ! Number paths between special points
GM  X                       ! Starting and ending special point
X   U                       !
K   GM                      !
GM  L                       !
```

Where the first line specify the Bravais family, followed by the number of points on each line segment and the number of line segments. Each line segment is specified by two labels. Use [crystal structure info](program/crystal_structure_info.md) to get a list of the possible labels and their coordinates. If this is not flexible enough, an arbitrary path can be specified

```
CUSTOM                      !
  100                       ! Number of points on each path
    4                       ! Number paths between special points
0.000 0.000 0.000   0.000 0.500 0.500 GM X
0.000 0.500 0.500   0.000 0.625 0.375 X  U
0.375 0.750 0.375   0.000 0.000 0.000 K  GM
0.000 0.000 0.000   0.000 0.500 0.000 GM L
```

This is the same path as above, but explicitly specified. Each segment is specified by 3+3 numbers and two labels. The coordinates are reciprocal fractional.

### <a name="infile.isotopes"></a> infile.isotopes

The phonon many-body perturbation theory considers scattering by isotopes. By default, they will use the natural distribution (tabulated in the code, taken from the symbol in [infile.ucposcar](#infile.ucposcar)). In case you want to specify some other distribution, you can:

```
1		  # number of isotopes for first atom in infile.ucposcar
1 28.0855 # concentration, mass, one line per isotope
2		  # number of isotopes for second atom
0.5 12.0  # concentration, mass
0.5 13.0  # concentration, mass
...
```

Per atom in the unit cell, you specify the number of isotopes, followed by the appropriate number of concentrations and masses (in atomic mass units).
