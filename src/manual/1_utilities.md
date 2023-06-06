title: Utility scripts

Part of the package is a set of small tools that help with preparing input for post-processing. What they do is to parse output files and prepare intput files according to this [format](files.html). If you understand the format (and it is not particularly complicated), it is easy to generate a parsing script for any code. The minimal requirements are positions, forces and energies. Other quantities are there mostly for reference, and can safely be set to zero.

#### VASP

The scripts `process_outcar.py` and `process_outcar_5.3.py` are used to extract information from [VASP](https://www.vasp.at) OUTCAR files. The usage is straightforward:

```
process_outcar_5.3.py OUTCAR
```

This will create the relevant [input files](files.hmtl). Alternatively, is you used some statistical sampling method and have many VASP simulations that need to be stitched together, use

```
process_outcar_5.3.py *OUTCAR
```

That is, the script takes any number of OUTCARs as argument. The script assumes you have used `ISIF = 2`, since without it forces don't get written.

#### Abinit

@todo Get Antoine to send me his script.

#### LAMMPS

LAMMPS does not write anything by default. To get usable output, I make sure to add the following to the input file:

```
units metal

variable          st    equal step
variable          tm    equal step
variable          Et    equal etotal
variable          Ep    equal pe
variable          Ek    equal ke
variable          tmp   equal temp
variable          pr    equal press/10000
variable          sxx   equal pxx/10000
variable          syy   equal pyy/10000
variable          szz   equal pzz/10000
variable          sxy   equal pxy/10000
variable          sxz   equal pxz/10000
variable          syz   equal pyz/10000

fix statdump all print 100 "${st} ${tm} ${Et} ${Ep} ${Ek} ${tmp} ${pr} ${sxx} ${syy} ${szz} ${sxy} ${sxz} ${syz}" screen no file dump.stat
dump posdump all custom 100 dump.positions xs ys zs
dump forcedump all custom 100 dump.forces fx fy fz
dump_modify posdump format "%20.15e %20.15e %20.15e"
dump_modify forcedump format "%20.15e %20.15e %20.15e"
```

The frequency of the dumps and so on can of course be adjusted. Make sure you use [generate structure](../program/generate_structure.html) to create the LAMMPS cell file, since it will also create structure input files in the format TDEP needs. These dump files are not exactly the correct format. A small script can rearrange it to the correct format:

```bash
#!/bin/bash

# figure out how many atoms there are
na=`head -n 4 dump.forces | tail -n 1`
# remove the header from the stat file
grep -v '^#' dump.stat > infile.stat
# figure out how many timesteps there are
nt=`wc -l infile.stat | awk '{print $1}'`

# create the positions and force files
[ -f infile.forces ] && rm infile.forces
[ -f infile.positions ] && rm infile.positions
for t in `seq 1 ${nt}`
do
    nl=$(( ${na}+9))
    nll=$(( ${nl}*${t} ))
    echo "t ${t} ${nl} ${nll}"
    head -n ${nll} dump.forces | tail -n ${na} >> infile.forces
    head -n ${nll} dump.positions | tail -n ${na} >> infile.positions
done
```

There are probably far better ways of doing this, happy to take any suggestions.
