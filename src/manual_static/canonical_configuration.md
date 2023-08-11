
### Short description

Use forceconstants or a Debye temperature to generate uncorrelated supercell configurations emulating a canonical ensemble. These configurations can be used to either start ensemble runs of molecular dynamics with negligible equilibration time, or be used to directly sample phase space.

### Command line options:




Optional switches:

* `--temperature value`, `-t value`  
    default value 300  
    Temperature to emulate

* `--nconf value`, `-n value`  
    default value 5  
    Number of configurations to generate

* `--quantum`  
    default value .false.  
    Use Bose-Einstein statistics instead of Maxwell-Boltzmann. That is, use $$ \sqrt{\frac{\hbar (2n+1) }{2 m \omega}} $$ as the mean normal mode amplitudes instead of the classical $$ \frac{1}{\omega}\sqrt{\frac{k_BT}{m}} $$

* `--output_format value`, `-of value`, value in: `1,2,3,4,5`  
    default value 1  
    Selects output format. 1 is VASP, 2 is Abinit, 3 is LAMMPS, 4 is FHI-Aims, 5 is Siesta. Default 1.

* `--mindist value`  
    default value -1  
    What is the smallest distance between two atoms allowed, in units of the nearest neighbour distance.

* `--debye_temperature value`, `-td value`  
    default value -1  
    Generate forceconstants that match a Debye temperature, and build displacements according to these. See details below.

* `--maximum_frequency value`, `-mf value`  
    default value -1  
    Generate forceconstants that match a maximum frequency (in THz), and build displacements according to these. See details below.

* `--help`, `-h`  
    Print this help message

* `--version`, `-v`  
    Print version
### Examples

`canonical_configuration -n 10 -t 300` 

`canonical_configuration -n 300 -t 0 --quantum` 

`canonical_configuration -n 20 -t 10 --quantum --debye_temperature 400` 
