
### Short description

Choose representative uncorrelated samples from an MD simulation. The samples are chosen to be approximately evenly spaced, and reproduce the average potential energy, average kinetic energy have the same standard deviation of potential and kinetic energy.

### Command line options:




Optional switches:

* `--nsamples value`, `-n value`  
    default value 50  
    Number of samples

* `--output_format value`, `-of value`, value in: `1,2,3,4`  
    default value 1  
    Output format. 1 is VASP, 2 Abinit, 3 LAMMPS, 4 FHI-AIMS.

* `--help`, `-h`  
    Print this help message

* `--version`, `-v`  
    Print version
### Examples

`samples_from_md -n 100` 

### Longer summary

Ab initio molecular dynamics are expensive calculations. There will be a tradeoff between numerical precision and the number of timesteps. To work around this, you can run the MD with rather low precision and gather statistics. Then, from the long simulation, choose a set of uncorrelated configurations and recalculate these with high precision.

From these low accuracy calculations we choose a set of $n$ uncorrelated samples and correct scalar parameter $a$ as

$$
\begin{equation}
a = <a^l> + \frac{1}{n}\sum_{i=1}^n a^h_i-a_i^l,
\end{equation}
$$

where $a^l$ are the low accuracy calculations and $a^h$ are calculations done with high accuracy. This exploits the fact that most omissions of numerical accuracy, such as basis set and k-point selection, lead to additive errors. This technique is well suited to determine the interatomic force constants and resulting thermodynamic/transport properties with high accuracy.

This code selects a choice of uncorrelated samples from BOMD via a Monte-Carlo algorithm, assuring the selection is not biased. We start with a calculation of average potential $E_p$, kinetic energies $E_k$, and their standard deviation. We check the distance between samples assuring that chosen samples are not temporally adjacent. The results of this procedure is written in output files (`outfile.stat_sample`). The average values and distance between selected points depend on the number of desired samples.

### Input files

* [infile.ucposcar](../files.md#infile.ucposcar)
* [infile.ssposcar](../files.md#infile.ssposcar)
* [infile.meta](../files.md#infile.meta)
* [infile.stat](../files.md#infile.stat)
* [infile.positions](../files.md#infile.positions)
* [infile.forces](../files.md#infile.forces)

### Output files

This code will generate a series of structures in the specified output format.


