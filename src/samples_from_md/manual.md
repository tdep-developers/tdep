author: Olle Hellman
display: none
graph: none
propname: samples from md
propnamelink: <a href="../program/samples_from_md.html">samples from md</a>
{!man/samples_from_md.md!}

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

* [infile.ucposcar](../page/files.html#infile.ucposcar)
* [infile.ssposcar](../page/files.html#infile.ssposcar)
* [infile.meta](../page/files.html#infile.meta)
* [infile.stat](../page/files.html#infile.stat)
* [infile.positions](../page/files.html#infile.positions)
* [infile.forces](../page/files.html#infile.forces)

### Output files

This code will generate a series of structures given in VASP POSCAR [format](../page/files.html#infile.ucposcar) with positions in fractional coordinates and velocities in Å/fs.
