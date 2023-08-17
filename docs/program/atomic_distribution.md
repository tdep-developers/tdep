
### Short description

Calculates properties of the atomic distribution from molecular dynamics, such as mean square displacement, pair distribution function, vector distribution functions and probability densities. Useful for analysing simulations close to instabilities/phase transitions to have some idea where the atoms are.

### Command line options:




Optional switches:

* `--cutoff value`, `-r value`  
    default value 8.0  
    Consider pairs up to this distance, in A.

* `--stride value`, `-s value`  
    default value 1  
    Use every N configuration instead of all.

* `--help`, `-h`  
    Print this help message

* `--version`, `-v`  
    Print version
### Examples

`atomic_distribution --cutoff 4.3` 

`atomic_distribution --cutoff 4.3 --notransform` 

### Long description

When using molecular dynamics you are not guaranteed to end up in the same crystal structure that you started with. For the TDEP method to work, you need to know what the crystal structure is, and be sure that it has not changed during the simulation. This code contains a few methods to determine if the structure stays stable during the simulation. So, before trying to get phonons and free energies and other things, it it useful to make sure that what is in your simulation box actually is what you think it is. A number of diagnostics is provided here.

#### Mean square displacement

The easiest measure is the mean square displacement:

$$
\textrm{msd}(t) = \frac{1}{N} \sum_i \left| \mathbf{r}_i(t)-\mathbf{r}_i(0) \right|^2
$$

where $\mathbf{r}_i(t)$ is the position of atom $i$ at time $t$. To be able to use the rest of the methods in this package the system needs to be solid. In a solid the mean square displacement fluctuates with time, but does not grow. If you find strange results, the mean square displacement is the first thing to check.

#### Radial distribution function

Having established that the system is at least still a solid we can look at the radial distribution function (or pair correlation function), defined as

$$
g(r) = \frac{ n(r) }{\rho 4 \pi r^2 dr}
$$

where $\rho$ is the mean particle density, and $n(r)$ the number of particles in an infinitesimal shell of width $dr$. Usually, this is averaged over all atoms in the system. In this code, I project it onto symmetrically equivalent pairs, yielding a projected pair distribution:

$$
g_i(r) = {\rho 4 \pi r^2 dr} \sum_i \delta\left( \left|r_i\right|-r \right)
$$

where the index $i$ corresponds to a coordination shell. The coordination shell is defined from the ideal lattice as set of pairs that can transform to each other via a spacegroup operation. Naturally, the sum over all projected PDFs yield the total.

<center>
<img src="/media/gan_pair_distribution.png" width="500" />
</center>


### Input files

* [infile.ucposcar](../files.md#infile.ucposcar)
* [infile.ssposcar](../files.md#infile.ssposcar)
* [infile.meta](../files.md#infile.meta)
* [infile.stat](../files.md#infile.stat)
* [infile.positions](../files.md#infile.positions)
* [infile.forces](../files.md#infile.forces)

### Output files

#### `outfile.mean_square_displacement.hdf5`

Below is a snippet that plots the mean square displacement.

```matlab
fn=outfile.mean_square_displacement.hdf5

figure(1); clf; hold on; box on;

n_atom=h5readatt(fn,'/','number_unique_atoms');

x=h5read(fn,'/time_axis');
y=h5read(fn,'/mean_square_displacement');
plot(x,y);
for i=1:n_atom
    y=h5read(fn,['/mean_square_displacement_atom_' num2str(i)]);
    plot(x,y)
end
```

Where the columns are time (in fs), mean square displacement (in Å<sup>2</sup>), followed by the partial mean square displacement per unique atom.

#### `outfile.pair_distribution_function.hdf5`

The hdf file is self-explainatory. Below is a matlab snippet that produces the plot above

```matlab
fn=outfile.pair_distribution.hdf5';

% fetch labels
lbl=strsplit(h5readatt(fn,'/','pair_shell_label'),'|')

% define some colors

clr=[
    0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840
];

figure(1); clf; hold on; box on;

for i=1:length(lbl)
    s=['/shell_' lbl{i} '/x']
    pdf.x{i}=h5read(fn,['/shell_' lbl{i} '/r']);
    pdf.y0{i}=h5read(fn,['/shell_' lbl{i} '/pdf']);

    plot(pdf.x{i},pdf.y0{i},'color',clr(i,:))
end

for i=1:length(lbl)
    s=['/shell_' lbl{i} '/x']
    xi=h5read(fn,['/shell_' lbl{i} '/ideal_peak_locations']);
    yi=interp1(pdf.x{i},pdf.y0{i},xi);
    plot(xi,yi,'color',clr(i,:),'marker','.','linestyle','none','markersize',15)
end

legend(lbl)
set(gca,'xminortick','on','yminortick','on')
```
.

