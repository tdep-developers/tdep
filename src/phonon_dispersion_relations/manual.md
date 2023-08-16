
### What does this code produce?

This code calculates phonon dispersion relations, usually presented like this:

<center><img src="/media/pbte_phonon_dispersions.png" width="600" /></center>

Where the color of the line indicates the projection onto species. The same data can also be presented as a density of states, like this:

<center><img src="/media/pbte_phonon_dos.png" width="500" /></center>

By default, minimal input is needed, but many options to tailor the output exist.

### Equations of motion in a harmonic crystal

As detailed in [extract forceconstants](extract_forceconstants.md) we have mapped the Born-Oppenheimer to an effective potential. A model Hamiltonian truncated at the second order (a harmonic Hamiltonian) is given by[^Born1998]

$$
\begin{equation}
\hat H=\sum_{\kappa} \frac{\mathbf{p}_{\kappa}^2}{2m_{\kappa}} + \frac{1}{2} \sum_{\kappa\lambda}\sum_{\alpha\beta}  \Phi_{\kappa\lambda}^{\alpha\beta} u^{\alpha}_{\kappa} u^{\beta}_{\lambda}\,
\end{equation}
$$

Where $\kappa,\lambda$ are indices to atoms, $\alpha,\beta$ Cartesian indices and $\mathbf{u}$ a displacement from equilibrium positions. To clarify, the instantaneous position of an atom is given by

$$
\begin{equation}
\mathbf{r}_{\kappa} = \mathbf{R}_{\mu} + \boldsymbol{\tau}_{i} + \mathbf{u}
\end{equation}
$$

such that $\mathbf{R} +\boldsymbol{\tau}$ denote the equilibrium position of the atom. The lattice vector $\mathbf{R}$ locate the unit cell, $\boldsymbol{\tau}$ locate the atom in the unit cell, and $\mathbf{u}$ how far the atom has moved from the reference positions. Ignoring surface effects, we can write the equations of motion:

$$
\begin{equation}
\ddot{\mathbf{u}}_{\mu i} m_{\mu i} =
-\sum_{\nu j} \mathbf{\Phi}_{\mu i,\nu j}\mathbf{u}_{\nu j} \,.
\end{equation}
$$

Here we expanded the atom index $\kappa$ to the double index indicating cell $\mu$ and position $i$. Since we deal with infinite crystals with periodic boundary conditions, the equations of motion will not depend on cell index $\mu$, and that index can be dropped.

$$
\begin{equation}\label{eq:eqmotion1}
\ddot{\mathbf{u}}_{i} m_{i} =
-\sum_{\nu j} \mathbf{\Phi}_{\mu i,\nu j}\mathbf{u}_{\nu j}.
\end{equation}
$$

I prefer a truncated notation:

$$
\begin{equation}
\ddot{\mathbf{u}}_{i} m_{i} =
-\sum_{\mathbf{R}} \mathbf{\Phi}_{ij}\left( \mathbf{R} \right)\mathbf{u}_{j}.
\end{equation}
$$

Where $\mathbf{R}$ implies the cell index corresponding to atom $j$. To solve these equations of motion, we can use a plane-wave ansatz:

$$
\begin{equation}
\mathbf{u}_{i}=\frac{1}{\sqrt{m_i}} \sum_{\mathbf{q}} A_{\mathbf{q}} \boldsymbol{\epsilon}_{\mathbf{q}}^{i} e^{ i \left(\mathbf{q} \cdot \mathbf{R} - \omega t \right) }.
\end{equation}
$$

Here the displacements are expressed as a sum of plane waves, or normal modes, each with wave vector $\mathbf{q}$ and frequency $\omega$. $\boldsymbol{\epsilon}$ is a polarisation vector, and $A_{\mathbf{q}}$ the normal mode amplitude. Substituting this into $\eqref{eq:eqmotion1}$ and exploiting the orthogonality of planes waves gives,

$$
\begin{split}
\sqrt{m_i}\sum_{\mathbf{q}} \omega^2 A_{\mathbf{q}} \boldsymbol{\epsilon}_{\mathbf{q}}^{i} e^{ i \left(\mathbf{q} \cdot \mathbf{R}_i - \omega t \right) }
= &
\sum_{\mathbf{R}}\sum_{\mathbf{q}}
\frac{A_{\mathbf{q}} \boldsymbol{\epsilon}_{\mathbf{q}}^{i}}{\sqrt{m_j}}
\mathbf{\Phi}_{ij}\left( \mathbf{R}_j \right)
e^{ i \left(\mathbf{q} \cdot \mathbf{R}_j - \omega t \right) } \\
%
%
\sum_{\mathbf{q}} \omega^2 A_{\mathbf{q}} \boldsymbol{\epsilon}_{\mathbf{q}}^{i} e^{ i \mathbf{q} \cdot \mathbf{R}_i }
= &
\sum_{\mathbf{R}}\sum_{\mathbf{q}}
\frac{A_{\mathbf{q}} \boldsymbol{\epsilon}_{\mathbf{q}}^{i}}{\sqrt{m_i m_j}}
\mathbf{\Phi}_{ij}\left( \mathbf{R}_j \right)
e^{ i \mathbf{q} \cdot \mathbf{R}_j } \\
%
%
\sum_{\mathbf{q}'} \omega^2 A_{\mathbf{q}'} \boldsymbol{\epsilon}_{\mathbf{q}'}^{i} e^{ i \mathbf{q}' \cdot \mathbf{R}_i }
e^{ -i \mathbf{q} \cdot \mathbf{R}_i }
= &
\sum_{\mathbf{R}}\sum_{\mathbf{q}'}
\frac{A_{\mathbf{q}'} \boldsymbol{\epsilon}_{\mathbf{q}'}^{i}}{\sqrt{m_i m_j}}
\mathbf{\Phi}_{ij}\left( \mathbf{R}_j \right)
e^{ i \mathbf{q}' \cdot \mathbf{R}_j }
e^{ -i \mathbf{q} \cdot \mathbf{R}_i } \\
%
%
\omega^2 A_{\mathbf{q}} \boldsymbol{\epsilon}_{\mathbf{q}}^{i}
= &
A_{\mathbf{q}} \boldsymbol{\epsilon}_{\mathbf{q}}^{i}
\sum_{\mathbf{R}}
\frac{\mathbf{\Phi}_{ij}\left( \mathbf{R}_j \right)}{\sqrt{m_i m_j}}
e^{ i \mathbf{q} \cdot \mathbf{R}_j } \\
%
%
\omega^2_{\mathbf{q}} \boldsymbol{\epsilon}_{\mathbf{q}}= & \mathbf{\Phi}(\mathbf{q}) \boldsymbol{\epsilon}_{\mathbf{q}}\,,
\end{split}
$$

where the periodic boundary conditions limit the choices of $\mathbf{q}$ to a wave vector $\mathbf{q}$ in the first Brillouin zone. The dynamical matrix $\mathbf{\Phi}(\mathbf{q})$ is given by

$$
\begin{equation}
\mathbf{\Phi}(\mathbf{q})=
\begin{pmatrix}
\mathbf{\Phi}_{11}(\mathbf{q}) & \cdots & \mathbf{\Phi}_{N 1}(\mathbf{q}) \\
\vdots & \ddots & \vdots \\
\mathbf{\Phi}_{1N}(\mathbf{q}) & \cdots & \mathbf{\Phi}_{N N}(\mathbf{q}) \\
\end{pmatrix}
\end{equation}
$$

This is the Fourier transform of the mass-weighted force constant matrix, where each $3 \times 3$ submatrix of the full $3N \times 3N$ is given by

$$
\begin{equation}
\mathbf{\Phi}_{ij}(\mathbf{q})=
\sum_{\mathbf{R}}
\frac{ \mathbf{\Phi}_{ij}(\mathbf{R}) }{\sqrt{m_i m_j}}
e^{i\mathbf{q}\cdot \mathbf{R}} \,.
\end{equation}
$$

The eigenvalues $\omega^2_{\mathbf{q}s}$ and eigenvectors $\boldsymbol{\epsilon}_{\mathbf{q}s}$ of the dynamical matrix denote the possible normal mode frequencies and polarizations of the system. The eigenvalues have the same periodicity as the reciprocal lattice; hence it is convenient to limit the solution to all vectors $\mathbf{q}$ in the first Brillouin zone. It is also worth defining the partial derivatives of the dynamical matrix:

$$
\begin{align}
\frac{\partial \mathbf{\Phi}_{ij}(\mathbf{q})}{\partial q_\alpha}
 & =
\sum_{\mathbf{R}}
iR_\alpha
\frac{ \mathbf{\Phi}_{ij}(\mathbf{R}) }{\sqrt{m_i m_j}}
e^{i\mathbf{q}\cdot \mathbf{R}}. \\
\frac{\partial^2 \mathbf{\Phi}_{ij}(\mathbf{q})}{\partial q_\alpha \partial q_{\beta}}
 & =
-\sum_{\mathbf{R}}
R_{\alpha} R_{\beta}
\frac{ \mathbf{\Phi}_{ij}(\mathbf{R}) }{\sqrt{m_i m_j}}
e^{i\mathbf{q}\cdot \mathbf{R}}.
\end{align}
$$

By calculating the eigenvalues and eigenvectors of the dynamical matrix over the Brillouin zone, all thermodynamic quantities involving the atomic motions can be determined, as far as the harmonic approximation is valid. If you are curious where the amplitudes of the normal modes disappeared, see [canonical configuration](canonical_configuration.md).

### Long-ranged interactions in polar materials

When determining the dynamical matrix, the sum over lattice vectors $\mathbf{R}$ should in principle be carried out over all pairs. In practice, we assume that the interactions are zero beyond some cutoff (see [extract forceconstants](extract_forceconstants.md) option `-rc2`), and truncate the sum. This is not a valid approach when dealing with polar material where the induced dipole-dipole interactions are essentially infinite-ranged. See for example Gonze & Lee.[^Gonze1994] [^Gonze1997] When dealing with these long-ranged electrostatics we start by defining the Born effective charge tensor:

$$
Z_{i}^{\alpha\beta} = \frac{\partial^2 U}{\partial \varepsilon^{\alpha} \partial u_{i}^{\beta}}
$$

That is the mixed derivative of the energy with respect to electric field $\varepsilon$ and displacement $\mathbf{u}$. Displacing an atom from its equilibrium position induces a dipole, given by

$$
\mathbf{d}_{i} = \mathbf{Z}_{i}\mathbf{u}_{i}
$$

The pairwise dipole-dipole interactions can be expressed as forceconstants:

$$
\begin{align}
\Phi^{\textrm{dd}}|^{\alpha\beta}_{ij} & = \sum_{\gamma\delta} Z_{i}^{\alpha\gamma}Z_{j}^{\beta\delta} \widetilde{\Phi}^{\gamma\delta}_{ij} \\
\widetilde{\Phi}^{\alpha\beta}_{ij} & =
\frac{1}{4\pi\epsilon_0}
\frac{1}{\sqrt{\det \epsilon}}
\left(
	\frac
	{ \widetilde{\epsilon}^{\alpha\beta} }
	{|\Delta_{ij}|_{\epsilon}^{3}}
	-3 \frac
	{\Delta_{ij}^{\alpha}\Delta_{ij}^{\beta}}
	{|\Delta_{ij}|_{\epsilon}^{5}}
\right)
\end{align}
$$

Here $\boldsymbol{\epsilon}$ is the dielectric tensor, $\widetilde{\boldsymbol{\epsilon}}=\boldsymbol{\epsilon}^{-1}$ its inverse and $\Delta$ realspace distances using the dielectric tensor as a metric:

$$
\begin{align}
\mathbf{r}_{ij} & = \mathbf{R}_{j}+\mathbf{\tau}_j-\mathbf{\tau}_i \\
\mathbf{\Delta}_{ij} & =\widetilde{\boldsymbol{\epsilon}}\mathbf{r}_{ij} \\
|\Delta_{ij}|_{\epsilon} & = \sqrt{\mathbf{\Delta}_{ij} \cdot \mathbf{r}_{ij}}
\end{align}
$$

This poses an issue when calculating the dynamical matrix: the interactions die off as $~1/r^{3}$ which makes it necessary to extend the sum over lattice vectors to infinity. This is remedied with the usual Ewald technique:

$$
\widetilde{\mathbf{\Phi}}_{ij} (\mathbf{q})=
\widetilde{\mathbf{\Phi}}^\textrm{r}+\widetilde{\mathbf{\Phi}}^\textrm{q}+\widetilde{\mathbf{\Phi}}^\textrm{c}
$$

Dividing the sum into a realspace part, a reciprocal part and a connecting part. The realspace part is given by

$$
\begin{align}
\widetilde{\mathbf{\Phi}}^\textrm{r}_{ij} & =
-\frac{\Lambda^3}{4\pi\epsilon_0\sqrt{\det \epsilon}}
\sum_{\mathbf{R}}
\mathbf{H}(\Lambda \Delta_{ij},\Lambda|\Delta_{ij}|_{\epsilon})
e^{i \mathbf{q} \cdot \mathbf{R}} \\
%%
\frac{\partial \widetilde{\mathbf{\Phi}}^\textrm{r}_{ij}}{\partial q_\alpha}  & =
-\frac{\Lambda^3}{4\pi\epsilon_0\sqrt{\det \epsilon}}
\sum_{\mathbf{R}}
iR_\alpha
\mathbf{H}(\Lambda \Delta_{ij},\Lambda|\Delta_{ij}|_{\epsilon})
e^{i \mathbf{q} \cdot \mathbf{R}}
\end{align}
$$

where

$$
H_{\alpha\beta}(\mathbf{x},y) = \frac{x_{\alpha}x_{\beta} }{y^2}
\left[
\frac{3\,\textrm{erfc}\,y}{y^3}
+
\frac{2 e^{-y^2}}{\sqrt{\pi}}\left(\frac{3}{y^2}+2 \right)
\right]
-\widetilde{\epsilon}_{\alpha\beta}
\left[
\frac{\textrm{erfc}\,y}{y^3} + \frac{2 e^{-y}}{\sqrt{\pi} y^2 }
\right]\,.
$$

In reciprocal space we have

$$
\begin{align}
%
% Not derivative
%
	\widetilde{\Phi}_\textrm{q}|_{ij}
	&=
	\sum_{\mathbf{K}=\mathbf{G}+\mathbf{q}}
	\chi_{ij}(\mathbf{K},\Lambda)
	\left( \mathbf{K}\otimes\mathbf{K} \right) \\
%
% X-direction
%
	\frac{\partial}{\partial q_x} \widetilde{\Phi}_\textrm{q}|_{ij}
	%
	&=
	\sum_{\mathbf{K}=\mathbf{G}+\mathbf{q}}
	\chi_{ij}(\mathbf{K},\Lambda)
	\left( \mathbf{K}\otimes\mathbf{K} \right)
	%
	\left(
	i\tau^x_{ij} - \left[\sum_{\alpha} K_\alpha \epsilon_{\alpha x}
	\right]
	\left[ \frac{1}{\|\mathbf{K}\|_{\epsilon}}+\frac{1}{4\Lambda^2} \right]
	\right) +
	%
	\chi_{ij}(\mathbf{K},\Lambda)
	\begin{pmatrix}
	2K_x & K_y & K_z \\
	K_y & 0 & 0 \\
	K_z & 0 & 0
	\end{pmatrix} \\
%
% Y-direction
%
	\frac{\partial}{\partial q_y} \widetilde{\Phi}_\textrm{q}|_{ij}
	%
	&=
	\sum_{\mathbf{K}=\mathbf{G}+\mathbf{q}}
	\chi_{ij}(\mathbf{K},\Lambda)
	\left( \mathbf{K}\otimes\mathbf{K} \right)
	%
	\left(
	i\tau^y_{ij} - \left[\sum_{\alpha} K_\alpha \epsilon_{\alpha y}
	\right]
	\left[ \frac{1}{\|\mathbf{K}\|_{\epsilon}}+\frac{1}{4\Lambda^2} \right]
	\right) +
	%
	\chi_{ij}(\mathbf{K},\Lambda)
	\begin{pmatrix}
	0 & K_x & 0 \\
	K_y & 2K_y & K_z \\
	0 & K_z & 0
	\end{pmatrix}\\
%
% Z-direction
%
	\frac{\partial}{\partial q_z} \widetilde{\Phi}_\textrm{q}|_{ij}
	%
	&=
	\sum_{\mathbf{K}=\mathbf{G}+\mathbf{q}}
	\chi_{ij}(\mathbf{K},\Lambda)
	\left( \mathbf{K}\otimes\mathbf{K} \right)
	%
	\left(
	i\tau^z_{ij} - \left[\sum_{\alpha} K_\alpha \epsilon_{\alpha z}
	\right]
	\left[ \frac{1}{\|\mathbf{K}\|_{\epsilon}}+\frac{1}{4\Lambda^2} \right]
	\right) +
	%
	\chi_{ij}(\mathbf{K},\Lambda)
	\begin{pmatrix}
	0 & 0 & K_x \\
	0 & 0 & K_y \\
	K_z & K_y & 2K_z
	\end{pmatrix}
\end{align}
$$

where

$$
\begin{align}
\chi_{ij}(\mathbf{K},\Lambda) &  =
\frac{1}{\Omega\epsilon_0}
\frac{
\exp\left(i \mathbf{q} \cdot \mathbf{\tau}_{ij} \right)
\exp\left( -\frac{\|\mathbf{K}\|_{\epsilon}}{4\Lambda^2}  \right)
}
{\|\mathbf{K}\|_{\epsilon}} \\
%%
\|\mathbf{K}\|_{\epsilon} & =\sum_{\alpha\beta}\epsilon_{\alpha\beta}K_\alpha K_\beta
\end{align}
$$

and finally a connecting part given by

$$
\widetilde{\mathbf{\Phi}}_\textrm{c}|_{ij}=
\delta_{ij}
\frac{\Lambda^3}{3 \epsilon_0 \pi^{3/2} \sqrt{ \det \boldsymbol{\epsilon} } }
\widetilde{\boldsymbol{\epsilon}}
$$

This expression is invariant for any $\Lambda>0$, which is chosen to make the real and reciprocal sums converge with as few terms as possible. Note that the Born charges are factored out from $\widetilde{\mathbf{\Phi}}$, and have to by multiplied in as described above. The derivatives with respect to $q$ are reproduced because they are useful when calculating the gradient of the dynamical matrix, needed when determining group velocities. In the long-wavelength limit, $\mathbf{q}\rightarrow 0$, the dipole-dipole dynamical matrix reduces to the familiar non-analytical term:

$$
\lim_{\mathbf{q}\rightarrow 0} \mathbf{\Phi}(\mathbf{q}) =
\frac{1}{\Omega\epsilon_0}
\frac{
\left(\mathbf{Z}_i\mathbf{q}\right) \otimes \left(\mathbf{Z}_j\mathbf{q}\right)
}{\mathbf{q}^T\boldsymbol{\epsilon}\mathbf{q}}
$$

### Separation into short- and long-range interactions



#### Other dispersive properties

We determine group velocities via the Hellman-Feynmann theorem:

$$
\begin{equation}
	\nabla_{\mathbf{q}} \omega^2_s = \left\langle \epsilon_{\mathbf{q}s} \right| \nabla_{\mathbf{q}} \Phi(\mathbf{q})  \left| \epsilon_{\mathbf{q}s} \right\rangle
\end{equation}
$$

where the derivatives of the dynamical matrix are defined above, and

$$
\begin{equation}
\mathbf{v}_{\mathbf{q}s} = \nabla_{\mathbf{q}} \omega_s(\mathbf{q})=
\frac{1}{2\omega_s} \nabla_{\mathbf{q}} \omega_s^2(\mathbf{q})
\end{equation}
$$

For degenerate modes this is ill-defined, we have to apply degenerate perturbation theory to resolve it. The practical procedure is as follows (this example is for a single component, but the same procedure works in general), first define

$$
\begin{equation}
    h_{ij} =
    \left\langle \epsilon_{\mathbf{q}i} \right|
    \partial_{q_{\alpha}} \Phi(\mathbf{q})  \left| \epsilon_{\mathbf{q}j} \right\rangle
\end{equation}
$$

where the indices run over the degenerate subspace. The gradient is then

$$
\begin{equation}
    \partial_{q_{\alpha}} \omega^{2}_{\mathbf{q}i} = \frac{1}{N} \sum_{i} \lambda_{i}
\end{equation}
$$

where $\lambda_i$ are the eigenvalues of $h_{ij}$.

The Hessian with respect to $q$ is given via

$$
\begin{align}
    \frac{\partial^2 \omega^{2}_{\lambda}}{\partial q_{\alpha} \partial q_{\beta}} & =
    2\left\langle \epsilon_{\lambda} \right|
        \partial_{q_{\alpha}}\partial_{q_{\beta}}\Phi(\mathbf{q})
    \left| \epsilon_{\lambda} \right\rangle +
    \left\langle \epsilon_{\lambda} \right| \partial_{q_{\alpha}}\Phi(\mathbf{q}) \left| \partial_{q_{\beta}} \epsilon_{\lambda} \right\rangle
    +
    \left\langle \epsilon_{\lambda} \right| \partial_{q_{\beta}}\Phi(\mathbf{q}) \left| \partial_{q_{\alpha}} \epsilon_{\lambda} \right\rangle
    %
    \\
    %
    \frac{\partial^2 \omega_{\lambda}}{\partial q_{\alpha} \partial q_{\beta}} & = \frac{1}{2\omega_{\lambda}}
    \left[
        \frac{\partial^2 \omega^{2}_{\lambda}}{\partial q_{\alpha} \partial q_{\beta}} -
        2 \frac{\partial \omega_{\lambda}}{\partial q_{\alpha}} \frac{\partial \omega_{\lambda}}{\partial q_{\beta}}
    \right]
\end{align}
$$

The mode Grüneisen parameters are a measure of the sensitivity of the vibrational frequencies to volume changes. They are given by

$$
\begin{equation}
\gamma_{\mathbf{q}s}=-\frac{V}{\omega_{\mathbf{q}s}}\frac{\partial \omega_{\mathbf{q}s}}{\partial V}
\end{equation}
$$

where $V$ is the volume and $\omega_{\mathbf{q}s}$ is the frequency of mode $s$ at wave vector $\mathbf{q}$. $\gamma_{\mathbf{q}s}$ can be obtained either by numerical differentiation of the phonon dispersion relations or from the third order force constants based on a perturbation approach:

$$
\begin{equation}
\gamma_{\mathbf{q}s}=-\frac{1}{6\omega_{\mathbf{q}s}^2}\sum_{ijk\alpha\beta\gamma}
\frac{\epsilon^{i\alpha\dagger}_{\mathbf{q}s} \epsilon^{j\beta}_{\mathbf{q}s}}
{\sqrt{m_i m_j}}
r_k^\gamma \Phi_{ijk}^{\alpha\beta\gamma}e^{i\mathbf{q}\cdot\ \mathbf{r}_j}
\end{equation}
$$

Here $\epsilon_{i\alpha}^{\mathbf{q}s}$ is component $\alpha$ associated eigenvector $\epsilon$ for atom $i$. $m_i$ is the mass of atom $i$, and $\mathbf{r}_i$ is the vector locating its position.

#### Density of states

Using `--dos` or will calculate the phonon density of states, given by

$$
g_s(\omega) = \frac{(2\pi)^3}{V} \int_{\mathrm{BZ}} \delta( \omega - \omega_{\mathbf{q}s}) d\mathbf{q}
$$

per mode $s$, summing those up yields the total. The site projected density of states for site $i$ is given by

$$
g_i(\omega) = \frac{(2\pi)^3}{V} \sum_{s} \int_{\mathrm{BZ}} \left|\epsilon^i_{\mathbf{q}s}\right|^2 \delta( \omega -  \omega_{\mathbf{q}s}) d\mathbf{q}
$$

that, also sums to the total density of states.

### Input files

* [infile.ucposcar](../files.md#infile.ucposcar)
* [infile.forceconstant](extract_forceconstants.md)

Optional files:

* [infile.qpoints_dispersion](../files.md#infile.qpoints_dispersion) (to specify a path for the phonon dispersions)
* [infile.forceconstant_thirdorder](extract_forceconstants.md) (for the Grüneisen parameter)

### Output files

Depending on options, the set of output files may differ. We start with the basic files that are written.

#### `outfile.dispersion_relations.hdf5`

This is a packed up format of all the data. Below is a short matlab snippet to plot parts of it, the hdf5 file is self-explanatory -- it contains all harmonic properties you can think of including some more.

```matlab
fn='outfile.dispersion_relations.hdf5';
x=h5read(fn,'/q_values');
xtck=h5read(fn,'/q_ticks');
xtckl=strsplit(h5readatt(fn,'/','q_tick_labels'));
y=h5read(fn,'/frequencies');
frequnit=h5readatt(fn,'/frequencies','unit');

figure(1); clf; hold on; box on;

plot(x,y)

set(gca,'xtick',xtck,'xticklabel',xtckl,'yminortick','on');
ylabel(['Frequencies (' frequnit ')'])
xlim([0 max(x)])
```

<a name="outfile.phonon_dos.hdf5"></a>
#### 'outfile.phonon_dos.hdf5'

Contains the phonon density of states. The file is self-documented, below is a sample matlab snippet that plots it.

=== "Python"

	```python
	hej
	```

=== "Matlab"
	```matlab
	fn=('outfile.phonon_dos.hdf5');
	unique_atom_labels=strsplit(h5readatt(fn,'/','unique_atom_labels'));
	energy_unit=h5readatt(fn,'/frequencies','unit');
	dos_unit=h5readatt(fn,'/dos','unit')

	omega=h5read(fn,'/frequencies');
	dos=h5read(fn,'/dos');
	dos_per_mode=h5read(fn,'/dos_per_mode');
	dos_per_site=h5read(fn,'/dos_per_site');
	dos_per_unique_atom=h5read(fn,'/dos_per_unique_atom');

	figure(1); clf;

	subplot(1,3,1); hold on; box on;

		title('Phonon density of states')
		for i=1:size(dos_per_mode,2)
			plot(omega,dos)
		end
		set(gca,'xminortick','on','yminortick','on')
		ylabel(['Phonon density of states (' dos_unit ')'])
		xlabel(['Frequency (' energy_unit ')'])

	subplot(1,3,2); hold on; box on;

		title('Phonon density of states per mode')
		for i=1:size(dos_per_mode,2)
			plot(omega,dos_per_mode(:,i))
		end
		set(gca,'xminortick','on','yminortick','on')
		ylabel(['Phonon density of states (' dos_unit ')'])
		xlabel(['Frequency (' energy_unit ')'])

	subplot(1,3,3); hold on; box on;

		title('Phonon density of states per unique atom')
		y=zeros(size(omega));
		for i=1:size(dos_per_unique_atom,2)
			y=y+dos_per_unique_atom(:,i)
			plot(omega,y)
		end

		l=legend(unique_atom_labels);
		set(l,'edgecolor','none','location','northwest')
		set(gca,'xminortick','on','yminortick','on')
		ylabel(['Phonon density of states (' dos_unit ')'])
		xlabel(['Frequency (' energy_unit ')'])
	```

#### `outfile.dispersion_relations`

This file contains a list of $q$ points (in 1/Å) according to the chosen path (default or specified by the user in `infile.qpoints_dispersions`) and the frequencies per mode (in the specified units via `--unit`) for the corresponding points.

<table class='table table-striped'>
<thead><tr>
	<th>Row</th>
	<th>Description</th>
</tr></thead>
<tbody>
<tr>
	<td>1</td>
	<td>
	\( q_1 \qquad \omega_1 \qquad \omega_2 \qquad \ldots \qquad \omega_{3N_a} \)
	</td>
</tr>
<tr>
	<td>2</td>
	<td>
	\( q_2 \qquad \omega_1 \qquad \omega_2 \qquad \ldots \qquad \omega_{3N_a} \)	</td>
</tr>
<tr>
	<td>...</td>
	<td>...</td>
</tr>
</tbody>
</table>

Not that the first column is for plotting purposes only. It serves to ensure that each line segment is scaled properly, since each segment contains a fix number of points.

#### `outfile.group_velocities`

This file contains the norm of the group velocities as a function of q. Format is identical to that of the dispersions:

<table class='table table-striped'>
<thead><tr>
	<th>Row</th>
	<th>Description</th>
</tr></thead>
<tbody>
<tr>
	<td>1</td>
	<td>
	\( q_1 \qquad |v_1| \qquad |v_2| \qquad \ldots \qquad |v_{3N_a}| \)
	</td>
</tr>
<tr>
	<td>2</td>
	<td>
	\( q_2 \qquad |v_1| \qquad |v_2| \qquad \ldots \qquad |v_{3N_a}| \)	</td>
</tr>
<tr>
	<td>...</td>
	<td>...</td>
</tr>
</tbody>
</table>

The units are in km/s.

#### `outfile.mode_gruneisen_parameters`

In case you used `--gruneisen`, the mode Grüneisen parameters will be written, in a format similar to the dispersions and group velocities:

<table class='table table-striped'>
<thead><tr>
	<th>Row</th>
	<th>Description</th>
</tr></thead>
<tbody>
<tr>
	<td>1</td>
	<td>
	$ q_1 \qquad \gamma_1 \qquad \gamma_2 \qquad \ldots \qquad \gamma_{3N_a} $
	</td>
</tr>
<tr>
	<td>2</td>
	<td>
	$ q_2 \qquad \gamma_1 \qquad \gamma_2 \qquad \ldots \qquad \gamma_{3N_a} $	</td>
</tr>
<tr>
	<td>...</td>
	<td>...</td>
</tr>
</tbody>
</table>

The units are in km/s.

<a name="outfile.phonon_dos"></a>
####  `outfile.phonon_dos`

The first column is the list of frequencies and the second column in the density of states for this frequency, given in arbitrary units. Note that integrated density of states is normalized. One can calculate partial density of states (`--projected_dos_mode`); if the user specifies this option, an extra column will be added into the `outfile.phonon_dos` file for each type of atom. With `--projected_dos_site` the same thing happens but with an extra column for each atom in the unit cell.

<table class='table table-striped'>
<thead><tr>
	<th>Row</th>
	<th>Description</th>
</tr></thead>
<tbody>
<tr>
	<td>1</td>
	<td>
	$ \omega_1 \qquad g(\omega_1) \qquad g_1(\omega_1) \qquad \ldots \qquad g_{N}(\omega_1) $
	</td>
</tr>
<tr>
	<td>2</td>
	<td>
	$ \omega_2 \qquad g(\omega_2) \qquad g_1(\omega_2) \qquad \ldots \qquad g_{N}(\omega_2) $
	</td>
</tr>
<tr>
	<td>...</td>
	<td>...</td>
</tr>
</tbody>
</table>

The units are in states per energy unit, depending on choice of `--unit`. The total DOS normalizes to 3N independent of choice of unit.

#### <a name="outfile.free_energy"></a> `outfile.free_energy`

If one chooses the option `--temperature_range` or `--temperature` then this file will display a list of temperatures with corresponding temperatures, vibrational free energies, vibrational entropies and heat capacities.

<table class='table table-striped'>
<thead><tr>
	<th>Row</th>
	<th>Description</th>
</tr></thead>
<tbody>
<tr>
	<td>1</td>
	<td>
	$ T_1 \qquad F_{\textrm{vib}} \qquad S_{\textrm{vib}} \qquad C_v $
	</td>
</tr>
<tr>
	<td>2</td>
	<td>
	$ T_2 \qquad F_{\textrm{vib}} \qquad S_{\textrm{vib}} \qquad C_v $
	</td>
</tr>
<tr>
	<td>...</td>
	<td>...</td>
</tr>
</tbody>
</table>

Temperature is given in K, \( F_{\textrm{vib}} \) in eV/atom, \( S_\textrm{vib} \) in eV/K/atom and heat capacity in eV/K/atom.

#### `outfile.grid_dispersions.hdf5`

Using option `--dumpgrid` writes all phonon properties for a grid in the BZ to an hdf5 file, that is self-documented.

[^Born1998]: Born, M., & Huang, K. (1964). Dynamical theory of crystal lattices. Oxford: Oxford University Press.

[^Hellman2011]: [Hellman, O., Abrikosov, I. A., & Simak, S. I. (2011). Lattice dynamics of anharmonic solids from first principles. Physical Review B, 84(18), 180301.](http://doi.org/10.1103/PhysRevB.84.180301)

[^Hellman2013a]: [Hellman, O., & Abrikosov, I. A. (2013). Temperature-dependent effective third-order interatomic force constants from first principles. Physical Review B, 88(14), 144301.](http://doi.org/10.1103/PhysRevB.88.144301)

[^Hellman2013]: [Hellman, O., Steneteg, P., Abrikosov, I. A., & Simak, S. I. (2013). Temperature dependent effective potential method for accurate free energy calculations of solids. Physical Review B, 87(10), 104111.](http://doi.org/10.1103/PhysRevB.87.104111)

[^Gonze1994]: [Gonze, X., Charlier, J.-C., Allan, D. C. & Teter, M. P. Interatomic force constants from first principles: The case of α-quartz. Phys. Rev. B 50, 13035–13038 (1994).](https://link.aps.org/doi/10.1103/PhysRevB.50.13035)

[^Gonze1997]: [Gonze, X. & Lee, C. Dynamical matrices, Born effective charges, dielectric permittivity tensors, and interatomic force constants from density-functional perturbation theory. Phys. Rev. B 55, 10355–10368 (1997).](http://link.aps.org/doi/10.1103/PhysRevB.55.10355)
