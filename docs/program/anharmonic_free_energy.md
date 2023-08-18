
### Short description

End world hunger!

### Command line options:




Optional switches:

* `--temperature value`  
    default value 300  
    Temperature used in the occupation numbers. Should be the same as the temperature the force constants where determined at.

* `--qpoint_grid value#1 value#2 value#3`, `-qg value#1 value#2 value#3`  
    default value 26 26 26  
    Density of q-point mesh for Brillouin zone integrations.

* `--integrationtype value`, `-it value`, value in: `1,2,3`  
    default value 2  
    Type of integration. 1 is Gaussian, 2 adaptive Gaussian and 3 Tetrahedron.

* `--sigma value`  
    default value 1.0  
    Global scaling factor for Gaussian/adaptive Gaussian smearing. The default is determined procedurally, and scaled by this number.

* `--temperature_range value#1 value#2 value#3`  
    default value -1 -1 -1  
    Evaluate thermodynamic phonon properties for a series of temperatures, specify min, max and the number of points.

* `--readiso`  
    default value .false.  
    Read the isotope distribution from file

* `--meshtype value`, value in: `1,2,3,4`  
    default value 1  
    Type of q-point mesh. 1 Is a Monkhorst-Pack mesh, 2 an FFT mesh and 3 my fancy wedge-based mesh with approximately the same density the grid-based meshes. 4 build the commensurate mesh of an approximately cubic supercell.

* `--readqmesh`  
    default value .false.  
    Read the q-point mesh from file. To generate a q-mesh file, see the genkpoints utility.

* `--help`, `-h`  
    Print this help message

* `--version`, `-v`  
    Print version
### Examples

`mpirun anharmonic_free_energy` 

This code calculates the anharmonic Helmholtz free energy. It includes the contributions from baseline shifts, renormalized phonons and higher order terms.

### Longer summary

Per the (quasi)harmonic approximation the free energy is determined as

$$
F(T,V) = U(V)+F_{ph}(T,V)
$$

where $U$ is the energy of the static lattice, and $F_{ph}$ is the free energy of the phonons. The only temperature dependence enters as occupation numbers in the phonon free energy expression. In the TDEP formalism it is only slightly more involved.

As established in [extract_forceconstants](extract_forceconstants.md) the TDEP free energy, to lowest order is given by

$$
F(T,V) = U_0(T,V)+F_{ph}(T,V)
$$

where $F_{ph}(T,V)$ is the free energy of the (effective) phonons, and $U_0$ is a renormalized baseline energy, given by

$$
U_0 = \left\langle
	U - \frac{1}{2}\sum_{ij} u_i\Phi_{ij}u_j -
	\frac{1}{2}\sum_{ij} u_i \Phi^{\text{polar}}_{ij}u_j
	\right\rangle
$$

that is, the potential energy from the molecular dynamics/stochastic sampling minus the potential energy of the force constant model. It's important to remember to also subtract the electrostatic energy since it is included in $F_{ph}$.

#### Higher order corrections

The free energy can be improved by explicitly including the contribution from the higher order terms in the Hamiltonian. The free energy then reads

$$
\begin{equation}
    F = U_0 + F_{\textrm{ph}} + \Delta F^{3\textrm{ph}} + \Delta F^{4\textrm{ph}}
\end{equation}
$$

where $F_{\textrm{ph}}$ is the usual phonon free energy and $U_0$ the renormalized baseline free energy given by

$$
\begin{equation}
    U_0 = \left\langle U -
    \frac{1}{2!}\sum_{\substack{ ij\\ \alpha\beta } }\overset{\textrm{lr}}{\Phi}_{ij}^{\alpha\beta}
u_i^\alpha u_j^\beta -
    \frac{1}{2!}\sum_{\substack{ ij\\ \alpha\beta } }\Phi_{ij}^{\alpha\beta}
u_i^\alpha u_j^\beta -
 \frac{1}{3!}
\sum_{\substack{ijk\\ \alpha\beta\gamma}}\Phi_{ijk}^{\alpha\beta\gamma}
u_i^\alpha u_j^\beta u_k^\gamma -
\frac{1}{4!}
	\sum_{\substack{
	ijkl\\
	\alpha\beta\gamma\delta
	}}
\Phi_{ijkl}^{\alpha\beta\gamma\delta}
u_i^\alpha u_j^\beta u_k^\gamma u_l^\delta
    \right\rangle
\end{equation}
$$

Here it is important to note that we have to subtract the long-ranged polar interactions to avoid double-counting them. The explicit anharmonic contributions are given via[^Leibfried1961][^Cowley1963][^wallace1998thermodynamics]

$$
\begin{equation}%\label{eq:deltaF3}
	\Delta F^{3\textrm{ph}} =
	-6
	\sum_{\lambda\lambda'\lambda''}
	\left|
		\Phi_{\lambda\lambda'\lambda''}
	\right|^2
	%
	\left(
	%\Bigg{\{}
	\frac{3n_{\lambda} n_{\lambda'} + 3n_{\lambda} + 1}
	{(\omega_{\lambda}+\omega_{\lambda'}+\omega_{\lambda''})_p}
	+
	\frac{ 6n_{\lambda} n_{\lambda''} - 3 n_{\lambda} n_{\lambda'} + 3n_{\lambda''}}
	{(\omega_{\lambda}+\omega_{\lambda'}-\omega_{\lambda''})_p}
	\right)
	%\Bigg{\}}
	%
	+9\Phi_{\lambda\bar{\lambda}\lambda''}\Phi_{\lambda'\bar{\lambda}'\bar{\lambda}''}
	\frac{4 n_{\lambda}( n_{\lambda'}+1)+1}
	{(\omega_{\lambda''})_p}\,,
\end{equation}
$$

and

$$
\begin{equation}%\label{eq:deltaF4}
	\Delta F^{4\textrm{ph}} =
	3\sum_{\lambda\lambda'}
	\Phi_{\lambda\bar{\lambda}\lambda'\bar{\lambda}'}(2n_{\lambda}+1)(2n_{\lambda'}+1)
\end{equation}
$$


[^Leibfried1961]: Leibfried, G. & Ludwig, W. (1961) Theory of Anharmonic Effects in Crystals. Solid State Phys 12, 275–444.

[^Cowley1963]: Cowley, R. A. (1963) The lattice dynamics of an anharmonic crystal. Adv Phys 12, 421–480.

[^wallace1998thermodynamics]: Wallace, D.C. Thermodynamics of Crystals. (Dover Publications).

