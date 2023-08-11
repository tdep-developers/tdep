author: Olle Hellman
display: none
graph: none
secret:
propname: project confused parrot
propnamelink: <a href="../program/anharmonic_free_energu.html">anharmonic free energy</a>
{!man/project_confused_parrot.md!}

This code calculates the anharmonic Helmholtz free energy. It includes the contributions from baseline shifts, renormalized phonons and higher order terms.

### Free energy

@todo Explain the basics of free energy, entropies and so on.

## Second order terms

The details an notation how to calculate phonon can be found [phonon_dispersion_relations](phonon_dispersion_relations.html), and I follow the same notation here. In the previous section, we treated the vibrations of atoms classically, by solving Newton's equations of motion. Quantum-mechanically, vibrational normal modes can be represented as quasi-particles called phonons, quanta of thermal energy. We note that our normal mode transformation is a sum over eigenfunctions of independent harmonic oscillators. This allows us to write the position and momentum operators in terms of creation and annihilation operators (without loss of generality, we can contract the notation for phonon mode $s$ at wave vector $\mathbf{q}$ to a single index $\lambda$):

$$
\begin{align}
\hat{u}_{i\alpha} = & \sqrt{ \frac{\hbar}{2N m_\alpha} }
\sum_\lambda \frac{\epsilon_\lambda^{i\alpha}}{ \sqrt{ \omega_\lambda} }
e^{i\mathbf{q}\cdot\mathbf{r}_i}
\left( \hat{a}^{\mathstrut}_\lambda + \hat{a}^\dagger_\lambda \right) \\
%
\hat{p}_{i\alpha} = & \sqrt{ \frac{\hbar m_\alpha}{2N} }
\sum_\lambda \sqrt{ \omega_\lambda } \epsilon_\lambda^{i\alpha}
e^{i\mathbf{q}\cdot\mathbf{r}_i-\pi/2}
\left( \hat{a}^{\mathstrut}_\lambda - \hat{a}^\dagger_\lambda \right)
%
\end{align}
$$

and their inverse

$$
\begin{align}
\hat{a}^{\mathstrut}_{\lambda} = & \frac{1}{\sqrt{2N\hbar}}
\sum_{i\alpha} \epsilon_\lambda^{i\alpha}
e^{-i\mathbf{q}\cdot\mathbf{r}_i} \left( \sqrt{m_i \omega_\lambda} \hat{u}_{i\alpha}-i \frac{\hat{p}_{i\alpha}}{ \sqrt{ m_i \omega_\lambda }} \right) \\
%
\hat{a}^\dagger_{\lambda} = & \frac{1}{\sqrt{2N\hbar}}
\sum_{i\alpha} \epsilon_\lambda^{i\alpha}
e^{i\mathbf{q}\cdot\mathbf{r}_i} \left( \sqrt{m_i \omega_\lambda} \hat{u}_{i\alpha}-i \frac{\hat{p}_{i\alpha}}{ \sqrt{ m_i \omega_\lambda }} \right)
\end{align}
$$

In terms of these operators, the vibrational Hamiltonian can be written as

$$
\begin{equation}
\hat{H}=\sum_{\lambda}\hbar\omega_\lambda \left( \hat{a}^\dagger_\lambda \hat{a}^{\mathstrut}_\lambda + \frac{1}{2}\right)\,.
\end{equation}
$$

Since $\hat{a}^\dagger_\lambda \hat{a}^{\mathstrut}_\lambda$ are commutative operators, the Hamiltonian is that of a sum of uncoupled harmonic quantum oscillators, each having the partition function

$$
\begin{equation}
Z_{\lambda}=\sum_{n=0}^{\infty}e^{-\beta (n +\frac{1}{2})\hbar\omega_{\lambda} } =
\frac{ e^{-\beta \hbar\omega_{\lambda}/2 } }{1-e^{-\beta \hbar\omega_{\lambda}}}
\end{equation}
$$

that gives the total partition function

$$
\begin{equation}
Z=\prod_{\lambda} \frac{ e^{-\beta \hbar\omega_{\lambda}/2 } }{1-e^{-\beta \hbar\omega_{\lambda}}}.
\end{equation}
$$

From this we can get the Helmholtz (phonon) free energy:

$$
\begin{equation}
F_{\textrm{ph}}= -k_B T \ln Z = \sum_{\lambda} \frac{\hbar \omega_{\lambda}}{2}+k_B T%
\ln \left( 1- \exp \left( -\frac{\hbar \omega_{\lambda}}{k_B T} \right) \right)
\end{equation}
$$

@todo copy-paste from thesis or something

## Higher order terms

@todo copy-paste from paper

## Baseline renormalization

In the conventional quasiharmonic approximation, the total free energy of the system (not considering any terms pertaining to magnetic or configurational degrees of freedom) can be expressed as

$$
\begin{equation}
	F = F_{\textrm{el}} + F_{\textrm{ph}}\,.
\end{equation}
$$

In the TDEP formalism it is not quite that simple.

#### <a name="sec_tdepthermo"></a> Determining the free energy with TDEP

In the TDEP formalism,[^Hellman2013]<sup>,</sup>[^Hellman2013a]<sup>,</sup>[^Hellman2011] with effective force constants, the phonon quasiparticles are different at each temperature. At fix temperature, they behave just like normal bosons, obeying Bose-Einstein statistics and so on. But changing the temperature will change both the occupation numbers and the states that are occupied. Moreover, in the harmonic approximation the baseline energy (with all atoms at their equilibrium positions) is that of the static lattice. With an effective Hamiltonian this baseline is a free parameter. The baseline shift is illustrated in the diagram below:

<center><img src="../media/explain_U0.png" width="600" /></center>

The density depicts the phase space samples used to fit the effective Hamiltonian. The reference energy, or baseline, has been shifted with respect to zero temperature. This baseline is determined by matching the potential energies of the samples of the Born-Oppenheimer surface, $U_{BO}$, and the potential energy of the TDEP model Hamiltonian:

$$
\begin{equation}
\begin{split}
	\left\langle U_{\textrm{BO}} - U_{\textrm{TDEP}} \right\rangle & =
	\left\langle U_{\textrm{BO}} - U_0 -\frac{1}{2} \sum_{ij} \sum_{\alpha\beta} \Phi_{ij}^{\alpha\beta} u^{\alpha}_i u^{\beta}_j  \right\rangle = 0 \\
	U_0 & = \left\langle U_{\textrm{BO}} - \frac{1}{2} \sum_{ij} \sum_{\alpha\beta} \Phi_{ij}^{\alpha\beta} u^{\alpha}_i u^{\beta}_j \right\rangle
\end{split}
\end{equation}
$$

Intuitively, this can be interpreted as that the TDEP force constants are determined by matching forces between the real and model system, in a similar manner we match the energies as well. The new baseline is conveniently expressed as a shift:

$$
\begin{equation}
	\Delta U = U_0-U_{\textrm{stat}}
\end{equation}
$$

Where $U_{\textrm{stat}}$ is the energy of the perfect lattice at 0K. The free energy is then (excluding configurational entropy, magnetic entropy etc.)

$$
\begin{equation}
	F = F_{\textrm{el}} + F_{\textrm{ph}} + \Delta U_0\,.
\end{equation}
$$

### Why does this code not output entropy?

@todo Refer to Cowley

@todo Simple explanation with partial derivatives

@todo Plot with heat capacity

### Harmonic thermodynamics


The entropy and heat capacity is not accessible from a single simulation. Since both the phonon free energy and baseline shift have non-trivial temperature dependencies, a series of calculations for different temperatures are needed, so that the entropy can be calculated via

$$
\begin{equation}
	S = -\left._{}\frac{dF}{dT}\right|_{V}
\end{equation}
$$

A series of calculations on a volume-temperature grid is required to calculate Gibbs free energy, with pressure explicitly calculated as

$$
\begin{equation}
	P = -\left._{}\frac{dF}{dV}\right|_{T}
\end{equation}
$$

Or something.

### Output files



Using option `--dumpgrid` writes all phonon properties for a grid in the BZ to an hdf5 file, that is self-documented.

[^Born1998]: Born, M., & Huang, K. (1964). Dynamical theory of crystal lattices. Oxford: Oxford University Press.

[^Hellman2011]: [Hellman, O., Abrikosov, I. A., & Simak, S. I. (2011). Lattice dynamics of anharmonic solids from first principles. Physical Review B, 84(18), 180301.](http://doi.org/10.1103/PhysRevB.84.180301)

[^Hellman2013a]: [Hellman, O., & Abrikosov, I. A. (2013). Temperature-dependent effective third-order interatomic force constants from first principles. Physical Review B, 88(14), 144301.](http://doi.org/10.1103/PhysRevB.88.144301)

[^Hellman2013]: [Hellman, O., Steneteg, P., Abrikosov, I. A., & Simak, S. I. (2013). Temperature dependent effective potential method for accurate free energy calculations of solids. Physical Review B, 87(10), 104111.](http://doi.org/10.1103/PhysRevB.87.104111)

[^Gonze1994]: [Gonze, X., Charlier, J.-C., Allan, D. C. & Teter, M. P. Interatomic force constants from first principles: The case of α-quartz. Phys. Rev. B 50, 13035–13038 (1994).](https://link.aps.org/doi/10.1103/PhysRevB.50.13035)

[^Gonze1997]: [Gonze, X. & Lee, C. Dynamical matrices, Born effective charges, dielectric permittivity tensors, and interatomic force constants from density-functional perturbation theory. Phys. Rev. B 55, 10355–10368 (1997).](http://link.aps.org/doi/10.1103/PhysRevB.55.10355)
