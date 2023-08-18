
### Short description

Calculates the lattice thermal conductivity from the iterative solution of the phonon Boltzmann equation. In addition, cumulative plots and raw data dumps of intermediate values are available.

### Command line options:




Optional switches:

* `--readiso`  
    default value .false.  
    Read the isotope distribution from `infile.isotopes`. The format is specified [here](../files.md#infile.isotopes).

* `--qpoint_grid value#1 value#2 value#3`, `-qg value#1 value#2 value#3`  
    default value 26 26 26  
    Density of q-point mesh for Brillouin zone integrations.

* `--integrationtype value`, `-it value`, value in: `1,2,3`  
    default value 2  
    Type of integration for the phonon DOS. 1 is Gaussian, 2 adaptive Gaussian and 3 Tetrahedron.

* `--sigma value`  
    default value 1.0  
    Global scaling factor for adaptive Gaussian smearing.

* `--threshold value`  
    default value 4.0  
    Consider a Gaussian distribution to be 0 after this many standard deviations.

* `--readqmesh`  
    default value .false.  
    Read the q-point mesh from file. To generate a q-mesh file, see the genkpoints utility.

* `--temperature value`  
    default value -1  
    Evaluate thermal conductivity at a single temperature.

* `--temperature_range value#1 value#2 value#3`  
    default value 100 300 5  
    Series of temperatures for thermal conductivity. Specify min, max and the number of points.

* `--logtempaxis`  
    default value .false.  
    Space the temperature points logarithmically instead of linearly.

* `--max_mfp value`  
    default value -1  
    Add a limit on the mean free path as an approximation of domain size.

* `--dumpgrid`  
    default value .false.  
    Write files with q-vectors, frequencies, eigenvectors and group velocities for a grid.

* `--noisotope`  
    default value .false.  
    Do not consider isotope scattering.

* `--help`, `-h`  
    Print this help message

* `--version`, `-v`  
    Print version
### Examples

`mpirun thermal_conductivity --temperature 300` 

`mpirun thermal_conductivity -qg 15 15 15 --temperature_range 200 600 50` 

`mpirun thermal_conductivity --integrationtype 2 -qg 30 30 30 --max_mfp 1E-6` 

### Longer summary

Heat transport can be determined by solving the inelastic phonon Boltzmann equation. By applying a temperature gradient $\nabla T_\alpha$ in direction $\alpha$, the heat current is given by the group velocities of phonon mode $\lambda$ and non-equilibrium phonon distribution function $\tilde{n}_\lambda$:[^peierls1955quantum]

$$
\begin{equation}
J_{\alpha}=\frac{1}{V}\sum_\lambda
\hbar \omega_\lambda v_{\lambda\alpha} \tilde{n}_{\lambda\alpha}.
\end{equation}
$$

Assuming the thermal gradient is small, the non-equilibrium distribution function can be linearised as,

$$
\tilde{n}_{\lambda\alpha} \approx n_{\lambda}-
v_{\lambda\alpha}
\tau_{\lambda\alpha}
\frac{d n_{\lambda}}{d T}
\frac{d T}{d \alpha} \, ,
$$

That is a linear deviation from the equilibrium distribution function $n_{\lambda}$. Inserting this into the equation 1, and exploiting the fact that the equilibrium occupation carries no heat, we arrive at,

$$
J_{\alpha}=\frac{1}{V}\sum_{\lambda}
\hbar \omega_{\lambda}
\frac{d n_{\lambda}}{d T}
v_{\lambda\alpha}
v_{\lambda\alpha}
\tau_{\lambda\alpha}
\frac{d T}{d \alpha}.
$$

Utilizing Fourier's law, $J=\kappa \nabla T$, and identifying the phonon heat capacity,

$$
c_{\lambda}=
\hbar \omega_\lambda
\frac{d n_{\lambda}}{d T},
$$

we arrive at,

$$
\kappa_{\alpha\beta}=\frac{1}{V} \sum_{\lambda}
c_{\lambda}
v_{\alpha \lambda}v_{\beta \lambda} \tau_{\beta \lambda},
$$

which can be interpreted as follows: the heat transported by each phonon will depend on how much heat it carries, how fast it travels, and how long it lives. The phonon-phonon induced lifetime can be determined from the self-energy $\Gamma_{\lambda}$. In addition, one must consider the scattering with mass impurities (isotopes), and the boundaries of the sample.

### Lifetimes

With the third order force constants we can calculate the phonon lifetimes needed as input to the thermal conductivity calculations. The lifetime due to phonon-phonon scattering is related to the imaginary part of the phonon self energy ( $\Sigma=\Delta+i\Gamma$ ).

$$
\frac{1}{\tau_{\lambda}}=2 \Gamma_{\lambda},
$$

where $\tau_{\lambda}$ is the lifetime phonon mode $\lambda$, and

$$
\begin{split}
\Gamma_{\lambda}=& \frac{\hbar \pi}{16} % _{\lambda'}
\sum_{\lambda'\lambda''}
\left|\Phi_{\lambda\lambda'\lambda''}\right|^2
\bigl[(n_{\lambda'}+n_{\lambda''}+1)
\delta(\omega_{\lambda}-\omega_{\lambda'}-\omega_{\lambda''}) \\
+ & 2(n_{\lambda'}-n_{\lambda''})
\delta(\omega_{\lambda}-\omega_{\lambda'}+\omega_{\lambda''}) \bigr]
\end{split}
$$

$n_{\lambda}$ is the equilibrium occupation number. The sum is over momentum conserving three-phonon processes, $\textbf{q}+\textbf{q}'+\textbf{q}''=\textbf{G}$, and the deltafunctions in frequency ensure energy conservation. The three-phonon matrix elements are given by

$$
\Phi_{\lambda\lambda'\lambda''} =
\sum_{ijk}
\sum_{\alpha\beta\gamma}
\frac{
\epsilon_{\lambda}^{i \alpha}
\epsilon_{\lambda'}^{j \beta}
\epsilon_{\lambda''}^{k \gamma}
}{
\sqrt{m_{i}m_{j}m_{j}}
\sqrt{
	\omega_{\lambda}
	\omega_{\lambda'}
	\omega_{\lambda''}}
}
\Phi^{\alpha\beta\gamma}_{ijk}
e^{i \mathbf{q}\cdot\mathbf{r}_i + i \mathbf{q}'\cdot\mathbf{r}_j+i \mathbf{q}''\cdot\mathbf{r}_k}
$$

where $m_i$ is the mass of atom $i$, $\epsilon_{\lambda}^{\alpha i}$ is component $\alpha$ of the eigenvector for mode $\lambda$ and atom $i$ and $\textbf{r}_i$ is the lattice vector associated with atom $i$.

Mass disorder, in the form of natural isotope distributions also cause thermal resistance. According to Tamura[^Tamura1983], if the isotopes are randomly distributed on the lattice sites then the strength of the isotope scattering can be given by a mass variance parameter $g$:

$$
g_i=\sum_j c_{i}^j \left(\frac{m_i^j-\bar{m_i}}{\bar{m_i}}\right)^2
$$

where $\bar{m_i}$ is the average isotopic mass( $\bar{m_i}=\sum_j c_i^j m_i^j$ ), $m^j_i$ is the mass of isotope $j$ of atom $i$ and $c^j_i$ is its concentration. The contribution to the imaginary part of the self-energy is

$$
\Gamma^{\textrm{iso}}_{\lambda}=
\frac{\pi}{4} \sum_{\lambda'}
\underbrace{\omega_{\lambda}\omega_{\lambda'} \sum_i g_i \left| \epsilon_{\lambda}^{i \dagger} \epsilon_{\lambda'}^{i} \right|^2}_{\Lambda_{\lambda\lambda'}}
\delta(\omega_{\lambda}-\omega_{\lambda'})
$$

Per default, the isotope distribution will be the natural distribution. In case some other distribution is desired, this can be specified.

Scattering by domain boundaries is implemented as

$$
\Gamma^{\textrm{boundary}}_{\lambda} = \frac{ v_{\lambda} }{2d}
$$

Where $d$ is a characteristic domain size.

### Beyond the relaxation time approximation

So far we have have considered the phonon heat conduction as an elastic process, whereas it is inelastic. This can be treated by iteratively solving the phonon boltzmann equation, formulated in terms of the (linear) deviations from equilibrium occupation numbers.[^peierls1929]<sup>,</sup>[^Omini1996]<sup>,</sup>[^Omini]<sup>,</sup>[^Broido2007]<sup>,</sup>[^Broido2005]

### Phonon scattering rates and the phonon Boltzmann equation

I always found it confusing how you arrived at most of these things. This is something I put together for myself, to clear it up a bit. Please bear in mind that this is not an attempt at a formal derivation whatsoever, just to make it a bit easier to interpret the different terms. There might be an arbitrary number of plusses and minuses and other things missing. Recall the transformation we introduced [earlier](phonon_dispersion_relations.md):

$$
\begin{equation}\label{eq:normalmodetransformation}
\hat{u}_{i\alpha} = \sqrt{ \frac{\hbar}{2N m_\alpha} }
\sum_\lambda \frac{\epsilon_\lambda^{i\alpha}}{ \sqrt{ \omega_\lambda} }
e^{i\mathbf{q}\cdot\mathbf{r}_i}
\left( \hat{a}^{\mathstrut}_\lambda + \hat{a}^\dagger_\lambda \right)
\end{equation}
$$

and consider the three-phonon process where two phonons combine into one:

$$
\begin{equation*}
\begin{split}
\mathbf{q} + \mathbf{q}' + \mathbf{q}'' & = \mathbf{G} \\
\omega + \omega' & = \omega''
\end{split}
\end{equation*}
$$

This process changes the state of the system:

$$
\begin{equation}
\underbrace{\left| \ldots , n_{\lambda},n_{\lambda'},n_{\lambda''} , \ldots \right\rangle}_{\left\vert i \right\rangle}
\rightarrow
\underbrace{\left| \ldots , n_{\lambda}-1,n_{\lambda'}-1,n_{\lambda''}+1, \ldots \right\rangle}_{\left\vert f \right\rangle}
\end{equation}
$$

that is, we lost one phonon at $\lambda$ and one at $\lambda'$, and created a phonon at $\lambda''$.
Mostly out of habit, we sandwich the Hamiltonian between the initial and final states:

$$
\begin{equation}\label{eq:sandwich}
{\left\langle f \middle\vert \hat{H} \middle\vert i \right\rangle} =
{\left\langle f \middle\vert \sum_i \frac{p^2_i}{2m}  +
\frac{1}{2!}\sum_{ij} \sum_{\alpha\beta}\Phi_{ij}^{\alpha\beta}
u_i^\alpha u_j^\beta +\frac{1}{3!}
\sum_{ijk} \sum_{\alpha\beta\gamma}\Phi_{ijk}^{\alpha\beta\gamma}
u_i^\alpha u_j^\beta u_k^\gamma \ldots
\middle\vert i \right\rangle}
\end{equation}
$$

and remember the rules for ladder operators, and that the eigenstates to the quantum harmonic oscillator are orthogonal:

$$
\begin{equation*}
\begin{split}
\hat{a}^\dagger \left\vert n \right\rangle & = \sqrt{n+1} \left\vert n + 1 \right\rangle \\
\hat{a} \left\vert n \right\rangle & = \sqrt{n} \left\vert n -1 \right\rangle \\
\left\langle i \middle\vert j \right\rangle & = \delta_{ij}
\end{split}
\end{equation*}
$$

Inserting eq \ref{eq:normalmodetransformation} into \ref{eq:sandwich} (and realising that the kinetic energy part and the second order part disappears), we end up with a pretty large expression, that we will deal with in steps, first identify

$$
\begin{equation}\label{eq:uprod}
\begin{split}
u^\alpha_{i}u^\beta_{j}u^\gamma_{k} & =
%
\left(\frac{\hbar}{2N}\right)^{3/2} \frac{1}{\sqrt{m_{i}m_{j}m_{k}}}
\sum_{\lambda\lambda'\lambda''}
\frac{
\epsilon_{\lambda}^{i \alpha}
\epsilon_{\lambda'}^{j \beta}
\epsilon_{\lambda''}^{k \gamma}
}{
\sqrt{
	\omega_{\lambda}
	\omega_{\lambda'}
	\omega_{\lambda''}}
}
e^{i \mathbf{q}\cdot\mathbf{r}_i + i \mathbf{q}'\cdot\mathbf{r}_j+i \mathbf{q}''\cdot\mathbf{r}_k}
 \left(a_{\lambda}+a_{\lambda}^\dagger \right)
\left(a_{\lambda'}+a_{\lambda'}^\dagger \right)
\left(a_{\lambda''}+a_{\lambda''}^\dagger \right)
\end{split}
\end{equation}
$$

as well as

$$
\begin{equation}
\begin{split}
& \sum_{\lambda\lambda'\lambda''}
\left\langle f \middle\vert
\left(a_{\lambda}+a_{\lambda}^\dagger \right)
\left(a_{\lambda'}+a_{\lambda'}^\dagger \right)
\left(a_{\lambda''}+a_{\lambda''}^\dagger \right)
\middle\vert i \right\rangle = \\
= & \sum_{\lambda\lambda'\lambda''} \left\langle f \middle\vert
\hat{a}_{\lambda}   \hat{a}_{\lambda'}   \hat{a}_{\lambda''} + \hat{a}_{\lambda}   \hat{a}_{\lambda'}   \hat{a}^{\dagger}_{\lambda''} + \hat{a}_{\lambda}   \hat{a}^{\dagger}_{\lambda'}   \hat{a}_{\lambda''} + \hat{a}_{\lambda}   \hat{a}^{\dagger}_{\lambda'}   \hat{a}^{\dagger}_{\lambda''} + \hat{a}^{\dagger}_{\lambda}   \hat{a}_{\lambda'}   \hat{a}_{\lambda''} + \hat{a}^{\dagger}_{\lambda}   \hat{a}_{\lambda'}   \hat{a}^{\dagger}_{\lambda''} + \hat{a}^{\dagger}_{\lambda}   \hat{a}^{\dagger}_{\lambda'}   \hat{a}_{\lambda''} + \hat{a}^{\dagger}_{\lambda}   \hat{a}^{\dagger}_{\lambda'}   \hat{a}^{\dagger}_{\lambda''}
\middle\vert i \right\rangle = \\
= & \sum_{\lambda\lambda'\lambda''} \left\langle f \middle\vert
a_{\lambda}a_{\lambda'}a^\dagger_{\lambda''}
\middle\vert i \right\rangle
 = 3 \sqrt{n_{\lambda}n_{\lambda'}(n_{\lambda''}+1)}
\end{split}
\end{equation}
$$

where the factor 3 comes from the multiplicity, to get at

$$
\begin{equation}
{\left\langle f \middle\vert \hat{H}_3 \middle\vert i \right\rangle} =
\frac{1}{2}
\sum_{ijk} \sum_{\alpha\beta\gamma}\Phi_{ijk}^{\alpha\beta\gamma}
\sqrt{n_{\lambda}n_{\lambda'}(n_{\lambda''}+1)}
%
\left(\frac{\hbar}{2N}\right)^{3/2}
\frac{
\epsilon_{\lambda}^{i \alpha}
\epsilon_{\lambda'}^{j \beta}
\epsilon_{\lambda''}^{k \gamma}
}{
\sqrt{m_{i}m_{j}m_{j}}
\sqrt{
	\omega_{\lambda}
	\omega_{\lambda'}
	\omega_{\lambda''}}
}
e^{i \mathbf{q}\cdot\mathbf{r}_i + i \mathbf{q}'\cdot\mathbf{r}_j+i \mathbf{q}''\cdot\mathbf{r}_k}
\end{equation}
$$

The initial factor 1/2 is the multiplicity cancelled by the 3! from the Hamiltonian. Here, as it happens, we can identify the three-phonon matrix elements and simplify a little bit more

$$
\begin{equation}
{\left\langle f \middle\vert \hat{H}_3 \middle\vert i \right\rangle} =
\frac{1}{2}
\sqrt{n_{\lambda}n_{\lambda'}(n_{\lambda''}+1)}
\left(\frac{\hbar}{2N}\right)^{3/2}
\Phi_{\lambda\lambda'\lambda''}
\end{equation}
$$

The probability of this particular three-phonon process can be estimated via the Fermi golden rule:

$$
\begin{equation}
\begin{split}
P_{\lambda\lambda'\rightarrow\lambda''} & =\frac{2\pi}{\hbar}
\left|{\left\langle f \middle\vert \hat{H}_3 \middle\vert i \right\rangle}\right|^2
\delta(E_f-E_i) =
\frac{\hbar^2\pi}{16N}
n_{\lambda}n_{\lambda'}(n_{\lambda''}+1)
\left| \Phi_{\lambda\lambda'\lambda''} \right|^2
\delta(E_f-E_i)
\end{split}
\end{equation}
$$

With near identical reasoning, we can also arrive at

$$
\begin{equation}\label{pplus}
P_{\lambda\rightarrow\lambda'\lambda''} =
\frac{\hbar^2\pi}{16N}
n_{\lambda}(n_{\lambda'}+1)(n_{\lambda''}+1)
\left| \Phi_{\lambda\lambda'\lambda''} \right|^2
\delta(E_f-E_i)
\end{equation}
$$

for the other kind of three-phonon processes, and

$$
\begin{equation}\label{pminus}
P_{\lambda\rightarrow\lambda'} =\frac{2\pi}{\hbar}\left|\langle f | H^{\textrm{iso}} | i \rangle \right|^2\delta(E_f-E_i) =
\frac{\pi\hbar}{2N}  n_{\lambda}(n_{\lambda'}+1) \Lambda_{\lambda\lambda'}\delta(E_f-E_i)
\end{equation}
$$

for the isotope scattering. I leave those derivations as an exercise. The phonon Boltzmann equation is stated as:

$$
\begin{equation}\label{eq:pbe}
\frac{\partial \tilde{n}_\lambda}{\partial T} \mathbf{v}_\lambda \cdot \nabla T =
\left. \frac{\partial \tilde{n}_\lambda }{\partial t} \right|_{\mathrm{coll}}
\end{equation}
$$

Where $\tilde{n}$ is the non-equilibrium occupation number. This is ridiculously complicated. To make life easier, we only consider the terms we outlined above as possible collisions. Gathering all possible events that involve mode $\lambda$ we get

$$
\begin{equation}\label{manyprob}
\begin{split}
\left. \frac{\partial n_{\lambda}}{ \partial t} \right|_{\mathrm{coll}}
= &  \sum_{\lambda'}
\left( P_{\lambda\rightarrow\lambda'}-P_{\lambda'\rightarrow\lambda } \right) +
\sum_{\lambda'\lambda''}
- P_{\lambda  \rightarrow \lambda' \lambda'' }
- P_{\lambda  \rightarrow \lambda''\lambda'  }
+ P_{\lambda' \rightarrow \lambda  \lambda'' }
+ P_{\lambda' \rightarrow \lambda''\lambda   }
+ P_{\lambda''\rightarrow \lambda  \lambda'  }
+ P_{\lambda''\rightarrow \lambda' \lambda   } \\
& - P_{\lambda   \lambda' \rightarrow \lambda'' }
- P_{\lambda   \lambda'' \rightarrow \lambda'  }
- P_{\lambda'  \lambda   \rightarrow \lambda'' }
+ P_{\lambda'  \lambda'' \rightarrow \lambda   }
- P_{\lambda'' \lambda   \rightarrow \lambda'  }
+ P_{\lambda'' \lambda'  \rightarrow \lambda   }
\end{split}
\end{equation}
$$

Which does not seem to make life easier. To make it slightly worse, we insert \ref{pplus} and \ref{pminus} into this, and at the same time say that the non-equilibrium distribution functions are the equilibrium distributions, plus a (small) deviation:

$$
\begin{equation}
\tilde{n}_{\lambda}\approx n_{\lambda}+\epsilon_{\lambda}
\end{equation}
$$

After some [hard work](https://reference.wolfram.com/language/ref/FullSimplify.html), and discarding terms of $\epsilon^2$ and higher, we get

$$
\begin{equation}
\begin{split}
\left. \frac{\partial n_{\lambda}}{ \partial t} \right|_{\mathrm{coll}}
= & \sum_{\lambda'\lambda''}
\frac{\hbar\pi}{8N}
\left| \Phi_{\lambda\lambda'\lambda''} \right|^2 \Big(
\left[
-n_{\lambda} \epsilon_{\lambda'} + n_{\lambda''} (\epsilon_{\lambda} + \epsilon_{\lambda'}) + \epsilon_{\lambda''} + n_{\lambda} \epsilon_{\lambda''} + n_{\lambda'} (-\epsilon_{\lambda} + \epsilon_{\lambda''})
\right]\delta(\omega_{\lambda}+\omega_{\lambda'}-\omega_{\lambda''}) + \\
& \left[
\epsilon_{\lambda'} + n_{\lambda} \epsilon_{\lambda'} + n_{\lambda''} (-\epsilon_{\lambda} + \epsilon_{\lambda'}) - n_{\lambda} \epsilon_{\lambda''} +
  n_{\lambda'} (\epsilon_{\lambda} + \epsilon_{\lambda''} )
\right]\delta(\omega_{\lambda}-\omega_{\lambda'}+\omega_{\lambda''}) - \\
& \left[(1 + n_{\lambda'} + n_{\lambda''})\epsilon_{\lambda} - n_{\lambda''}\epsilon_{\lambda''} - n_{\lambda'} \epsilon_{\lambda''} + n_{\lambda} (\epsilon_{\lambda'} + \epsilon_{\lambda''} )\right]
\delta(\omega_{\lambda}-\omega_{\lambda'}-\omega_{\lambda''}) \Big)
\end{split}
\end{equation}
$$

Which does not seem like a lot of help. If we make another substitution, and say that the deviation from equilibrium behaves sort of like the equilibrium (with no loss of generality, just to make life easier):

$$
\begin{equation}
\epsilon_{\lambda} =
\frac{\partial n_{\lambda} }{\partial \omega_\lambda}
\frac{k_B T}{\hbar} \zeta_{\lambda}=-n_{\lambda}(n_{\lambda}+1) \zeta_{\lambda}
\end{equation}
$$

Inserting this, and more tedious algebra, we get

$$
\begin{equation}
\begin{split}
\left. \frac{\partial n_{\lambda}}{ \partial t} \right|_{\mathrm{coll}}
=& \frac{\hbar\pi}{4N}
\sum_{\lambda'\lambda''}
\left| \Phi_{\lambda\lambda'\lambda''} \right|^2 \Big(
n_{\lambda} n_{\lambda'} (n_{\lambda''}+1) \delta(\omega_{\lambda}+\omega_{\lambda'}-\omega_{\lambda''} )
\left( \zeta_{\lambda} + \zeta_{\lambda'} - \zeta_{\lambda''} \right) + \\
& \frac{1}{2} n_{\lambda} (n_{\lambda'}+1) (n_{\lambda''}+1) \delta(\omega_{\lambda}-\omega_{\lambda'}-\omega_{\lambda''})
\left( \zeta_{\lambda} - \zeta_{\lambda'} -\zeta_{\lambda''} \right) \Big)
\end{split}
\end{equation}
$$

If we add the isotope term again, that I forgot at some point between the beginning and here, we can rearrange this in terms of scattering rates that should look familiar (using strange relations for occupation numbers that only hold when the deltafunctions in energy are satisfied):

$$
\begin{equation}
\left. \frac{\partial n_{\lambda}}{ \partial t} \right|_{\mathrm{coll}} =
\sum_{\lambda'\lambda''}
\tilde{P}^{+}_{\lambda\lambda'\lambda''}
\left( \zeta_{\lambda}+\zeta_{\lambda'}-\zeta_{\lambda''} \right)
+\frac{1}{2}\tilde{P}^{-}_{\lambda\lambda'\lambda''}
\left( \zeta_{\lambda}-\zeta_{\lambda'}-\zeta_{\lambda''} \right)+
\sum_{\lambda'}
\tilde{P}^\textrm{iso}_{\lambda\lambda'} \left( \zeta_{\lambda}-\zeta_{\lambda'} \right)
\end{equation}
$$

where

$$
\begin{align}
\tilde{P}^{+}_{\lambda\lambda'\lambda''}&=
\frac{\hbar \pi}{4 N}
n_{\lambda}n_{\lambda'}(n_{\lambda''}+1)\left|\Phi_{\lambda\lambda'\lambda''}\right|^2
\delta(\omega_{\lambda}+\omega_{\lambda'}-\omega_{\lambda''})
\\
\tilde{P}^{-}_{\lambda\lambda'\lambda''}&=
\frac{\hbar \pi}{4 N}
n_{\lambda}(n_{\lambda'}+1)(n_{\lambda''}+1)\left|\Phi_{\lambda\lambda'\lambda''}\right|^2
\delta(\omega_{\lambda}-\omega_{\lambda'}-\omega_{\lambda''})
\\
\tilde{P}^\textrm{iso}_{\lambda\lambda'} &=
\frac{\pi}{2N} n_{\lambda}(n_{\lambda'}+1) \Lambda_{\lambda\lambda'}
\delta(\omega_{\lambda}-\omega_{\lambda})
\end{align}
$$

What we have done here is to rearrange the transition propabilities to scattering rates. If we let

$$
\begin{equation}
\zeta_{\lambda}=\frac{\hbar}{k_B T} \mathbf{F}_{\lambda} \cdot \nabla T
\end{equation}
$$

and combine everything we end up with

$$
\begin{equation}
\begin{split}
-\frac{\omega_{\lambda}}{T}n_{\lambda}(n_{\lambda}+1)\mathbf{v}_{\lambda} \cdot \nabla T = &
 \sum_{\lambda'}
\tilde{P}^\textrm{iso}_{\lambda\lambda'}
\left(\mathbf{F}_{\lambda}-\mathbf{F}_{\lambda'}\right)\cdot\nabla T +
\sum_{\lambda'\lambda''}
\tilde{P}^{+}_{\lambda\lambda'\lambda''}
\left(\mathbf{F}_{\lambda}+\mathbf{F}_{\lambda'}-\mathbf{F}_{\lambda''}\right)\cdot\nabla T+
\tilde{P}^{-}_{\lambda\lambda'\lambda''}
\left(\mathbf{F}_{\lambda}-\mathbf{F}_{\lambda'}-\mathbf{F}_{\lambda''}\right)\cdot\nabla T =
\\ = &
\mathbf{F}_{\lambda}\cdot\nabla T
\left(
\sum_{\lambda'}
\tilde{P}^\textrm{iso}_{\lambda\lambda'}
+
\sum_{\lambda'\lambda''}
\tilde{P}^{+}_{\lambda\lambda'\lambda''}+
\frac{1}{2}\tilde{P}^{-}_{\lambda\lambda'\lambda''}
\right)- \\
& - \sum_{\lambda'}
\tilde{P}^\textrm{iso}_{\lambda\lambda'}\mathbf{F}_{\lambda'}\cdot\nabla T +
\sum_{\lambda'\lambda''}
\tilde{P}^{+}_{\lambda\lambda'\lambda''}
\left(\mathbf{F}_{\lambda'}-\mathbf{F}_{\lambda''}\right)\cdot\nabla T-
\frac{1}{2}\tilde{P}^{-}_{\lambda\lambda'\lambda''}
\left(\mathbf{F}_{\lambda'}-\mathbf{F}_{\lambda''}\right)\cdot\nabla T
\end{split}
\end{equation}
$$

Where we can identify

$$
\begin{equation}
Q_{\lambda}=\sum_{\lambda'}
\tilde{P}^\textrm{iso}_{\lambda\lambda'}
+
\sum_{\lambda'\lambda''}
\tilde{P}^{+}_{\lambda\lambda'\lambda''}+
\frac{1}{2}\tilde{P}^{-}_{\lambda\lambda'\lambda''}
\end{equation}
$$

And rearrange terms

$$
\begin{equation}
\mathbf{F}_{\lambda}=
\frac{\omega_{\lambda} \bar{n}_{\lambda}(\bar{n}_{\lambda}+1)\mathbf{v}_{\lambda} }{T Q_{\lambda}}
+
\frac{1}{Q_{\lambda}}\left[
\sum_{\mathbf{q}'\mathbf{q}''}\sum_{s's''}
\tilde{P}^{+}_{\lambda\lambda'\lambda''}
\left( \mathbf{F}_{\lambda'}-\mathbf{F}_{\lambda''} \right)-
\frac{1}{2}\tilde{P}^{-}_{\lambda\lambda'\lambda''}
\left( \mathbf{F}_{\lambda'}-\mathbf{F}_{\lambda''} \right)
\right]
\end{equation}
$$

And we have a set of equations for $F$ that we can solve self-consistently. Previously, we used the imaginary part of the self-energy to get a phonon lifetime. What we got here, from Fermi golden rule, is related:

$$
\sum_{\lambda'} \tilde{P}^\textrm{iso}_{\lambda\lambda'} =
\frac{\pi}{2N} n_{\lambda}(n_{\lambda}+1) \sum_{\lambda'}  \Lambda_{\lambda\lambda'}
\delta(\omega_{\lambda}-\omega_{\lambda}) = 2 n_{\lambda}(n_{\lambda}+1) \Gamma^{\textrm{iso}}_{\lambda}
$$

This can also be done for the three-phonon terms:

$$
\begin{equation}
\begin{split}
\sum_{\lambda'\lambda''} \tilde{P}^{+}_{\lambda\lambda'\lambda''}+
\frac{1}{2}\tilde{P}^{-}_{\lambda\lambda'\lambda''} & =
\frac{\hbar \pi}{8 N}
\sum_{\lambda'\lambda''} \left|\Phi_{\lambda\lambda'\lambda''}\right|^2
\left[
n_{\lambda}(n_{\lambda'}+1)(n_{\lambda''}+1) \delta(\omega_{\lambda}-\omega_{\lambda'}-\omega_{\lambda''})+
2n_{\lambda}n_{\lambda'}(n_{\lambda''}+1)\delta(\omega_{\lambda}+\omega_{\lambda'}-\omega_{\lambda''})
\right] \\
& = n_{\lambda}(n_{\lambda}+1) \frac{\hbar \pi}{8 N}
\sum_{\lambda'\lambda''} \left|\Phi_{\lambda\lambda'\lambda''}\right|^2
\left[
\frac{n_{\lambda}(n_{\lambda'}+1)(n_{\lambda''}+1)}{n_{\lambda}(n_{\lambda}+1)}
\delta(\omega_{\lambda}-\omega_{\lambda'}-\omega_{\lambda''})
+
\frac{2n_{\lambda}n_{\lambda'}(n_{\lambda''}+1)}{n_{\lambda}(n_{\lambda}+1)}
\delta(\omega_{\lambda}+\omega_{\lambda'}-\omega_{\lambda''})
\right] \\
& =
n_{\lambda}(n_{\lambda}+1) \frac{\hbar \pi}{8 N}
\sum_{\lambda'\lambda''} \left|\Phi_{\lambda\lambda'\lambda''}\right|^2
\left[
(n_{\lambda'}+n_{\lambda''}+1)
\delta(\omega_{\lambda}-\omega_{\lambda'}-\omega_{\lambda''})
+
(n_{\lambda'}-n_{\lambda''})
\delta(\omega_{\lambda}+\omega_{\lambda'}-\omega_{\lambda''})
\right] \\
& = 2 n_{\lambda}(n_{\lambda}+1) \Gamma_{\lambda}
\end{split}
\end{equation}
$$

Where the second to last step seems a little impossible, but with $\hbar\omega/k_BT = x$, you get

$$
\begin{equation}
\frac{  n_{\lambda}(n_{\lambda'}+1)(n_{\lambda''}+1)  }{ n_{\lambda}(n_{\lambda}+1) } -
\left( n_{\lambda'} + n_{\lambda''} + 1 \right)
=
\frac{
1-\exp[x'+x''-x]
}{
\left( \exp[x'] -1 \right) \left( \exp[x''] -1 \right)
}
\end{equation}
$$

which comes out to 0 when $x=x'+x''$, which the deltafunction ensures. In the same way

$$
\begin{equation}
\frac{  n_{\lambda}n_{\lambda'}(n_{\lambda''}+1)  }{ n_{\lambda}(n_{\lambda}+1) } -
\left( n_{\lambda'} - n_{\lambda''} \right)
=
\frac{
\exp[-x]\left(\exp[x+x']-\exp[x''] \right)
}{
\left( \exp[x'] -1 \right) \left( \exp[x''] -1 \right)
}
\end{equation}
$$

comes out to 0 when $x''=x+x'$. We can directly relate the relaxation time lifetime

$$
\begin{equation}
\tau_{\lambda} = \frac{1}{2\Gamma_{\lambda}} = \frac{ n_{\lambda}(n_{\lambda}+1) }{Q_{\lambda}}
\end{equation}
$$

to an initial guess

$$
\mathbf{F}^0_{\lambda} =
\frac{\tau_{\lambda} \omega_{\lambda} \mathbf{v}_{\lambda} }{T}
$$

and iteratively solve

$$
\begin{equation}
\mathbf{F}^{i+1}_{\lambda}=
\mathbf{F}^0_{\lambda}
+
\frac{1}{Q_{\lambda}}\left[
\sum_{\lambda'\lambda''}
\tilde{P}^{+}_{\lambda\lambda'\lambda''}
\left( \mathbf{F}^{i}_{\lambda'}-\mathbf{F}^{i}_{\lambda''} \right)-
\frac{1}{2}\tilde{P}^{-}_{\lambda\lambda'\lambda''}
\left( \mathbf{F}^{i}_{\lambda'}-\mathbf{F}^{i}_{\lambda''} \right)
\right]
\end{equation}
$$

to arrive at the non-equilibrium distributions. The thermal conductivity tensor is then given as

$$
\begin{equation}
\kappa_{\alpha\beta} =
\frac{1}{V}
\sum_{\lambda}
\frac{T c_{\lambda} v_{\lambda}^\alpha F_{\lambda}^\beta}{\omega_{\lambda}}
\end{equation}
$$

### Cumulative kappa

@todo Check code snippets

@todo Spectral kappa, links to things.

Experimentally, the cumulative thermal conductivity with respect to phonon mean free path,

$$
l_{\lambda} = \left| v_{\lambda} \right| \tau_{\lambda} \,,
$$

can be measured.[^Minnich2012] The cumulative thermal conductivity can then be computed as a sum of the fraction of heat that is carried by phonons with mean free paths smaller than $l$:

$$
\kappa_{\alpha\beta}^{\textrm{acc}}(l)=
\frac{1}{V} \sum_{\lambda}
C_{\lambda} v^{\alpha}_{\lambda} v^{\beta}_{\lambda} \tau_{\lambda} \Theta(l- l_{\lambda} ) \,,
$$

where $\Theta$ is the Heaviside step function.

One can also define a spectral thermal conductivity as

$$
\kappa_{\alpha\beta}(\omega)=
\frac{1}{V} \sum_{\lambda}
C_{\lambda} v^{\alpha}_{\lambda} v^{\beta}_{\lambda} \tau_{\lambda} \delta(\omega- \omega_{\lambda} )
$$

which is a measure which frequencies contribute most to thermal transport.

### Thin film scattering

Constrained geometries will incur additional scattering from domain boundaries. For a thin film (thin, but thick enough that the interior of the film is accurately described by bulk phonons) one can estimate the suppression due to film thinkness.[^Minnich2015] Assyming the cross-plane direction of the film is in the $y$-direction, and the thermal gradient is applied in the $z$-direction, the in-plane thermal conductivity $\kappa_{zz}$ is supressed as:

$$
\kappa_{zz}(d)=A+B+C,
$$

where

$$
\begin{split}
x_{\lambda} = & \frac{\hbar\omega_{\lambda}}{V}
\frac{\partial n_{\lambda}} {\partial T}
v_{\lambda}^z l_{\lambda}^z
  \\
A = & -\frac{1}{d} \sum_{v_y>0}
x_{\lambda}
\left(  -l^{y}_{\lambda} \exp\left[\frac{d}{l^{y}_{\lambda}}\right]+l^{y}_{\lambda}-d  \right) \\
B = & -\frac{1}{d} \sum_{v_y<0}
x_{\lambda}
\left(
l^{y}_{\lambda}
\exp\left[  -\frac{d}{l^{y}_{\lambda}}  \right]
-l^{y}_{\lambda}-d
\right) \\
C = & \sum_{v_y=0} x_{\lambda}
\end{split}
$$

where $v_y$ and $v_z$ are the components of the phonon group velocity along the $y$ and $z$ directions, $\tau_{\lambda}$ is the phonon relaxation time. $l^{y}_{\lambda}$ is the $y$ component of the MFP and $d$ is the thickness of the film in $y$-direction.

### Input files

These files are necesarry:

* [infile.ucposcar](../files.md#infile.ucposcar)
* [infile.forceconstant](extract_forceconstants.md#infile.forceconstant)
* [infile.forceconstant_thirdorder](extract_forceconstants.md#infile.forceconstant_thirdorder)

and these are optional:

* [infile.isotopes](../files.md#infile.isotopes) (for non-natural isotope distribution)

### Output files

Depending on options, the set of output files may differ. We start with the basic files that are written after running this code.

#### `outfile.thermal_conductivity`

This file contains components of the thermal conductivity tensor $\kappa_{\alpha \beta}$ for each temperature.

<table class='table table-striped'>
<thead><tr>
	<th>Row</th>
	<th>Description</th>
</tr></thead>
<tbody>
<tr>
	<td>1</td>
	<td>
	\( T_1 \qquad \kappa_{xx} \quad \kappa_{yy} \quad \kappa_{zz} \quad \kappa_{xz} \quad	\kappa_{yz} \quad	\kappa_{xy} \quad \kappa_{zx} \quad \kappa_{zy} \quad	\kappa_{yx} \)
	</td>
</tr>
<tr>
	<td>2</td>
	<td>
	\( T_2 \qquad \kappa_{xx} \quad \kappa_{yy} \quad \kappa_{zz} \quad \kappa_{xz} \quad	\kappa_{yz} \quad	\kappa_{xy} \quad \kappa_{zx} \quad \kappa_{zy} \quad	\kappa_{yx} \)
	</td>
</tr>
<tr>
	<td>...</td>
	<td>...</td>
</tr>
</tbody>
</table>

#### `outfile.cumulative_kappa.hdf5`

This file is self-explainatory. It contains the different cumulative plots described above, at a series of temperatures. Below is a matlab snippet that plots part of the output.

```matlab
figure(1); clf; hold on; box on;

% filename
fn='outfile.cumulative_kappa.hdf5';
% which temperature?
t=1;

subplot(1,3,1); hold on; box on;

    % read in cumulative kappa vs mean free path from file
    x=h5read(fn,['/temperature_' num2str(t) '/mean_free_path_axis']);
    xunit=h5readatt(fn,['/temperature_' num2str(t) '/mean_free_path_axis'],'unit');
    y=h5read(fn,['/temperature_' num2str(t) '/cumulative_kappa_vs_mean_free_path_total']);
    % projections to modes and/or atoms
    z=h5read(fn,['/temperature_' num2str(t) '/cumulative_kappa_vs_mean_free_path_per_atom']);
    %z=h5read(fn,['/temperature_' num2str(t) '/cumulative_kappa_vs_mean_free_path_per_mode']);

    yunit=h5readatt(fn,['/temperature_' num2str(t) '/cumulative_kappa_vs_mean_free_path_total'],'unit');

    % plot
    plot(x,y)
    plot(x,z)

    % set a legend
    lgd{1}='Total';
    for i=1:size(z,2)
        lgd{i+1}=['Atom ' num2str(i)];
    end
    l=legend(lgd);
    set(l,'edgecolor','none','location','northwest');

    % some titles
    title('Cumulative kappa vs mean free path');
    ylabel(['Cumulative \kappa (' yunit ')']);
    xlabel(['Mean free path (' xunit ')']);

    % get some reasonable ranges
    minx=x(max(find(y<max(y*1E-3))));
    maxx=x(min(find(y>max(y*0.9999))))*2;
    xlim([minx maxx]);
    set(gca,'xscale','log','yminortick','on');

subplot(1,3,2); hold on; box on;

    % read in spectral kappa vs frequency
    x=h5read(fn,['/temperature_' num2str(t) '/frequency_axis']);
    y=h5read(fn,['/temperature_' num2str(t) '/spectral_kappa_vs_frequency_total']);
    z=h5read(fn,['/temperature_' num2str(t) '/spectral_kappa_vs_frequency_per_mode']);
	%z=h5read(fn,['/temperature_' num2str(t) '/spectral_kappa_vs_frequency_per_atom']);

	plot(x,y)
    plot(x,z)

    set(gca,'xminortick','on','yminortick','on')
    xlabel('Frequency (THz)')
    ylabel('Spectral \kappa (W/m/K/THz)')
	title('Spectral kappa vs frequency')

subplot(1,3,3); hold on; box on;

    % read in cumulative kappa vs mean free path from file
    x=h5read(fn,['/temperature_' num2str(t) '/boundary_scattering_lengths']);
    y=h5read(fn,['/temperature_' num2str(t) '/boundary_scattering_kappa']);
    % grab only kxx
    y=squeeze(y(1,1,:));
	plot(x,y)

    set(gca,'xscale','log','yminortick','on')
    xlabel('Domain size (m)')
    ylabel('Kappa (W/mK)')
	title('Kappa vs boundary scattering')

    % get a reasonable range in x
    minx=max(x(find(y<max(y*1E-2))))
    maxx=x(min(find(y>max(y*0.9999))))*2
    xlim([minx maxx])
```

#### `outfile.grid_thermal_conductivity.hdf5`

Option `--dumpgrid` produces this self-explainatory file. It will not get written if you use more than one temperature, the reason is that this file can get uncomfortably large, nearly all quantities on the full q-grid are written. Below is a matlab snippet that plots a subset:

```matlab

% file to read from
fn='outfile.grid_thermal_conductivity.hdf5';
% convert units to THz from Hz?
toTHz=1/1E12/2/pi;

figure(1); clf; hold on;

subplot(1,3,1); hold on; box on;

    x=h5read(fn,'/frequencies');
    y=h5read(fn,'/linewidths');

    for i=1:size(x,1)
        plot(x(i,:)*toTHz,y(i,:)*toTHz,'marker','.','linestyle','none','markersize',8)
    end
    set(gca,'xminortick','on','yminortick','on')
    xlabel('Frequency (THz)')
    ylabel('Linewidth (THz)')

subplot(1,3,2); hold on; box on;

    x=h5read(fn,'/frequencies');
    y=h5read(fn,'/lifetimes');

    for i=1:size(x,1)
        plot(x(i,:)*toTHz,y(i,:),'marker','.','linestyle','none','markersize',8)
    end
    set(gca,'yscale','log','xminortick','on')
    xlabel('Frequency (THz)')
    ylabel('Lifetime (s)')

subplot(1,3,3); hold on; box on;

    x=h5read(fn,'/frequencies');
    y=h5read(fn,'/mean_free_paths');

    for i=1:size(x,1)
        plot(x(i,:)*toTHz,y(i,:),'marker','.','linestyle','none','markersize',8)
    end
    set(gca,'yscale','log','xminortick','on')
    xlabel('Frequency (THz)')
    ylabel('Mean free paths (m)')

```

[^peierls1929]: Peierls, R. E. (1929). Annalen der Physik, 3, 1055

[^peierls1955quantum]: [Peierls, R. E. (1955). Quantum Theory of Solids. Clarendon Press.](https://books.google.com/books?id=WvPcBUsSJBAC)

[^Minnich2012]: [Minnich, A. J. (2012). Determining phonon mean free paths from observations of quasiballistic thermal transport. Physical Review Letters, 109(20), 1–5.](http://doi.org/10.1103/PhysRevLett.109.205901)

[^Minnich2015]: [Minnich, A. J. (2015). Thermal phonon boundary scattering in anisotropic thin films. Applied Physics Letters, 107(18), 8–11.](http://doi.org/10.1063/1.4935160)

[^Tamura1983]: [Tamura, S. (1983). Isotope scattering of dispersive phonons in Ge. Physical Review B, 27(2), 858–866.](http://doi.org/10.1103/PhysRevB.27.858)

[^Omini1996]: [Omini, M., & Sparavigna, A. (1996). Beyond the isotropic-model approximation in the theory of thermal conductivity. Physical Review B, 53(14), 9064–9073.](http://doi.org/10.1103/PhysRevB.53.9064)

[^Omini]: [Omini, M., & Sparavigna, A. (1997). Heat transport in dielectric solids with diamond structure. Nuovo Cimento Della Societa Italiana Di Fisica D, 19D, 1537–63.](http://www.sif.it/riviste/ncd/econtents/1997/019/10/article/5)

[^Broido2007]: [Broido, D. A., Malorny, M., Birner, G., Mingo, N., & Stewart, D. A. (2007). Intrinsic lattice thermal conductivity of semiconductors from first principles. Applied Physics Letters, 91(23), 231922.](http://doi.org/10.1063/1.2822891)

[^Broido2005]: [Broido, D. A., Ward, A., & Mingo, N. (2005). Lattice thermal conductivity of silicon from empirical interatomic potentials. Physical Review B, 72(1), 1–8.](http://doi.org/10.1103/PhysRevB.72.014308)
