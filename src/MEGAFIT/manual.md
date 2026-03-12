author: Olle Hellman
display: none
graph: none
propnamelink: <a href="../program/MEGAFIT.html">MEGAFIT</a>
propname: MEGAFIT
{!man/MEGAFIT.md!}

### Longer summary

This code calculates interatomic force constants as a function of intensive quantities. In short, it gives you the tools to create a continuous Hamiltonian:

$$
H = U(\boldsymbol{x})
+\frac{1}{2}\sum_{ij} \Phi_{ij}^{\alpha\beta}(\boldsymbol{x})
u_{i}^{\alpha}(\boldsymbol{x})
u_{j}^{\beta}(\boldsymbol{x})
+\frac{1}{3}\sum_{ij} \Phi_{ijk}^{\alpha\beta\gamma}(\boldsymbol{x})
u_{i}^{\alpha}(\boldsymbol{x})
u_{j}^{\beta}(\boldsymbol{x})
u_{k}^{\gamma}(\boldsymbol{x})
+\ldots
$$

Where $x$ could be volume, temperature $c/a$-ratio or anything really. It means you can go from simulations at a set of discrete points in $V-T$ space

_IMAGE_

to a continuous function like this:

_IMAGE_

This is very handy for thermodynamics, in the sense that any derivative in $x$ is well defined, so that quantitities like pressure, entropy and heat capacity can be evaluated to high accuracy.

#### Interpolating force constants

Ever since my first project using TDEP, I have used this idea.

#### Special variables

It is generally a good idea not to treat all possible $x$ equally. First up we check volume

###### Volume

If we inspect the well-known Birch-Murnaghan equation of state, it can be written as

$$
\begin{equation}\label{eq:birchmur_energy}
E = A_0 + A_1 \nu + A_2 \nu^2 + A_3 \nu^3
\end{equation}
$$

Where

$$
\begin{equation}
\nu = \left(\frac{V_0}{V}\right)^{2/3}
\end{equation}
$$

This is generally a rather good functional form for interpolating energy-volume curves. The second derivative with respect to volume of


@todo Make Tscale automatic

@todo Add max order per dimension

@todo Add contraction to 2D EOS

@todo Add convex hull test when evaluating

@todo Sanity check of infile.evalpoints early in the code

@todo Figure out a linspace that always includes zero, for pressure

@todo Get dispersions vs T printout

@todo Get the structure interpolation working ok. I want it do be written in terms of a symmetry-reduced strain, and a symmetry-reduced reference position and shift.
