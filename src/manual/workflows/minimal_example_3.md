title: Minimal example III
author: Olle Hellman

### Fast stochastic calculations

Initially, the TDEP method was based on sampling the correct region of phase space via ab initio molecular dynamics. Since then I have also been using stochastic methods of sampling phase space, essentially self-consistent phonons.

@note I do have a formal derivation that minimizing the difference in forces and minimizing the error in free energy is equivalent, and that the former is just a very efficient way of realizing the latter. Should probably publish it in a paper before I publish it here though.

The algorithm is simple:

* Create a set of supercells that sample phase space from a set of phonons
* Determine the effective harmonic Hamiltonian
* With new phonons, create new supercells and repeat until convergence

There are some subtleties here though. The first issue is where to start: since it is a self-consistent procedure the speed will depend on how close to the final answer the initial guess is.
The second issue is about efficiency: how many configurations do you choose in each iteration, and is there a scheme to speed up convergence? If done right, stochastic calculations can be extremely efficient. If done wrong, not so much.

In fact, for weakly anharmonic systems stochastic sampling can be faster than finite difference calculations or DFPT. The cost is linear in the number of Wyckoff positions, approximately.

## Geometric series

I will start by describing the simplest possible scheme. It might not always be the fastest, but empirically I have found it to be very robust. For starters, prepare an empty folder that contain the unit cell. My suggestion is to use bcc Li as an example, mostly because DFT calculations of Li are fast, but also because it is anharmonic enough to not be completely trivial.



(I suggest practicing on a simple metal, bcc Li is a great example) Build a supercell of around 100 atoms



##

@todo Make example with classical potential?
