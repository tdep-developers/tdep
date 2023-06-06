author: Olle Hellman
display: none
graph: none
propname: generate structure
propnamelink: <a href="../program/generate_structure.html">generate structure</a>
{!man/generate_structure.md!}

### Longer summary

This code builds a supercell with user-specified dimensions from the unit cell given in `infile.ucposcar`. All positions are in fractional coordinates.

The code handles diagonal and non-diagonal cells. The diagonal cells are simple repetitions of the unit cell, making it $N_1 \times N_2 \times N_3$ larger. The non-diagonal cells can be useful when the unit cell has an awkward shape. The supercell will have the lattice vectors

$$
\begin{pmatrix}
\mathbf{A}_1 \\
\mathbf{A}_2 \\
\mathbf{A}_3 \\
\end{pmatrix}
= \mathbf{M}
\begin{pmatrix}
\mathbf{a}_1 \\
\mathbf{a}_2 \\
\mathbf{a}_3
\end{pmatrix}
$$

where $\det( \mathbf{M} )$ is a positive integer. Some useful transformations are, for example, with `-nd -n1 n1 n1 n2 -n2 n2 n3 n3 -n3`, from a primitive fcc lattice to the conventional cubic

$$
\begin{pmatrix}
n_1 a & 0 & 0 \\
0 & n_2 a & 0 \\
0 & 0 & n_3 a
\end{pmatrix}
=\begin{pmatrix}
-n_1 & n_1 & n_1 \\
n_2 & -n_2 & n_2 \\
n_3 & n_3 & -n_3
\end{pmatrix}
\begin{pmatrix}
0 & a/2 & a/2 \\
a/2 & 0 & a/2 \\
a/2 & a/2 & 0
\end{pmatrix}
$$

Similarly, from a primitive bcc to conventional cubic:

$$
\begin{pmatrix}
n_1 a & 0 & 0 \\
0 & n_2 a & 0 \\
0 & 0 & n_3 a
\end{pmatrix}
=\begin{pmatrix}
0 & n_1 & n_1 \\
n_2 & 0 & n_2 \\
n_3 & n_3 & 0
\end{pmatrix}
\begin{pmatrix}
-a/2 & a/2 & a/2 \\
a/2 & -a/2 & a/2 \\
a/2 & a/2 & -a/2
\end{pmatrix}
$$

There are some non-obvious ones as well, such as converting a hexagonal lattice to an orthorhombic:

$$
\begin{pmatrix}
n_1 a & 0 & 0 \\
0 & \sqrt{3} n_2 a  & 0 \\
0 & 0 & n_3 c
\end{pmatrix}
=\begin{pmatrix}
n_1 & n_1 & 0 \\
-n_2 & n_2 & 0 \\
0 & 0 & n_3
\end{pmatrix}
\begin{pmatrix}
a/2 & -a\sqrt{3}/2 & 0 \\
a/2 & a \sqrt{3}/2 & 0 \\
0 & 0 & c
\end{pmatrix}
$$

or rhombohedral (lattice parameter $a$, angle $\alpha$) to hexagonal:

$$
\begin{equation*}
\begin{split}
\begin{pmatrix}
0 & n_1 a \sqrt{2-2 \cos\alpha} & 0  \\
-n_2 a\sqrt{ \frac{2-2 \cos\alpha}{3} } & -n_2 a \sqrt{2-2 \cos\alpha}  & 0 \\
0 & 0 & n_3 a \sqrt{3+6\cos\alpha}
\end{pmatrix}
= \\ \\ \begin{pmatrix}
n_1 & n_1 & 0 \\
-n_2 & n_2 & 0 \\
0 & 0 & n_3
\end{pmatrix}
\begin{pmatrix}
a\sqrt{2-2\cos\alpha} & 0 & a\sqrt{\frac{1+2\cos\alpha}{3}} \\
-a\sqrt{\frac{1-\cos\alpha}{2}} & a\sqrt{3\frac{1-\cos\alpha}{2}} & a\sqrt{\frac{1+2\cos\alpha}{3}} \\
-a\sqrt{\frac{1-\cos\alpha}{2}} & -a\sqrt{3\frac{1-\cos\alpha}{2}} & a\sqrt{\frac{1+2\cos\alpha}{3}}
\end{pmatrix}
\end{split}
\end{equation*}
$$

There are too many ways to build supercells to list here. The examples are suggestions so that you can pick a supercell that is as cubic as possible. Converging results with respect to supercell size goes approximately as the size of the largest sphere you can fit in the cell.

#### Automatic generation

The option `-na` will try to determine a non-diagonal supercell that is as cubic as possible, with the number of atoms close to the specified value. It does it by generating a massive amount of possible supercell transformation matrices, and comparing the ratio of the cell volume and the volume of the inscribed sphere. By optimizing with respect to this ratio, a reasonably cubic cell should be the result. Very useful when you have a unitcell with a tricky shape and it's non-trivial to figure out a meaningful transformation matrix.

To clarify, a cell defined by vectors $\mathbf{a}$, $\mathbf{b}$ and $\mathbf{c}$ will have the volume

$$
V=\left|\det
\begin{pmatrix}
  a_x & a_y & a_z \\
  b_x & b_y & b_z \\
  c_x & c_y & c_z
\end{pmatrix}\right|
$$

In the same cell, the radius of the largest sphere that can be inscribed is given by

$$
2r=\min\left\{
\left|\frac{ \mathbf{b} \times \mathbf{c} }{\left| \mathbf{b} \times \mathbf{c} \right|}\cdot \mathbf{a}\right| \,,\quad
\left|\frac{ \mathbf{a} \times \mathbf{b} }{\left| \mathbf{a} \times \mathbf{b} \right|}\cdot \mathbf{c}\right| \,,\quad
\left|\frac{ \mathbf{a} \times \mathbf{c} }{\left| \mathbf{a} \times \mathbf{c} \right|}\cdot \mathbf{b}\right|
\right\}  
$$

I define the ratio

$$
f^3 = \frac{r^3 4 \pi}{3V} \le 1
$$

As the function to be maximized (it is the ratio of the filling ratio of a cell and the filling ratio of a cube with the same volume, a value of 1 indicates a perfect cube). The algorithm works by searching the space of non-diagonal supercell matrices that produce approximately the desired number of atoms and simultaneously maximize the filling ratio. This algorithm produces rather non-intuitive cells. Consider Bi<sub>2</sub>Te<sub>3</sub> defined in the primitive rhombohedral cell:

```
Bi2Te3
  10.314046162
 0.243360208934  0.000000000000  0.969935981757
-0.121680104467  0.210756123207  0.969935981757
-0.121680104467 -0.210756123207  0.969935981757
Bi Te
2 3
direct coordinates
 0.599898812406  0.599898812406  0.599898812406 site: 1
 0.400101187594  0.400101187594  0.400101187594 site: 2
 0.791308614612  0.791308614612  0.791308614612 site: 3
 0.208691385388  0.208691385388  0.208691385388 site: 4
 0.000000000000  0.000000000000  0.000000000000 site: 5
```

Using this algorithm I found the supercell matrix

$$
M=\begin{pmatrix}
3 & 2 & -4 \\
4 & -3 & -2 \\
-2 & 4 & -3
\end{pmatrix}
$$

that produce the lattice vectors

```
Bi2Te3 upercell
      10.314046161996
    0.97344083573636     1.26453673924247     0.96993598175736
    1.58184135807159    -0.21075612320708    -0.96993598175736
   -0.60840052233523     1.47529286244955    -0.96993598175736
 Bi Te
 86 129
```
This cell has a filling ratio corresponding to 98.4% of the ideal cube.

### Input files

* [infile.ucposcar](../page/files.html#infile.ucposcar)

### Output files

* [outfile.ssposcar](../page/files.html#infile.ssposcar)
