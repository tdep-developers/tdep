
### Short description

Calculate the three-phonon scattering surface.

### Command line options:




Optional switches:

* `--qpoint value#1 value#2 value#3`  
    default value 0 0 0  
    Specify the q-point for which to calculate the scattering surface.

* `--highsymmetrypoint value`  
    default value none  
    Same as above, but you can specify the label of a high-symmetry point instead, e.g. "X" or "L".

* `--modespec value#1 [value#2...]`  
    default value 1 2 3 0 1 100  
    Specify which surfaces to look at. 6 integers per surface: band1 band2 band3 plus/minus color opacity. Modes go from 1 to number of bands. plus/minus is 0 for plus events 1 for minus. Color from 1-3, opacity from 0-100.

* `--povray`  
    default value .false.  
    Write POV-Ray script for rendering.

* `--intensities`  
    default value .false.  
    Calculate intensities (matrix elements) as well.

* `--pv_theta_phi value#1 value#2`  
    default value 20 30  
    Povray specification, theta and phi angles (in degrees) deteriming camera location.

* `--pv_quality value`  
    default value 9  
    Povray rendering quality, 1-9 where 1 is worst.

* `--nband value`  
    default value -1  
    Specify first band index

* `--qpoint_grid value#1 value#2 value#3`, `-qg value#1 value#2 value#3`  
    default value 26 26 26  
    Density of q-point mesh for Brillouin zone integrations.

* `--readqmesh`  
    default value .false.  
    Read the q-point mesh from file. To generate a q-mesh file, see the genkpoints utility.

* `--help`, `-h`  
    Print this help message

* `--version`, `-v`  
    Print version
### Examples

`phasespace_surface --highsymmetrypoint X` 

`phasespace_surface --qpoint 0.1 0.2 0.3` 
