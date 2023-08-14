
### Short description

The main algorithm of the TDEP method. Starting with a symmetry analysis, this code finds the irreducible representation of interatomic forceconstants and extracts them from position and force data.

### Command line options:




Optional switches:

* `--secondorder_cutoff value`, `-rc2 value`  
    default value 5.0  
    Cutoff for the second order force constants

* `--thirdorder_cutoff value`, `-rc3 value`  
    default value -1
    mutually exclude "--thirdorder_njump"  
    Cutoff for the third order force constants

* `--fourthorder_cutoff value`, `-rc4 value`  
    default value -1
    mutually exclude "--fourthorder_njump"  
    Cutoff for the fourth order force constants

* `--polar`  
    default value .false.  
    Add dipole-dipole corrections for polar materials.

* `--stride value`, `-s value`  
    default value 1  
    Use every N configuration instead of all. Useful for long MD simulations with linearly dependent configurations.

* `--firstorder`  
    default value .false.  
    Include the first order force constants. These can be used to find the finite temperature equilibrium structure.

* `--readforcemap`  
    default value .false.  
    Read `infile.forcemap.hdf5` from file instead of calculating all symmetry relations. Useful for sets of calculations with the same structure.

* `--readirreducible`  
    default value .false.  
    Read the irreducible forceconstants from `infile.irrifc_*` instead of solving for them. This option requires an `infile.forcemap.hdf5`, as above.

* `--potential_energy_differences`, `-U0`  
    default value .false.  
    Calculate the difference in potential energy from the simulation and the forceconstants to determine U0. As referenced in the thermodynamics section of [phonon dispersion relations](phonon_dispersion_relations.html) this is the renormalized baseline for the TDEP free energy: $$U_0= \left\langle U^{\textrm{BO}}(t)-\frac{1}{2} \sum_{ij}\sum_{\alpha\beta} \Phi_{ij}^{\alpha\beta} \mathbf{u}^{\alpha}_i(t) \mathbf{u}^{\beta}_j(t) \right\rangle$$ This number should be added to the appropriate phonon free energy.

* `--printforcemap`  
    default value .false.  
    Print `outfile.forcemap.hdf5` for reuse.

* `--temperature value`  
    default value -1  
    Temperature for self-consistent solver.

* `--help`, `-h`  
    Print this help message

* `--version`, `-v`  
    Print version
### Examples

`extract_forceconstants -rc2 5.1` 

`extract_forceconstants -rc2 4.5 -rc3 3.21` 
