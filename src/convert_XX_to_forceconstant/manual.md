author: Olle Hellman
display: none
graph: none
propname: convert DFT dynamical matrices to TDEP
propnamelink: <a href="../program/convert_XX_to_forceconstant.html">convert DFT dynamical matrices to TDEP</a>
{!man/convert_XX_to_forceconstant.md!}

### Longer summary

This is a utility to convert DFT or DFPT inputs with dynamical matrices (e.g. Abinit DDB files) to the TDEP forceconstant format. Future adaptation for QE or phonopy inputs is foreseen. The main usage will be to seed optimal canonical configurations without having to make an MLIP or do some fake forceconstant model.

A short disclaimer should be that it is not trivial to figure out what information is written to the abinit DDB files by parsing the DDB file alone. The conversion script has been tested on hundreds of Abinit calculations, but your mileage might vary.

The script will try to figure out the q-mesh used in the DFPT calculations automatically, but in case you used a shifted mesh, you might want to specify it manually, since the information about the q-mesh is not printed as such in the DDB file.

@todo Add io
