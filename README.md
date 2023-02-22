# RAMSES_Dipolar_Dark_Matter
My patched version of the N-body code RAMSES to simulate Dipolar Dark Matter (DDM).
RAMSES' documentation can be found at https://bitbucket.org/rteyssie/ramses/src/master/. Most changes to the main version are highlighted by a comment ```dipolar``` in the code.

The main idea of this project is to customize RAMSES in order to take into account a special type of dark matter which gravitates differently. Instead of the standard Poisson equation <img src="https://latex.codecogs.com/svg.latex?\Delta&space;\phi&space;=4&space;\pi&space;G&space;\rho" title="\Delta \phi =4 \pi G \rho" />, dipolar dark matter features a new term that physically can be interpreted as a gravitational polarization <img src="https://latex.codecogs.com/svg.image?\vec{\Pi}" titre="pi" />:

<img src="https://latex.codecogs.com/svg.image?\Delta&space;\phi&space;=4&space;\pi&space;G&space;\rho&space;-&space;\vec{\nabla}&space;\cdot&space;\vec{\Pi}" titre="Modified Poisson" />

<img src="https://latex.codecogs.com/svg.image?\vec{\Pi}" titre="pi" /> is a gravity quantity that lives on the RAMSES grid and has its own dynamics that has also been implemented in this patch. DDM also introduce a new *internal* force that affects only this very special type of dark matter but not the baryons. To implement it, dark matter is differentiated from the baryons with the tag of RAMSES. (Dipolar) Dark matter carries the tag 0, while the baryons have the tag 1.  

At the moment, the code has only been tested with isolated simulations of spherical dwarf galaxies. More details can be found in https://arxiv.org/abs/2209.07831

In order to run the code, one needs to provide initial conditions for this specific dynamics, we have done it for a spherical sphere of dark matter (see [here](https://github.com/cspotz/RAMSES_Bi-Poisson/blob/main/postprocess/postprocess.md) for an implementation in Mathematica), and a [King profile](https://github.com/GFThomas/MOND).
Once the initial conditions are provided by the user, one can run RAMSES in the standard way by running commands of the type ``bin/ramses3d namelist/ramses.nml``
