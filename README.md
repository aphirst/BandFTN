BandFTN
=======

Electronic bandstructure calculations using the empirical pseudopotential method.

Implemented in Fortran 2008, using some object-oriented approaches, to ease development and extension. Effort has been made to ensure that the code is either self-documenting, or sufficiently well-commented so that any interested physicist can understand the mechanics.

Licensed under GPLv3, to ensure that everyone is free to learn from the code, and that the code itself remains free-as-in-freedom.


How to Use
----------

Currently there is only support for GNU/Linux, and also required are an up-to-date version of gfortran, a LAPACK implementation (I personally recommend the multithreaded [OpenBLAS](http://www.openblas.net/)), and gnuplot.

Edit the Makefile to set sane build and object directories (defaults to inside ```/tmp```), run ```make```, and if compilation is successful the program can simply be run with ```./BandFTN```. Note that for the time being, the parameters of the calculation are hard-coded; rather than document their usage I shall instead prioritise implementing the ability to specify them at runtime.

The generated ```.gnu``` file is a gnuplot command file; running gnuplot and ```load```ing it should display the desired bandstructure graph, the raw data for which is contained in an accompanying generated ```.dat``` file.


Current Features
----------------

* Generates an EPM bandstructure for FCC crystals along lines of high symmetry
* Anti-symmetric structure/form-factors are also taken into account for multi-element crystals, e.g. GaAs
* Zeros energy axis at valence band maximum, in accordance with convention
* Irreducible Brillouin Zone k-point mesh generation, and density of states calculations
* Outputs a data file and a Gnuplot command file


Planned Features
----------------

* Determination & evaluation of band gaps
* *Plotting* of density of states data
* Material/crystal property input/job files
* Additional cubic lattice types
* Charge density calculations
* Support for other operating systems and compilers

Informal References
-------------------

* Two brief project reports posted online (one by [Aaron J. Danner](http://www.ece.nus.edu.sg/stfpage/eleadj/pseudopotential.htm), another by [Kevin D. Welsher](http://large.stanford.edu/courses/2007/ap272/welsher1/)); which clarified the calculation and usage of the asymmetric form-factors
* A paper by [W. Setyawan and S. Curtarolo](http://arxiv.org/abs/1004.2974); providing an excellent, consistent resource for the high-symmetry points, lines, and irreducible Brillouin zone wedges
