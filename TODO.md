To-Do Notes for BandFTN
=======================

1. Re-organise file structure, one class per file if possible

2. Implement symmetry-operation and Monkhorst-Pack generation as a module
   * Use matrix method rather than "fold shifters", which are nigh-on-impossible to understand or maintain
   * Needs to be generalisable to Cubic symmetries with Oh < 48

3. Determine sane data structure for Energies and "weights", which can be easily-parsed

4. Verify equations (e.g. [Prof. Clark's PhD](http://cmt.dur.ac.uk/sjc/thesis_prt/node104.html)) for Density of States, and total energy

5. Decide whether to unify DoS output data with bandstructure data

6. Determine how to separate "bandstructure" and "DOS" parts of gnuplot-file generation

7. Book up on ```plplot```
