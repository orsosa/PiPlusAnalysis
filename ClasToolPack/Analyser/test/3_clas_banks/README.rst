CLASTool and Analyser example programs
======================================

This source code is used to analyse the ``clas_42011_01_1.pass2.root`` file,
or any similar file.


The ``create_hist.cxx`` program reads the ROOT file and creates the *betta*
versus *momentum* histograms (TH2F and TProfile) for every particle, and saves
them in two files, ``particle_histograms.root`` and
``particle_profiles.root``.  The TProfile objects are fitted with the betta
function.  Now you can use the macros ``show_hist.C`` to see the TH2F
histograms, ``show_prof.C`` to see the TProfile histograms and ``show_fit.C``
to see the fitting of the profiles with the betta function (which is in the
``betta.C`` file).


An alternative way is provided by the ``write_tree.cxx`` program. It creates a
TTree with branches for *betta*, *moment* and *particle_id* values of every
particle. This tree is saved in the ``particle_data.root`` file.  Then you can
use the ``getFit.C`` macro to show the fitting of the particle given as
parameter to the macro.


The ``mass_distribution.cxx`` program creates the *mass* versus *momentum*
histograms (TH2F and TProfile) for every particle with positive charge and
``StatSC > 1``.  The histograms are saved in the ``particle_mass_dist.root``
file. You can use the macro ``show_mass_dist.C`` to show them.


The three programs listed above are portable, i.e., you can compile them using
the provided Makefile, or you can use the CINT interpreter and run them as
normal macros. But, you can not use ACLiC to compile them.
