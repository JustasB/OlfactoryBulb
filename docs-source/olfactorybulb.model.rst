olfactorybulb.model
=========================================

The main class to build and run the network model. The class constructor builds the model using one of the
specified `parameter classes <olfactorybulb.paramsets.html>`__. If `autorun==True`, the model will be simulated
after building. Otherwise, the `run(tstop) <#olfactorybulb.model.OlfactoryBulb.run>`__ is used to run the simulation.

MPI/multi-core simulations can be performed using `initslice.py` found in the repo root. For example, the following
command will build and run the network model using the
`GammaSignature parameter set <olfactorybulb.paramsets.html#olfactorybulb.paramsets.case_studies.GammaSignature>`_
using 16 cores::

    # From repo root, run this shell command:
    mpiexec -np 16 python initslice.py -paramset GammaSignature -mpi

Simulation results will be stored under `[repo]/olfactorybulb/results/GammaSignature`. The last folder in the path uses
the name of the parameter class.

See the
`LFP Wavelet Analysis.ipynb <https://github.com/JustasB/OlfactoryBulb/blob/master/notebooks/LFP%20Wavelet%20Analysis.ipynb>`_
notebook for examples of how the results are analyzed.

.. automodule:: olfactorybulb.model
    :members:
    :show-inheritance: