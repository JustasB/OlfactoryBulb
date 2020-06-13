"""
This is a helper file for running multi-core/MPI simulations. For example, this shell command:

mpiexec -np 16 python initslice.py -paramset GammaSignature -mpi

Will run the network model using 16 cores using the olfactorybulb.paramsets.GammaSignature parameters
"""

import sys
if '-mpi' in sys.argv:
    from mpi4py import MPI

from olfactorybulb.model import OlfactoryBulb as OB

if '-paramset' in sys.argv:
    paramset = sys.argv[sys.argv.index("-paramset")+1]
    ob = OB(paramset)

else:
    ob = OB()
