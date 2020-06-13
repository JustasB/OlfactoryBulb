"""
This file is used to sequentially run the network model using different sets of parameters

The paramsets array should contain class names found in: [repo]/olfactorybulb/paramsets/*.py
"""

import os

paramsets = [
    "GammaSignature",
    "GammaSignature_NoInhibition",
    "GammaSignature_NoTCGJs",
    "GammaSignature_NoMCGJs",
    "GammaSignature_EqualTCMCInputs"
]


for i, params in enumerate(paramsets):
    print('Starting paramset: ' + params + ' (%s/%s)...' % (i+1, len(paramsets)))
    os.system('mpiexec -np 16 python initslice.py -paramset '+params+' -mpi')
