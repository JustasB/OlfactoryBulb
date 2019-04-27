#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This program creates a Migliore & Shepherd 2008 gran cell model along with tables to pull data.
# Only the two biggest compartments are modelled.

import os
import sys
import math

sys.path.extend(["..","../channels","../neuroml","../simulations"])

from moose_utils import *

from load_channels import *
from MorphML_reader import *
from simset_inhibition import *

SETTLETIME = 50e-3 # s
RUNTIME = REALRUNTIME + SETTLETIME

from pylab import * # part of matplotlib that depends on numpy but not scipy

class granTest:

    def __init__(self):
        load_channels()
        MML = MorphML()
        gran_dict = MML.readMorphMLFromFile('../cells/granule_granadityaMS2007_neuroML_L1_L2_L3.xml',{})
        self.context = moose.PyMooseBase.getContext()
        granCellId = self.context.deepCopy(self.context.pathToId('/library/granule'),self.context.pathToId('/'),"granule")
        self.granCell = moose.Cell(granCellId)
        # note that soma is not just /granule/soma, it'll also have its segid appended - thus the general purpose function below:
        self.granSoma = moose.Compartment(get_matching_children(self.granCell, ['Soma','soma'])[0]) # take the first [0] available soma!!
        self.somaVm = setupTable('soma_vm', self.granSoma,'Vm')

if __name__ == "__main__":
    gran = granTest()
    iclamp = setup_iclamp(gran.granSoma, 'gran_iclamp', SETTLETIME, 100e-3, 500e-12) # 500pA for 100ms to get a single AP
    resetSim(gran.context, SIMDT, PLOTDT)
    gran.context.step(RUNTIME)
    timevec = arange(0.0,RUNTIME+1e-10,SIMDT)
    plot(timevec, gran.somaVm,'r,-', label='gran soma Vm (V)')
    legend()
    xlabel('time (s)')
    show()

