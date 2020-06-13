# Run with: mpiexec -n 2 python testmpi.py -mpi

import sys

if '-mpi' in sys.argv:
    from mpi4py import MPI
    print('MPI: ON')
else:
    print('MPI: OFF')

from time import sleep
from neuron import h, gui, load_mechanisms
import os, sys

load_mechanisms('prev_ob_models/Birgiolas2020/Mechanisms')

pc = h.ParallelContext()
id = pc.id()
nhost = pc.nhost()

print('On:', id, nhost)

gg = 0.01

if id == 0:
    print('Creating 1')
    soma1 = h.Section(name='soma1')
    soma1.insert('pas')
    soma1.insert('hh')
    soma1.L = soma1.diam = 10
    soma1.nseg = 3

    gap1 = h.GapJunction(0.5, sec=soma1)
    gap1.g = gg



    ic = h.IClamp(soma1(0.5))
    ic.amp = 2
    ic.dur = 0.1
    ic.delay = 20

    # pc.set_gid2node(0, 0)
    pc.source_var(soma1(0.5)._ref_v, 10)
    pc.target_var(gap1._ref_v_other, 20)

    vv1 = h.Vector()
    vv1.record(soma1(0.5)._ref_v)

if id == 1:
    print('Creating 2')
    soma2 = h.Section(name='soma2')
    soma2.insert('pas')
    soma2.insert('hh')

    gap2 = h.GapJunction(0.5, sec=soma2)
    gap2.g = gg

    # pc.set_gid2node(1, 1)
    pc.source_var(soma2(0.5)._ref_v, 20)
    pc.target_var(gap2._ref_v_other, 10)

    vv2 = h.Vector()
    vv2.record(soma2(0.5)._ref_v)

if id in [0,1]:
    tv = h.Vector()
    tv.record(h._ref_t)


pc.setup_transfer()
pc.set_maxstep(1)
h.stdinit()
h.dt = 0.1
h.tstop = 50
pc.psolve(h.tstop)


if id == 0:
    from matplotlib import pyplot as plt
    plt.plot(tv.to_python(), vv1.to_python())
    plt.title(id)
    plt.show()

if id == 1:
    from matplotlib import pyplot as plt
    plt.plot(tv.to_python(), vv2.to_python())
    plt.title(id)
    plt.show()
