#!/usr/bin/env python

# This program creates a Davison mitral cell model along with tables to pull data
import sys
import math
# The PYTHONPATH should contain the location of moose.py and _moose.so
# files.  Putting ".." with the assumption that moose.py and _moose.so
# has been generated in ${MOOSE_SOURCE_DIRECTORY}/pymoose/ (as default
# pymoose build does) and this file is located in
# ${MOOSE_SOURCE_DIRECTORY}/pymoose/examples
# sys.path.append('../..')
try:
    from neuron import *
except ImportError:
    print "ERROR: Could not import neuron. Please add the directory containing neuron.py in your PYTHONPATH"
    import sys
    sys.exit(1)

from pylab import *

h.load_file("stdlib.hoc")
h.load_file("stdrun.hoc") # for h.stop ,etc.

h.xopen("/home/aditya/nrn-7.1/share/nrn/lib/hoc/celbild/celtopol.hoc")
h.xopen("/home/aditya/nrn-7.1/share/nrn/lib/hoc/celbild/celgeom.hoc")
h.xopen("/home/aditya/nrn-7.1/share/nrn/lib/hoc/celbild/inhomofn.hoc")
h.xopen("/home/aditya/nrn-7.1/share/nrn/lib/hoc/celbild/psubset.hoc")
h.xopen("/home/aditya/nrn-7.1/share/nrn/lib/hoc/celbild/celset.hoc")
h.xopen("/home/aditya/nrn-7.1/share/nrn/lib/hoc/celbild/celmemb.hoc")
h.xopen("/home/aditya/nrn-7.1/share/nrn/lib/hoc/celbild/celmang.hoc")
h.xopen("/home/aditya/nrn-7.1/share/nrn/lib/hoc/celbild/celbild1.hoc")

# USING PHYSIOLOGICAL UNITS AS NEURON USES THEM
TSTOP = 1600.0 # ms # 500ms delay, 600ms current injection, 500ms decay time

class PGTest:

    def __init__(self):        
        print "NeuroML import of Cell..."

        cellb = h.CellBuild() # Cell Builder # This requires all the above xopen loading of files!

        neuromlFile = "../PG_aditya2010_neuroML_L1_L2.xml"
        #neuromlFile = "Level2.xml"

        cellb.manage.neuroml(neuromlFile)

        cellb.continuous = 1
        cellb.cexport() # export the cell from Cell Builder to neuron top level

        print "Loaded NeuroML from: ",  neuromlFile
        self.PGsoma = h.soma
        ##### NOTE: soma is mapped from a cable, hence it is a section, not a segment! Hence soma.pas does not work.
        ##### soma(0.5) refers to a segment, hence soma(0.5).pas works.
        ##### Can also do for seg in soma.allseg(), but I typically ensure we have only one segment per section.
        
        self.stim = h.IClamp(0.5, sec=self.PGsoma) # current clamp in soma
        
        self.somavecV = h.Vector()
        self.vect = h.Vector()
        
    def run(self):
        self.stim.amp = 0.01 # nA i.e 10pA
        self.stim.delay = 500 #ms
        self.stim.dur = 600 #ms
        h.tstop = TSTOP
        h.dt = 0.001 # ms
        self.somavecV.record(self.PGsoma(0.5)._ref_v)
        self.vect.record(h._ref_t)
        h.run()
                
if __name__ == "__main__":
    pg = PGTest()
    print "soma diameter = ",pg.PGsoma.diam," microns."
    print "soma length = ",pg.PGsoma.L," microns."
    print "Specific axial resistance RA = ",pg.PGsoma.Ra," Ohm-cm."
    soma_crosssectionA = 3.14159*(pg.PGsoma.diam/2.0)**2 * 1e-8 # cm^2 from micron^2
    soma_surfaceA = 3.14159*pg.PGsoma.diam*pg.PGsoma.L * 1e-8 # cm^2 from micron^2
    print "soma Rm = ",(1.0/pg.PGsoma(0.5).pas.g)/soma_surfaceA," Ohms."
    print "soma Cm = ",pg.PGsoma.cm*soma_surfaceA*1e-6," Farads." # Farads from microFarads
    print "soma Ra = ",pg.PGsoma.Ra * pg.PGsoma.L*1e-4/soma_crosssectionA," Ohms."
    print "soma Na gmax = ",pg.PGsoma(0.5).Na_rat_ms.gmax*soma_surfaceA," Siemens."
    print "soma K gmax = ",pg.PGsoma(0.5).KDR_ms.gmax*soma_surfaceA," Siemens."
    pg.run()
    plot(pg.somavecV,',-r')
    show()
