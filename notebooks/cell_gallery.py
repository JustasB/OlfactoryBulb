# Start this script with python cell_gallery.py
# Then open Blender and import all cells
# Creates all 15 cell models in NEURON, each showing a spike
# Open Blender with BlenderNEURON and import groups with activity
# Reposition cells in Blender as needed

from blenderneuron import neuronstart
from prev_ob_models.Birgiolas2020.isolated_cells import *
from neuron import h, gui

cells = [MC1(),MC2(), MC3(), MC4(), MC5(), 
         TC1(), TC2(), TC3(), TC4(), TC5(), 
         GC1(), GC2(), GC3(), GC4(), GC5()]

ics = [h.IClamp(0.5, sec=cell.soma) for cell in cells]

h.tstop = 20
h.celsius = 35

# MCs spike
for ic in ics[0:5]: ic.delay=5; ic.amp = 1; ic.dur=10;

# TCs
for ic in ics[5:10]: ic.delay=5; ic.amp = 0.3; ic.dur=10;

# GCs
for ic in ics[10:15]: ic.delay=5; ic.amp = 0.1; ic.dur=10;

