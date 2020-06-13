# Plots current clamp responses of each cell type
# Uncomment the desired cell type line below and then
# Start this script with python cell_current_responses.py


from prev_ob_models.Birgiolas2020.isolated_cells import *
from neuron import h, gui

#stim_set = ([MC1(),MC2(), MC3(), MC4(), MC5()],  [0.25, -0.1, 0])
#stim_set = ([TC1(), TC2(), TC3(), TC4(), TC5()], [0.25, -0.1, 0])
stim_set = ([GC1(), GC2(), GC3(), GC4(), GC5()], [0.06, -0.02, 0])

cells = stim_set[0]
amps = stim_set[1]

ics = [h.IClamp(0.5, sec=cell.soma) for cell in cells]


h.cvode_active(0)
delay = 200
dur   = 700
h.tstop = 1000
h.celsius = 35
h.steps_per_ms = 10
h.dt = 1.0 / h.steps_per_ms

# Plots
for cell in cells:
    h.newPlotV()
    h.Graph[-1].erase_all()
    h.Graph[-1].addvar(cell.__class__.__name__ + '[0].soma.v(0.5)')
    h.Graph[-1].exec_menu("Keep Lines")


for amp in amps:
    for ic in ics: 
        ic.delay=delay; ic.amp = amp; ic.dur=dur;

    h.run()








