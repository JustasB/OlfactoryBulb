from prev_ob_models.utils import RunInClassDirectory

class MC():
    def __init__(self):
        with RunInClassDirectory(MC):
            from neuron import h, gui

            h.xopen("mathslib.hoc")
            h.xopen("tabchannels.hoc")
            h.xopen("parameters_fig1.hoc")
            h.load_file('mitral.tem')

            self.cell = h.Mit(100)  # 100 is used in bulb.hoc

            h.define_shape()

            h.celsius = 6.3

            h.cvode_active(1)

            self.h = h
            self.soma = h.Mit[0].soma
