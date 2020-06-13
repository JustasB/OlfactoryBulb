from prev_ob_models.utils import RunInClassDirectory, IsolatedCell

class BaseMC(IsolatedCell):
    def __init__(self):
        with RunInClassDirectory(BaseMC):
            from neuron import h,gui

            self.h = h
            self.cell = self.get_cell()

            h.celsius=35

            self.soma = self.cell.soma

            h.cvode_active(1)

class MC1(BaseMC):
    def get_cell(self):
        self.h.load_file("Cell1.hoc")
        return self.h.Cell1("1","1","1")

class MC2(BaseMC):
    def get_cell(self):
        self.h.load_file("Cell2.hoc")
        return self.h.Cell2("2","2","2")

class MC3(BaseMC):
    def get_cell(self):
        self.h.load_file("Cell3.hoc")
        return self.h.Cell3("3","3","3")

class MC4(BaseMC):
    def get_cell(self):
        self.h.load_file("Cell4.hoc")
        return self.h.Cell4("4","4","4")

class MC5(BaseMC):
    def get_cell(self):
        self.h.load_file("Cell5.hoc")
        return self.h.Cell5("5","5","5")

class MC6(BaseMC):
    def get_cell(self):
        self.h.load_file("Cell6.hoc")
        return self.h.Cell6("6","6","6")



