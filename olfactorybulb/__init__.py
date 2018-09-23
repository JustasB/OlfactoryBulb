from database import Layer, GlomeruliLocation, Parameter, CellLocation

class OlfactoryBulb:
    def __init__(self):
        # Load parameters
        self.params = list(Parameter.select())

        # Load layers
        self.layers = list(Layers.select())

        # Load glom positions
        self.glomeruli = list(GlomeruliLocation.select())

        # Load cell positions
        self.glomeruli = list(CellLocation.select())


    def build(self):
        # Synapses
        # Inputs
        # Recordings
        # LFP
        pass

    def run(self):
        pass

if __name__ == "__main__":
    model = OlfactoryBulb()
    model.run()