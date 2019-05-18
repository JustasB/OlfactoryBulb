from database import *

class OlfactoryBulb:

    def __init__(self):
        from netpyne import specs as netpyne_specs, sim as netpyne_sim

        # Load parameters

        # import pydevd
        # pydevd.settrace('192.168.0.34', port=4200, stdoutToServer=True, stderrToServer=True)

        self.properties = {param.name: eval(param.value) for param in Property.select()}

        # Load layers
        self.layers = {layer.id: layer for layer in Layer.select()}

        # Load layer vertices
        self.verts = list(LayerVertex.select())

        # Load glom positions
        self.glomeruli = list(Glomerulus.select().where(self.within_modeled_region(Glomerulus)))

        # Load cell types
        self.cell_types = {type.id: type for type in CellType.select()}

        # Load cell positions
        self.cells = list(Cell.select().where(self.within_modeled_region(Cell)))

        self.lfp_electrodes = list(LfpElectrode.select())

    def build(self):
        self.netParams = netpyne_specs.NetParams()

        #import pydevd
        #pydevd.settrace('192.168.0.34', port=4200, stdoutToServer=True, stderrToServer=True)

        netpyne_cells = { 'MC': [], 'GC': [], 'TC': [], 'PGC': [] }

        for i, cell in enumerate(self.cells):
            netpyne_cell = {
                'cellLabel': cell.type.id + str(i).zfill(6),
                'x': cell.x,
                'y': cell.y,
                'z': cell.z

            }
            netpyne_cells[cell.type.id].append(netpyne_cell)


        self.netParams.popParams['MC'] =  {'cellType': "MC", 'cellModel': 'HH', 'cellsList': netpyne_cells['MC']}
        self.netParams.popParams['TC'] =  {'cellType': "TC", 'cellModel': 'HH', 'cellsList': netpyne_cells['TC']}
        self.netParams.popParams['GC'] =  {'cellType': "GC", 'cellModel': 'HH', 'cellsList': netpyne_cells['GC']}
        self.netParams.popParams['PGC'] = {'cellType': "PGC", 'cellModel': 'HH', 'cellsList': netpyne_cells['PGC']}

        self.netParams.cellParams['MC'] = {
            'conds': { 'cellType': 'MC' },
            'secs': {
                'soma': {
                    'geom': {
                        'diam': 18.8,
                        'L': 18.8,
                        'Ra': 123.0
                    },
                    'mechs': {
                        'hh': {
                            'gnabar': 0.12,
                            'gkbar': 0.036,
                            'gl': 0.003,
                            'el': -70.0
                        }
                    },
                }
            }
        }
        self.netParams.cellParams['TC'] = {
            'conds': { 'cellType': 'TC' },
            'secs': {
                'soma': {
                    'geom': {
                        'diam': 10,
                        'L': 10,
                        'Ra': 123.0
                    },
                    'mechs': {
                        'hh': {
                            'gnabar': 0.12,
                            'gkbar': 0.036,
                            'gl': 0.003,
                            'el': -70.0
                        }
                    },
                }
            }
        }

        self.netParams.cellParams['GC'] = {
            'conds': { 'cellType': 'GC' },
            'secs': {
                'soma': {
                    'geom': {
                        'diam': 5,
                        'L': 5,
                        'Ra': 123.0
                    },
                    'mechs': {
                        'hh': {
                            'gnabar': 0.12,
                            'gkbar': 0.036,
                            'gl': 0.003,
                            'el': -70.0
                        }
                    },
                }
            }
        }

        self.netParams.cellParams['PGC'] = {
            'conds': { 'cellType': 'PGC' },
            'secs': {
                'soma': {
                    'geom': {
                        'diam': 5,
                        'L': 5,
                        'Ra': 123.0
                    },
                    'mechs': {
                        'hh': {
                            'gnabar': 0.12,
                            'gkbar': 0.036,
                            'gl': 0.003,
                            'el': -70.0
                        }
                    },
                }
            }
        }

        self.netParams.synMechParams['exc'] = {
            'mod': 'Exp2Syn',
            'tau1': 0.1,
            'tau2': 5.0,
            'e': 0,
        }

        self.netParams.synMechParams['inhib'] = {
            'mod': 'Exp2Syn',
            'tau1': 0.1,
            'tau2': 5.0,
            'e': -70,
        }

        self.netParams.stimSourceParams['bkg'] = {
            'type': 'NetStim',
            'rate': 10,
            'noise': 0.5
        }

        self.netParams.stimTargetParams['bkg->MC'] = {
            'source': 'bkg',
            'conds': { 'cellType': 'MC' },
            'weight': 0.01,
            'delay': 5,
            'synMech': 'exc'
        }

        self.netParams.stimTargetParams['bkg->TC'] = {
            'source': 'bkg',
            'conds': { 'cellType': 'TC' },
            'weight': 0.01,
            'delay': 5,
            'synMech': 'exc'
        }

        self.netParams.connParams['MC->GC'] = {
            'preConds': {'pop': 'MC'},
            'postConds': {'pop': 'GC'},
            'probability': 1.0,
            'weight': 0.01,
            'delay': 5,
            'sec': 'soma',
            'loc': 0.5,
            'synMech': 'exc'
        }

        self.netParams.connParams['GC->MC'] = {
            'preConds': {'pop': 'GC'},
            'postConds': {'pop': 'MC'},
            'probability': 1.0,
            'weight': 0.01,
            'delay': 1,
            'sec': 'soma',
            'loc': 0.5,
            'synMech': 'inhib'
        }

        self.netParams.connParams['TC->GC'] = {
            'preConds': {'pop': 'TC'},
            'postConds': {'pop': 'GC'},
            'probability': 1.0,
            'weight': 0.01,
            'delay': 5,
            'sec': 'soma',
            'loc': 0.5,
            'synMech': 'exc'
        }

        self.netParams.connParams['GC->TC'] = {
            'preConds': {'pop': 'GC'},
            'postConds': {'pop': 'TC'},
            'probability': 1.0,
            'weight': 0.01,
            'delay': 1,
            'sec': 'soma',
            'loc': 0.5,
            'synMech': 'inhib'
        }

        self.netParams.connParams['MC->PGC'] = {
            'preConds': {'pop': 'MC'},
            'postConds': {'pop': 'PGC'},
            'probability': 1.0,
            'weight': 0.01,
            'delay': 5,
            'sec': 'soma',
            'loc': 0.5,
            'synMech': 'exc'
        }

        self.netParams.connParams['PGC->MC'] = {
            'preConds': {'pop': 'PGC'},
            'postConds': {'pop': 'MC'},
            'probability': 1.0,
            'weight': 0.01,
            'delay': 1,
            'sec': 'soma',
            'loc': 0.5,
            'synMech': 'inhib'
        }

        self.netParams.connParams['TC->PGC'] = {
            'preConds': {'pop': 'TC'},
            'postConds': {'pop': 'PGC'},
            'probability': 1.0,
            'weight': 0.01,
            'delay': 5,
            'sec': 'soma',
            'loc': 0.5,
            'synMech': 'exc'
        }

        self.netParams.connParams['PGC->TC'] = {
            'preConds': {'pop': 'PGC'},
            'postConds': {'pop': 'TC'},
            'probability': 1.0,
            'weight': 0.01,
            'delay': 1,
            'sec': 'soma',
            'loc': 0.5,
            'synMech': 'inhib'
        }

    def setup_simulation(self):
        self.simConfig = netpyne_specs.SimConfig()

        self.simConfig.recordLFP = [[e.x, e.y, e.z] for e in self.lfp_electrodes]

        self.simConfig.duration = 1000
        self.simConfig.dt = 1 / 400.0
        self.simConfig.verbose = False
        self.simConfig.recordTraces = {
            'V_soma': {
                'sec': 'soma',
                'loc': 0.5,
                'var': 'v'
            }
        }
        self.simConfig.recordStep = 0.5
        self.simConfig.filename = 'model_output'
        self.simConfig.savePickle = False

        self.simConfig.analysis['plotRaster'] = True
        self.simConfig.analysis['plotTraces'] = {'include': [0, 1]}
        self.simConfig.analysis['plot2Dnet'] = True
        self.simConfig.analysis['plotLFP'] = True

    def run(self):
        netpyne_sim.createSimulateAnalyze(self.netParams, self.simConfig)


    def within_modeled_region(self, target):
        return (
            (target.x >= self.properties['modeled_region_xyz_start'][0]) &
            (target.x <= self.properties['modeled_region_xyz_end'][0]) &
            (target.y >= self.properties['modeled_region_xyz_start'][1]) &
            (target.y <= self.properties['modeled_region_xyz_end'][1]) &
            (target.z >= self.properties['modeled_region_xyz_start'][2]) &
            (target.z <= self.properties['modeled_region_xyz_end'][2])
        )

if __name__ == "__main__":
    model = OlfactoryBulb()
    model.build()
    model.setup_simulation()
    model.run()