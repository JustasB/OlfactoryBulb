def build_slice():
    """
    To build the model of the slice, most of the work is performed in Blender.

    NEURON is used to instantiate cells, which are exported to Blender, where
    they are positioned and their morphologies modified. These modifications
    are saved to files that NEURON can load later, to run the simulation.

    This file serves as the launcher of NEURON+Blender. It starts by launching
    NEURON with its part of the BlenderNEURON. A few helper methods are added
    to NEURON that can be called from Blender.

    Then, once NEURON is running, in parallel, Blender is started with the BlenderNEURON
    addon. Blender imports cells instantiated in NEURON and uses Blender functions
    to manipulate their morphology.

    Once the cells are positioned, they are saved into files that NEURON can use to
    load the slice model.

    Once the script starts, monitor the console output for progress. After all the cells are
    positioned, connected, and saved, the Blender window will open, showing the model.
    """


    import os
    from olfactorybulb.slicebuilder.nrn import SliceBuilderNRN

    # Start NRN and the addon
    sbn = SliceBuilderNRN()

    # Start Blender and build the model
    os.system("blender blender-files/ob-gloms-fast.blend --python olfactorybulb/slicebuilder/blender.py")

if __name__ == '__main__':
    build_slice()