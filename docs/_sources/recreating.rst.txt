*********************************************
Recreating the Model
*********************************************

The model was created in stages: cell models, network models, odor input, LFP output, simulations, and experiments.

==========================================================================================================================
Mitral, tufted, and granule cell morphologies cleaned, morphology metrics computed, and representative morphologies chosen
==========================================================================================================================

**Morphology Archives**

Morphologies for each cell type (mitral, tufted, granule) were obtained from `NeuroMorpho.org <http://neuromorpho.org/>`_, by performing a `'Metadata search' <http://neuromorpho.org/MetaData.jsp>`_ of reconstruction archives. The archives were filtered to include cells with the following parameters: mouse, adult, olfactory bulb region, control condition, and mitral, tufted, or granule cell type.

**Inspection and Editing of the SWC Files**

Each resulting archive was downloaded and organized into folders which can be found under `[repo]/morphology-data <https://github.com/JustasB/OlfactoryBulb/tree/master/morphology-data>`_. The 3D renderings of the .SWC files under each archive were manually examined using the `neuTube <https://www.neutracing.com/>`_ software. neuTube functions were used to fix simple and obvious flaws (e.g. parts of dendritic trees shifted laterally, or all sections mislabeled as 'basal' dendrites). Archives with cells missing important stereotypical features of each cell type were excluded. Besides, neuTube, other `neuronal morphology/SWC editing software <https://neurojustas.com/2019/03/10/tools-for-editing-swc-files-neuron-morphology/>`_ was evaluated, but ultimately not used.

Each SWC archive was labeled for whether it was a 'thin' (parts of lateral dendrites missing) or 'full' (the whole cell was reconstructed) slice. Whether the cells included a soma reconstruction, and whether the reconstructed dendrites included width/radius information (e.g. some archives used constant radius).

The fixed SWC files were saved back into the archive folders and archive categorizations were stored in the `[repo]/notebooks/morphology.ipynb Jupyter Notebook <https://github.com/JustasB/OlfactoryBulb/blob/master/notebooks/morphology.ipynb>`_.

**Computing Reconstructed Cell Population Morphology Metrics**

For each reconstruction, NeuroMorpho.org displays a set of whole cell morphology metrics (e.g. see 'Measurements' section of `this reconstruction <http://neuromorpho.org/neuron_info.jsp?neuron_name=1-09-TD2b>`_). For cells in this model, those metrics were computed for soma, apical dendrites, and lateral dendrites of each reconstruction.

The metrics were computed using `pyLMeasure tool <https://github.com/justasb/pylmeasure>`_ (which is a Python wrapper of the `L-Measure <http://cng.gmu.edu:8080/Lm/>`_ tool). Each metric was formalized as tests of the `NeuronUnit framework <https://github.com/scidash/neuronunit>`_. The test definitions have been incorporated into NeuronUnit and can be found under `[neuronunit]/tests/morphology.py <https://github.com/scidash/neuronunit/blob/dev/neuronunit/tests/morphology.py>`_

The morphology metrics of each cell were computed using suites of NeuronUnit tests (see: `[repo]/notebooks/morphology.ipynb <https://github.com/JustasB/OlfactoryBulb/blob/master/notebooks/morphology.ipynb>`_) and results stored in `the model's SQLite Database <https://github.com/JustasB/OlfactoryBulb/blob/master/olfactorybulb/model-data.sqlite>`_ in the ``measurement`` table. The database tables and records can be viewed with an SQLite editor like `SQLite Expert <http://sqliteexpert.com/>`_.

The overall means and standard deviations of each morphology metric were stored in the database in the ``property`` table. These overall means and standard deviations were used to select five stereotypical morphologies of each cell type. Reconstructions whose morphology metrics were closest to the cell type population means were selected to be used for cell models in this model. They can be found under `[repo]/prev_ob_models/Birgiolas2020/SWCs <https://github.com/JustasB/OlfactoryBulb/tree/master/prev_ob_models/Birgiolas2020/SWCs>`_.

The validation of representativeness of the selected morphologies can be found in the `morphology-validation.ipynb notebook <https://github.com/JustasB/OlfactoryBulb/blob/master/notebooks/morphology-validation.ipynb>`_

**Morphology Standardization to Prepare for Network Embedding**

To facilitate the embedding of the cell models within reconstructed network layers, the cell morphologies were further edited as follows. The cells were reoriented so that their apical dendrites were aligned with the Z-axis. The natural curvature of lateral dendrites of mitral and tufted cells was planarized so that the lateral dendrites were essentially perpendicular to the apical dendrites. This prepared the lateral dendrites for an algorithm to confine them to the appropriate cell-type-specific laminar regions in the network model. These standardized morphologies can be found under `[repo]/prev_ob_models/Birgiolas2020/morphology <https://github.com/JustasB/OlfactoryBulb/tree/master/prev_ob_models/Birgiolas2020/morphology>`_.

====================================================================
Cell electrical properties and behavior identified and characterized
====================================================================

 - Ion channels placed onto chosen cell morphologies and conductances fitted to experimental distributions
 - Olfactory bulb layers were reconstructed from sagittal and coronal slices
 - Fitted cell somas were placed within each cell type's stereotypical laminar location
 - Apical dendrites were rotated towards their stereotypical terminations
 - Mitral and tufted cell lateral dendrites were aligned with the curvature of reconstructed olfactory bulb layers
 - Reciprocal synapses were formed based on dendritic proximity between principal and granule cell dendrites
 - Gap junctions were formed between glomerular sibling principal cell tufted dendrites
 - Glomeruli stimulated during odor experiments were mapped onto the model glomeruli
 - A model of glomerular input spikes was created to stimulate the model
 - An extracellular local field potential electrode was placed into the granule cell layer
 - NEURON+MPI simulations were performed and LFP signal analyzed using wavelet transform
 - Network parameters were explored until the two-cluster gamma fingerprint was reproduced
 - Computational experiments were performed to demonstrate the mechanisms underlying the gamma fingerprint

 - cells
    - morphology
        - SWC archives
        - quality / cleaning
        - morphology selection
        - validation

    - electrophysiology properties and database
    - ion channels
    - electrophysiology property tests
    - fitting
    - comparison to other models

 - layers
    - reconstruction
    - mesh simplification

 - cell placement
    - placement within layers
    - orientation
    - dendritic alignment

 - synapses
    - M/TC <-> GC Reciprocal Synapses
        - AMPA/NMDA
        - GABA
        - spines
        - dendritic proximity rule

    - M/TC Gap Junctions

 - odor input
    - Migliore14 glomerular activation maps
    - Glomerular registration
    - Glomerular intensity to spikes model
    - Input connections

 - simulation
    - single thread
    - mpi
        - cell rank assignment
        - gap junction considerations
        - synaptic connections

 - recordings
    - somatic
    - extracellular lfp
    - wavelet analysis

 - experiments
    - silent network
    - Only MCs or TCs
    - Added GCs
    - Added Gap Junctions