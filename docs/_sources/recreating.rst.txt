*********************************************
Recreating the Model
*********************************************

The model was created in stages: cell models, network models, odor input, LFP output, simulations, and experiments. While these stages are described in the `dissertation <https://repository.asu.edu/attachments/223567/content/Birgiolas_asu_0010E_19503.pdf>`_, the sections below describe technical details and include links to specific files, tools, and resources. The sections are organized in a chronological order, which can be followed to re-create this or a similar model.

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

For each reconstruction, NeuroMorpho.org displays a set of whole-cell morphology metrics (e.g. see 'Measurements' section of `this reconstruction <http://neuromorpho.org/neuron_info.jsp?neuron_name=1-09-TD2b>`_). For cells in this model, those metrics were computed for soma, apical dendrites, and lateral dendrites of each reconstruction.

The metrics were computed using the `pyLMeasure tool <https://github.com/justasb/pylmeasure>`_ (which is a Python wrapper of the `L-Measure <http://cng.gmu.edu:8080/Lm/>`_ tool). Each metric was formalized as tests of the `NeuronUnit framework <https://github.com/scidash/neuronunit>`_. The test definitions have been incorporated into NeuronUnit and can be found under `[neuronunit]/tests/morphology.py <https://github.com/scidash/neuronunit/blob/dev/neuronunit/tests/morphology.py>`_

The morphology metrics of each cell were computed using suites of NeuronUnit tests (see: `[repo]/notebooks/morphology.ipynb <https://github.com/JustasB/OlfactoryBulb/blob/master/notebooks/morphology.ipynb>`_) and results stored in `the model's SQLite Database <https://github.com/JustasB/OlfactoryBulb/blob/master/olfactorybulb/model-data.sqlite>`_ in the ``measurement`` table. The database tables and records can be viewed online with tools like `SQLite online <https://sqliteonline.com/>`_ or offline with an SQLite editor like `SQLite Expert <http://sqliteexpert.com/>`_.

The overall means and standard deviations of each morphology metric were stored in the database in the ``property`` table. These overall means and standard deviations were used to select five stereotypical morphologies of each cell type. Reconstructions whose morphology metrics were closest to the cell type population means were selected to be used for cell models in this model. They can be found under `[repo]/prev_ob_models/Birgiolas2020/SWCs <https://github.com/JustasB/OlfactoryBulb/tree/master/prev_ob_models/Birgiolas2020/SWCs>`_.

The validation of the representativeness of the selected morphologies and comparison with previously published models can be found in the `morphology-validation.ipynb notebook <https://github.com/JustasB/OlfactoryBulb/blob/master/notebooks/morphology-validation.ipynb>`_

**Morphology Standardization to Prepare for Network Embedding**

To facilitate the embedding of the cell models within reconstructed network layers, the cell morphologies were further edited as follows. The cells were reoriented so that their apical dendrites were aligned with the Z-axis. The natural curvature of lateral dendrites of mitral and tufted cells was planarized so that the lateral dendrites were essentially perpendicular to the apical dendrites. This prepared the lateral dendrites for an algorithm to confine them to the appropriate cell-type-specific laminar regions in the network model. These standardized morphologies can be found under `[repo]/prev_ob_models/Birgiolas2020/morphology <https://github.com/JustasB/OlfactoryBulb/tree/master/prev_ob_models/Birgiolas2020/morphology>`_.

====================================================================
Cell electrical properties and behavior identified and characterized
====================================================================

Once cell morphologies were selected, they were turned into active models by inserting ion channel mechanisms. To understand which ion channels should be included, a survey of electrical behaviors of each cell type was conducted via a literature search.

**Generic Test and Publication Python Classes**

Electrical properties of mitral, tufted, and granule cells from main olfactory bulbs of adult mice were aggregated from literature. These properties included passive (e.g. input resistance and resting voltage) and active membrane properties (e.g. spike half-width, spike accommodation, and inhibitory rebound). The protocols used to measure and compute each property were formalized using a set of NeuronUnit tests. Each publication used slightly different protocols and algorithms to compute the properties (e.g. different slice temperatures, different delays to measure input resistance, different ways to compute action potential threshold). To fairly evaluate how a model reproduces a property, these protocol differences were formalized into the NeuronUnit tests as well. A given property reported in a publication was formalized as a NeuronUnit test that consisted of 1) a generic electrophysiology property Python class (e.g. `RestingVoltageTest <https://github.com/JustasB/OlfactoryBulb/blob/master/olfactorybulb/neuronunit/tests/tests.py#L53>`_) and 2) publication-specific modifications of the property measurement protocol (e.g. `BurtonUrban2014 <https://github.com/JustasB/OlfactoryBulb/blob/master/olfactorybulb/neuronunit/tests/publications.py#L34>`_).

A NeuronUnit test that computes a property using protocol variations from a specific publication is declared by combining the generic test class with a publication class using Python's multiple inheritance (e.g. ``class InputResistanceUrbanBurton2014Test(UrbanBurton2014, InputResistanceTest)`` ). Due to the large number of property and publication combinations, these combined classes were created at run time using `a helper method <https://github.com/JustasB/OlfactoryBulb/blob/master/olfactorybulb/neuronunit/tests/__init__.py#L73>`_.

The measurements from each publication are recorded in the SQLite database, ``measurement`` table. The NeuronUnit class to use for each measurement is defined in the ``property`` table and ``test_class_generic`` column. At run-time, the helper method combined the generic test class with the publication class to create a class that computed the property using the protocol parameters specified in each publication.

=======================================================================================================
Ion channels placed onto chosen cell morphologies and conductances fitted to experimental distributions
=======================================================================================================

**Previously Published Olfactory Models**

A literature search to identify each cell type's electrical properties revealed a set of ion channels likely responsible for the behavior. Furthermore, a separate literature search was performed to identify previous computational models of each cell type. `ModelDB <https://senselab.med.yale.edu/ModelDB/default>`_ was utilized to locate runnable NEURON simulator models. Results from the initial search of all olfactory bulb or cell models were manually inspected by following the instructions to run the models and then inspecting which cell types were modeled. If a model included multi-compartment morphologies, the morphology was extracted into .SWC files using the `hoc2swc <https://github.com/JustasB/hoc2swc>`_ tool. The set of previous olfactory models can be found under `[repo]/prev_ob_models <https://github.com/JustasB/OlfactoryBulb/tree/master/prev_ob_models>`_.

**Selection of Ion Channel Models**

In addition to morphologies, the previous models included ion channel .MOD files. Many of the .MOD files were nearly identical copies of files from earlier publications. After tracing the publication geneology of each .MOD file, the most recent versions were identified. By examining the reasons why the channels were included, and comparing the electrical behavior seen in experimental literature, channels hypothesized to be responsible for the experimentally observed behaviors were selected for inclusion in this model. If an ion channel model did not include temperature sensitivity, equations to utilize the Q10 coefficient were added to the .MOD file. The ion channel models used in this model can be found under `[repo]/prev_ob_models/Birgiolas2020/Mechanisms <https://github.com/JustasB/OlfactoryBulb/tree/master/prev_ob_models/Birgiolas2020/Mechanisms>`_.

NEURON Cell Builder feature was used to import morphology SWC files and then insert channel mechanisms. The initial, unoptimized cell models can be found under `[repo]/prev_ob_models/Birgiolas2020/Cells <https://github.com/JustasB/OlfactoryBulb/tree/master/prev_ob_models/Birgiolas2020/Cells>`_.

**Parameter Fitting and Validation**

A genetic optimization algorithm was used to identify the combinations of ion channel parameters that best reproduced experimentally observed electrical behaviors (minimized the error between model electrical property values and experimentally observed property means). The algorithm is defined in `fitting.py <https://github.com/JustasB/OlfactoryBulb/blob/master/prev_ob_models/Birgiolas2020/fitting.py>`_ and Jupyter notebooks used for fitting can be found under `[repo]/notebooks <https://github.com/JustasB/OlfactoryBulb/tree/master/notebooks>`_ (files that start with ``fitting``).

The validation of the cell models against experimental data and comparison to previously published models can be found in `all-model-validation-results.ipynb notebook <https://github.com/JustasB/OlfactoryBulb/blob/master/notebooks/all-model-validation-results.ipynb>`_.

The final, optimized models and their parameters are defined in `[repo]/prev_ob_models/Birgiolas2020/isolated_cells.py <https://github.com/JustasB/OlfactoryBulb/blob/master/prev_ob_models/Birgiolas2020/isolated_cells.py>`_ (see MC1..5, GC1..5, TC1..5 classes).

=========================================================================
Olfactory bulb layers were reconstructed from sagittal and coronal slices
=========================================================================

====================================================================================
Fitted cell somas were placed within each cell type's stereotypical laminar location
====================================================================================

======================================================================
Apical dendrites were rotated towards their stereotypical terminations
======================================================================

===============================================================================================================
Mitral and tufted cell lateral dendrites were aligned with the curvature of reconstructed olfactory bulb layers
===============================================================================================================

=========================================================================================================
Reciprocal synapses were formed based on dendritic proximity between principal and granule cell dendrites
=========================================================================================================

====================================================================================
Gap junctions were formed between glomerular sibling principal cell tufted dendrites
====================================================================================

=================================================================================
Glomeruli stimulated during odor experiments were mapped onto the model glomeruli
=================================================================================

=====================================================================
A model of glomerular input spikes was created to stimulate the model
=====================================================================

=======================================================================================
An extracellular local field potential electrode was placed into the granule cell layer
=======================================================================================

=====================================================================================
NEURON+MPI simulations were performed and LFP signal analyzed using wavelet transform
=====================================================================================

=======================================================================================
Network parameters were explored until the two-cluster gamma fingerprint was reproduced
=======================================================================================

=======================================================================================================
Computational experiments were performed to demonstrate the mechanisms underlying the gamma fingerprint
=======================================================================================================

=======
Outline
=======

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