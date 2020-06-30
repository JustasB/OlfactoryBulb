*********************************************
Recreating the Model
*********************************************

The model was created in stages: cell models, network models, odor input, LFP output, simulations, and experiments.

==========================================================================================================================
Mitral, tufted, and granule cell morphologies cleaned, morphology metrics computed, and representative morphologies chosen
==========================================================================================================================


 - Cell electrical properties and behavior identified and characterized
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