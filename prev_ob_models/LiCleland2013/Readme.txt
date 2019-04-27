NEURON CODE FOR IMPLEMENTING SIMULATIONS OF CHOLINERGIC
NEUROMODULATION IN OLFACTORY BULB

				Guoshi Li
				Computational Physiology Lab
				Department of Psychology
				Cornell University
				Ithaca, NY 14853

A full description of the model may be found in:

Li G and Cleland TA(2013) A two-layer biophysical model of cholinergic
neuromodulation in olfactory bulb.  Journal of Neuroscience
33:3037-3058.

For questions, please email: gl275@cornell.edu
	
The OB network model contains 25 mitral cells (MCs), 25 periglomerular
cells (PGs) and 100 granule cell (GCs)
 
The package contains five folders:

"celldata" folder to store data from single cell simulations
"data" folder to store data from network simulations 
"SP" folder to store the timing of random background spikes to MCs
"Connection" folder to store the connectivity information between MCs
     and GCs
"Input" folder contains two pre-generated data files read by the
     program to set the pre-odor and steady-state odor values

For network simulation, run mosinit.hoc
For MC single-cell simulation, run MC_Stim.hoc
For GC single-cell simulation, run GC_Stim.hoc
For PG single-cell simulation, run PG_Stim.hoc

The simulation step is set to 0.001 ms for results presented in the
Journal of Neuroscience Paper (Li and Cleland 2013)
Simulation time of the full model is about 3.5 hours for one run in a
workstation
Small time step is used to ensure accuracy of results
To reproduce identical results shown in the paper, make sure a time
step of 0.001 ms is used
Larger step could be used for testing purpose

Some key parameters:
NTCE = 0: Full model with EPL simulation
NTCE = 1: Glomerular model without EPL (i.e., no granule cells)
NICOTIN  = 0: nAChRs inactive
NICOTIN  = 1: nAChRs activated
MUSCARIN = 0: mAChRs inactive
MUSCARIN = 1: mAChRs activated	 

The default setting is the full model under the Control case (nAChRs
and mAChRs inactive)

The data saved in the "data" folder after simulation is analyzed using
the following custom Matalb scripts

PlotV.m: plot cell voltages
PlotG.m: plot GABAa conductances
Raster.m: generate raster plots of spikes
LFP.m: frequency analysis of the sLFP
Phase.m: generate phase distribution plot and raster plot of spike phases
