NEURON mod files from the paper:
M. Migliore, M. Hines and GM Shepherd
The role of distal dendritic gap junctions in synchronization
of mitral cell axonal output, J. Comput. Neurosci., 18:151-161 (2005).

In a realistic model of two electrically connected mitral cells (M1, M2),
the paper shows that the somatically-measured experimental properties 
of GJs may correspond to a variety of different local coupling strengths 
and dendritic distributions of GJs in the tuft. The model suggests 
that the propagation of the GJ-induced local tuft depolarization
is a major mechanim for intraglomerular synchronization of mitral cells.

The main effect is demonstrated here by reproducing the simulations
in Fig.6A-B-C of the paper.

In each case the somatic membrane potential of M1 and M2
and the gap current are shown for the first 300ms of simulation.
At the end of each run the cross-correlation function for 
the last 150ms of simulation is calculated and plotted in an interactive window.

Under unix systems:
to compile the mod files use the command 
nrnivmodl 
and run the simulation hoc file with the command 
nrngui forfig6-modeldb.hoc

Under Windows systems:
to compile the mod files use the "mknrndll" command.
A double click on the simulation file
forfig6-modeldb.hoc 
will open the simulation window.

Questions on how to use this model
should be directed to michele.migliore@pa.ibf.cnr.it


