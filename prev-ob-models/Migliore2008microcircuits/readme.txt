NEURON mod files from the paper:

Michele Migliore and Gordon M Shepherd
Dendritic action potentials connect distributed dendrodendritic microcircuits
J. Comp. Neurosci. (2007) (online first).

The simulation file can be used to reproduce the bottom left histogram
in Fig.3,representing the M2 mitral cell firing in the presence of the
GCs network.  In the simulation window, "odor" is the number of odor
input (3-11), and traces for the 3 Mitral cells are plotted using
different colors.

To start, auto-launch from ModelDB or

Under unix systems:
to compile the mod files use the command 
nrnivmodl 
and run the simulation hoc file with the command 
nrngui mosinit.hoc

Under Windows systems:

to compile the mod files use the "mknrndll" program on the folder
expanded from this archive, then double click on the simulation file
mosinit.hoc to start the simulation.

Under MAC OS X:

Expand the MC-GC.zip archive after downloading it.  Drag the MC-GC
folder onto the mknrndll icon in the NEURON application folder.  Drag
the mosinit.hoc file in the MC-GC folder onto the nrngui icon in the
NEURON application folder.

For each platform:

Enter the odor number next to the odor button, press Run, and count
the number of M2 APs (red trace).  Repeat for odors 3 through 11.

Questions on how to use this model
should be directed to michele.migliore@pa.ibf.cnr.it
