NEURON mod files from the paper:
Activity-Dependent Adjustments of the Inhibitory Network 
in the Olfactory Bulb Following Early Postnatal Deprivation,
by Saghatelyan, A., Roux, P., Migliore, M., Rochefort, C., Desmaisons,
D., Charneau, P., Shepherd, G.M., and Lledo, P.M. 
Neuron 46:103-116, 2005.

The paper investigates the role of experience in 
sculpting neuronal network in the Olfactory bulb.
It was found that nostril closure decreased the number of
newborn granule cells of the MOB, the complexity of their dendritic arborization,
and spine density. However, due to a compensatory increase in the newborn 
granule cells excitability, the action potential-dependent GABA release was dramatically 
enhanced, counteracting thus the reduction in the spine density and leading to an 
unaltered synchronization of mitral cells firing activity.
The model supports the experimental findings, and shows that a -10mV shift in the
Na activation or a reduction in the dendritic lenght of newborn GC
could independenlty explain the observed increase in excitability. 

The effect is demonstrated here by reproducing the simulations
in Figs.6 and 8 of the paper.
3 simulations can be run for different kind of model setup:
- control conditions 
- shift in the GC Na activation. 
- reduced dendritic length of GC

Under unix systems:
to compile the mod files use the command 
nrnivmodl 
and run the simulation hoc file with the command 
nrngui occlusion.hoc

Under Windows systems:
to compile the mod files use the "mknrndll" command.
A double click on the simulation file
occlusion.hoc 
will open the simulation window.

Questions on how to use this model
should be directed to michele.migliore@pa.ibf.cnr.it


