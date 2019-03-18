NEURON mod files for the I-A and I-K currents from the paper:
Wang, et al., Whole-cell K+ currents in identified olfactory 
bulb output neurones of rats. J Physiol. 1996 490:63-77.

Running the kinetics.hoc simulation file will show 
the activation and inactivation steady-states, the time constants, 
and a family of curves generated modeling the same protocols 
used for Figs.6B-C of the paper.

Markers show the experimental values reported in the paper (Table.1).
For the I-A current, two parameters (shi_kamt and sha_kamt) 
have been introduced to take into account the shift in the kinetics (Fig.4B) 
observed using different bathing solutions 
(with or without Ca++, see legend of Fig.4). 
Default values for shi_kamt (inactivation) and sha_kamt (activation)
are 5.7 and 9.9mV, respectively. 
shi_kamt = sha_kamt = 0 was used here (solid lines in Fig.4B).
 
Under unix systems:
to compile the mod files use the command 
nrnivmodl 
and run the simulation hoc file with the command 
nrngui kinetics.hoc

Under Windows using NEURON 5.1:
to compile the mod files use the "mknrndll" command.
A double click on the simulation file
kinetics.hoc 
will open the simulation window.

Questions on the model parameters should be directed to the 
authors.

Questions on how to use this model with NEURON
should be directed to michele.migliore@pa.ibf.cnr.it
