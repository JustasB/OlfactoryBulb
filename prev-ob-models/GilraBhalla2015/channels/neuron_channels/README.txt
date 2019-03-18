-> These .mod files have been modified:
1) Their suffix (channel name) has been changed to match those in my (Aditya's) olfactory bulb model.
2) The inf-s and tau-s have been put into the form xinf anf xtau, etc. eg. minf and mtau.
3) Ensure that in the NEURON {} block, xinf and xtau are defined as GLOBAL.

-> I also need to load mit_memb.hoc to insert the tables from files required for kslow and kfast.
1) I kept only the table loading part of the file.
2) I modified the table_tabkinf_kslowtab to table_tabkinf_kslowtab, and similarly for kfasttab.

Check out NeuronSimulatorChannelTest.py and MOOSEChannelTest.py for plots of inf-s and tau-s.
First run nrnivmodl in the neuron_channels directory.
