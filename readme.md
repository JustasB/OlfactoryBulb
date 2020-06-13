# Folders

**Folders needed to run the network model**

 - `olfactorybulb` Classes and database defining the network model
 - `prev_ob_models` Cell and network models developed by others. The models are compared against experimental data and to each other. Also contains the cell models (under `Birgiolas2020`) used in this network model
 - `blender-files` Stores .blend used in network construction (e.g. layer coordinates)
 - `initslice.py` Runs the network model using a provided parameter set
 - `docker` Docker scripts to run the model using Docker 
   
   
**Folders used to construct the network model and cell models**

 - `digitized-figures` Extracted figures that contained experimental data used in the model
 - `morphology-data` Subfolders with .SWC morphology archives from [NeuroMorpho.org](http://neuromorpho.org/) of the three cell types
 - `neuronunit` [NeuronUnit](https://github.com/scidash/neuronunit/) classes that define tests used to validate cell models
 - `notebooks` Jupyter notebooks used to validate, fit, and simulate the cell and network models
 - `worksheets` Excel worksheets used to derive experimental data properties when they were not directlya vailable
  
**Other folders**
 - `notes` Mostly notes and temporary model validation files
 - `dissertation-figures-tables` Excel spreadsheets of some dissertation tables
 - `media` Images of dissertation figures and videos used in dissertation defense slides
  
# Other Files

 - `prev_ob_models/compile_mod.sh` Compiles all .mod files 
 - `runbatch.py` Allows specifying different parameter sets and runs the model with each set
 - `testmpi.py` A test of NEURONs MPI-based network
 