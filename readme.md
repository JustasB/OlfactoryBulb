# Getting Started

The following steps will run the simulations showing the gamma signature
and the results of experiments that elucidate the mechanisms underlying it.
 
This requires NEURON+Python+MPI. See this [tutorial to install NEURON+Python+MPI](https://neurojustas.com/2018/03/27/tutorial-installing-neuron-simulator-with-python-on-ubuntu-linux/).
 
 1. Clone the repository and install required packages
 ```
 git clone https://github.com/justasb/olfactorybulb
 cd olfactorybulb
 pip install -r requirements.txt
 ```
 
 2. Compile the .mod files
 ```
 cd olfactorybulb/prev_ob_models
 ./compile_mod.sh
 cd ..
 ```
 
 3. Edit the [repo]/runbatch.py file, replacing 16 on the 
 last line with the number of CPU cores on your machine.
 
 4. Start Jupyter Notebook
 ```
 jupyter lab "notebooks/LFP Wavelet Analysis.ipynb"
 ```
 Run all cells of notebook to run the simulations and see the analysis plots
 
 
# To run a single simulation of a parameter set
 
 Select the name of a parameter set class from a .py file under
 [repo]/olfactorybulb/paramsets/ (e.g. 'GammaSignature' from case_studies.py)
 
 Then run the simulation of the parameter set with:
 ```
 mpiexec -np 16 python initslice.py -paramset GammaSignature -mpi
 ```
 
 Should see output like:
 ```
 numprocs=16
 Rank Complexity min: 551, mean: 662.625, max: 791
 Starting simulation...
 Time: 13.0 ms
 ```
 
 When finished, the results are stored under [repo]/results/[parameter set class]

 Then run the following command in the above Jupyter notebook:
 ```
 show_plots('GammaSignature', sniff_count=8)
 ``` 
 
 
# Folders

**Folders needed to run the network model**

 - `olfactorybulb` Classes and database defining the network model
 - `prev_ob_models` Cell and network models developed by others, but cell models under `Birgiolas2020` are used in this network model. The models are compared against experimental data and to each other. 
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
 
 
