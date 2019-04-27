## Readme for the juxtaglomerular models associated with the paper *Role of intraglomerular circuits in shaping temporally structured responses to dynamic olfactory input*

Ryan M Carey, William Erik Sherwood, Alla Borisyuk, Michael T. Shipley, Matt Wachowiak

### Requirements
All of the files included in this archive are MATLAB `.m` files (models) and `.mat` files (data). Some functions require the MATLAB Signal Processing Toolbox to run properly (mainly those for processing the raw calcium imaging data into ORN firing rates).

### Directory Organization

`startup.m`  simply adds all of the subdirectories of the ModelDB folder to the MATLAB path. Doing so is necessary for many of the functions to work properly.

`example.m` is a script with an example of how to run one of the models with default parameters and default inputs. It uses `doloop.m` to run the specified model once for each of the single glomerulus input traces in the specified ORN input data file. 

#### Data

`ModelDB/Data/` contains input data to be fed into the juxtaglomerular models. Subdirectories: 

* `Configuration/` contains initial conditions, gL, vL parameters for the grid sweep used to calibrate the ET cell baseline excitability. (This is a subset of the ET cell calibration file described below.)
* `Input/` contains input data for:
	* `ET_Calibration/` initial conditions for ET cell calibration 
	* `ORN_Input_Data/` raw ORN calcium imaging data (NAME) and ORN input signals preprocessed with bleaching removal and simulated ORN depression (NAME)
	* `Reference_Spikes/` voltage traces of single ET and PG action potentials used for calibrating synapse parameters. Note that these data are stored as Nx2 vectors with time values in the first column and voltage values in the second.
* `Parameters/` a collection of `.mat` files containing default parameters for individual components of the various juxtaglomerular models, including MC, MC with RI (MCGC), PG, and ET cells; excitatory (fast) and inhibitory (fast, slow) synapses with and without depression. SS140_pars.mat, SS170_pars.mat, etc. contain parameter values for the slow inhibitory PG-MC  synapse with decay times of 140 ms, 170 ms, and so forth.

#### DataProcessing
`ModelDB/DataProcessing/` contains the code for preprocessing the raw calcium signals contained in `raw_sniff_playback_orn_signals.mat`

#### Vfields
`ModelDB/Vfields/` contains `.m` files the vector fields describing each of the juxtaglomerular models, as well as each component neuron model individually. These vector fields are passed to the MATLAB integrators as function handles. 

#### Aux
  
The subdirectories of `ModelDB/Aux` (`ETaux/`, `MCaux/`, `PGaux/`, and `Synaux/`) contain `.m` files with "helper" functions for the vector fields in `Vfields`. These helper functions describe the channel kinetics for ion channels in the ET, MC, and PG cells, as well as activation functions for the various synapse models. They must be on the MATLAB path in order for the models to run properly.

#### Models

`ModelDB/Models/` contains `.m` files for functions that set up and integrate the various juxtaglomerulus models. In general, these functions take as input 
	
1. an input trace, 
2. the input trace sampling rate, 
3. a variable length argument list which is used to change specific model parameters. 

The functions call MATLAB's ode integrator on the appropriate vector field, with an event function set which will detect MC spike times (if applicable). The model files are grouped into subdirectories as follows:

* `ORN-ET` a single ET cell receiving ORN input
* `ORN-PG` a single PG cell receiving ORN input
* `ORN-MCRI` a single MC cell receiving ORN input, with optional recurrent inhibition. Used for the ORN-MC and ORN-MCRI models.
*  `ORN-ET-MC` contains two different models: `ORN-ET-MCRI.m` implements a feedforward excitatory circuit mediated by ET cells (only). `ORN-ET-MCRIpexcite.m` implements a direct ORN-MC input as well as the ET-MC connection.
*  `ORN-PG-MC` ORN-MCRI subcircuit with PG cell-mediated feedforward inhibition.
*  `Full` the full model combining the parallel excitatory feedforward circuit of `ORN-ET-MCRIpexcite.m` with one ORN-PG and one ORN-ET-PG inhibitory circuit.


#### Loops

`ModelDB/Loops` contains two helper function files:

* `doloop.m` contains code to run a selected model for on a collection of input traces (given as a single file name) with a given set of parameter settings, saving the output to a file. 
* `doloop_grid.m` is similar to `doloop.m`, but it runs the given model on the input traces, varying a specified pair of parameters. The values of the varying parameters must be given and should form a rectangular grid. 

The each of the subdirectories contains functions for using specific juxtaglomerular models with either `doloop.m` or `doloop_grid.m`. The naming convention is the same as for `ModelDB/Models`.

### Example
The file `example.m` demonstrates how to use the model files in `ModelDB/Models` along with the `doloop` function. Use of the `doloop_grid` function is very similar.