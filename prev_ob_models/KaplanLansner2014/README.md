OlfactorySystem
===============

This repository holds code that has been used for the publication:
Kaplan BA and Lansner A (2014) A spiking neural network model of self-organized pattern recognition in the early mammalian olfactory system. Front. Neural Circuits 8:5. doi: 10.3389/fncir.2014.00005

The paper is open access and can be found here: http://journal.frontiersin.org/Journal/10.3389/fncir.2014.00005/abstract

The simulation code for the paper has been run on a Cray super computer and will very likely not run on any personal computer (in the near future) due to memory and computation time requirements.


Olfactory system simulation 

  The neuron code is in the neuron_files subfolder.
  Before running a simulation run nrnivmodl in neuron_files folder.

  For the distribution of Olfactory Receptor affinities:
    The data in the Haddad_data folder is from http://www.nature.com/nmeth/journal/v5/n5/extref/nmeth.1197-S3.xls [Haddad 2008 "A metric for odorant comparison", Nature Methods] 
    The script cluster_odorant_space.py computes the distance between the
    real-world data and the virtual olfactory receptors (=centroids after
    k-means clustering) for many trials and for various numbers of centroids.
    
    Distances between ORs and odorants is pooled by
    average_OR_affinity_distributions.py over many trials.
    It tries to fit a distribution to the data and writes the fit parameters
    to a file, which can be displayed by plot_OR_placement_fit_params.py.
    
     
  ORN response curves: run_epth_response_curve.py (calls
  prepare_epth_response_curve.py and requires SetOfCurvesPlotter.py, and
  MergeSpikefiles.py and the NEURON files --> start_file_epth_response_curve.hoc)

  OB response curve: similar to ORN response curve measurements, but with ORNs
  projecting to the OB 
  
  
  For the system simulations:
  First, run prepare_epth_ob_prelearning.py, then run_epth_ob_prelearning.py.
  
  After that, create a new folder e.g. use the old folder name ("Folder") with
  the post fix "postLearning" (--> "Folder_postLearning"), i.e. in
  simulation_parameters.py set_folder_name function:
    params['folder_name'] = 'Folder_postLearning'
  and after that, create a new folder structure by doing:
    python CreateObOcConnections.py new
  Make sure the new folders are created as you expect.
    
  Then, copy all files storing the number of spikes fired by mitral cells to the new folder:
    cp Folder/NumberOfSpikes/mit_nspikes_* Folder_postLearning/NumberOfSpikes/
    
  Then, do
    python CreateObOcConnections.py 
  and 
    python run_full_system_recognition_task.py
  
  
  To run the system without resimulating the epithelium and the OB, again in
  simulation_parameters.py set 
    params['folder_name'] = 'Folder_OcOnly'
    
  Then, copy all files storing the merged mitral cell spiketimes to the new folder:
  Merge all the spike time files:
    python MergeSpikefiles.py Folder/ mit
  Copy to new folder:
    cp Folder/Spiketimes/mit_spiketimes_merged* Folder_OcOnly/Spiketimes/
    cp Folder/NumberOfSpikes/mit_nspikes_merged* Folder_OcOnly/NumberOfSpikes/
  Prepare and create connections:
    python CreateObOcConnections.py 
    
  

Lindgren specific information

  Compilation:
  bkaplan@emil-login2:/cfs/klemming/nobackup/b/bkaplan/OlfactorySystem/neuron_files>
  /cfs/klemming/nobackup/b/bkaplan/neuron-7.3/nrn-mpi/bin/nrnivmodl




