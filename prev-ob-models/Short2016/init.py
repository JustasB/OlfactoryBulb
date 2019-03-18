#!/usr/bin/python

"""init.py 
Relies on pre_init.py having been run and then the archive compressed
(zip'ed) and uploaded to the NSG.  This file will start a parallel
context, forms three relevant file/folder variables (the relative
paths to and then loads) the num_of_columns.hoc and parameters.hoc,
and the path to tdt2mat_data). Then the build_net_Shep_NSG.hoc which runs
the simulations.  Subsequently the tanks are written in the
appropriate location.

"""

from mpi4py import MPI
from neuron import h, gui
pc = h.ParallelContext()
id = int(pc.id())
nhost = int(pc.nhost())
print "I am", id, "of", nhost

# form the three relevant files/folders based on id
if id<18: # 2888:
  run_folder = "run_%d"%(id) # helper folder name variable
  parameters_dot_hoc_file = run_folder+"/parameters.hoc"
  num_of_columns_dot_hoc_file = run_folder+"/num_of_columns.hoc"
  #tank_folder = "tdt2mat_data/"
  tank_folder = run_folder+"\/tdt2mat_data\/"
  
  # make these variables exist in NEURON so the approp. files and can be loaded/saved
  # The first two are used in buil_net_Shep_NSG.hoc and the tank_folder is used
  # in this NSG version of tdt2mat_data.hoc
  
  h("strdef parameters_dot_hoc_file")
  h("strdef num_of_columns_dot_hoc_file")
  h("strdef tank_folder")
  
  h('parameters_dot_hoc_file="%s"'%(parameters_dot_hoc_file))
  h('num_of_columns_dot_hoc_file="%s"'%(num_of_columns_dot_hoc_file))
  h('tank_folder="%s"'%(tank_folder))
  
  # create the simulation model structure, run and store
  h.load_file("build_net_Shep_NSG.hoc")
  
  
