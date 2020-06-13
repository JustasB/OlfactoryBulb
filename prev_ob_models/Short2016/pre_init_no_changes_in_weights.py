
#!/usr/bin/python

"""pre_init.py will create multiple folders for simple_circuit model
runs with names like run_X where X goes from 1 to the number of jobs
and places scripts and subfolders in them.

Workflow:
In general pre_init.py is run with a command like ./pre_init.py on my
local computer (laptop).  Then the archive is zipped up and sent to
NSG where it is run from the file init.py. Finally the results are
uploaded to louise for subsequent matlab analysis.

Only the files create_array.py and this file, pre_init.py need to be
edited to the parameters of the simulation.

phase 1 the run_X folders are created . In addition a tdt2mat_data subfolder is created.
This phase might be easiest done before the parallel job is submitted. 

phase 2 the num_of_columns.hoc and parameters.hoc files are created in
each run_x folder so that they are all ready to run with their
different assigned values.  A parameters[] list prepares the contents of
each of these files before they are written.

phase 3 The parallel job is uploaded to NGS: the mod files are
compiled and the jobs are started by each run worker executing
build_net_Shep.hoc and loading the appropriate run_X/parameters.hoc
and run_X/num_of_columns.hoc scripts. These run_workers automatically
save the neuron tank data (and optionally spike time data used to
create subsequent graphs).

phase 4 After the simulation completes, the simulation data is
uploaded to louise where the matlab tanks are created

phase 5 On louise, the tanks are analyzed to create and save matlab
polar plot figures

phase 6 raster plots of spike activity are created and
saved if desired.
"""


#************************************************************************
print "phase 1"
#************************************************************************
# phase 1 the run_X folders are created . In addition a tdt2mat_data subfolder is created.
# This phase might be easiest done before the parallel job is submitted. 

# the task of create arrays is to create "both" which has both the breathing peak rate
# and the stimulus (light) peak rate in list of (B, S) tuples.  Once this is created
# we can determine the number of jobs by multiplying the length of both by the number
# of network models, e.g. three for three different types of network models (full_net,
# pg_net, gc_net)
execfile("create_arrays.py")
jobs_per_nn=len(both)
num_of_columns_tested=len(add_columns) # add_columns supplied by create_arrays.py
net_type=['full_net','pg_net','gc_net']
num_of_nn_types=len(net_type) # just pg and noinhib rerun # 4 # pg, gc, full, noinhib
num_of_jobs=jobs_per_nn*num_of_nn_types*num_of_columns_tested
print "num_of_jobs: "+ str(num_of_jobs)

# http://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python

import os, errno

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            print "the folder ",path," seems to already exist - it will be cleared ********"
            pass
        else: raise

for folder_num in range(num_of_jobs):
  new_folder="run_"+str(folder_num)
  mkdir_p(new_folder)

# now make sure the folders are empty
# http://stackoverflow.com/questions/185936/delete-folder-contents-in-python

import os, shutil
for folder_num in range(num_of_jobs):
  # folder = '/path/to/folder'
  folder="run_"+str(folder_num)
  for the_file in os.listdir(folder):
    file_path = os.path.join(folder, the_file)
    try:
        if os.path.isfile(file_path):
            os.unlink(file_path)
        elif os.path.isdir(file_path): shutil.rmtree(file_path)
    except Exception, e:
        print e

# now add the tdt2mat_data subfolder
for folder_num in range(num_of_jobs):
  new_folder="run_"+str(folder_num)
  tdt_folder=new_folder+"/tdt2mat_data"
  mkdir_p(tdt_folder)


#************************************************************************
print "phase 2"
#************************************************************************
# In phase 2 the num_of_columns.hoc and parameters.hoc files are created in each run_x folder so that they
# are all ready to run with their different assigned values.

# num_of_columns.hoc setting
# set the total_num_of_columns_master to the value desired to be used for each of all the simulations

# use below if for simple case of two columns total - otherwise rely on columns supplied from create_arrays.py
#total_num_of_columns_master=2
#num_of_additional_columns = total_num_of_columns_master - 1

# note that since the hoc code uses num_of_columns to set the number of columns in addition to the "recorded"
# mitral cell column, the setting of num_of_columns should be to num_of_additional_columns

# instead of using the below code assign num_of_columns[] list to write the 
# num_of_columns.hoc files when the parameters.hoc files are being written,
# that way, it can be coordinated with the network types

#for folder_num in range(num_of_jobs):
#  folder="run_"+str(folder_num)
#  fid=open(folder+"/num_of_columns.hoc","w")
#  fid.write("n = "+str(num_of_additional_columns)+" // n easier to type than num_of_cols\n")
#  fid.close()

# parameters.hoc settings
# includes both setting parameters and running functions that copy parameters to all the columns
# functions:
#   adjust_netcons_from_top() copies all netcons from [X][0] to [X][Y>0]
#   toggle_gc_connection() toggles all columns gc netcons 0/gc_on
#   toggle_pg_connection() toggles all columns pg netcons 0/pg_on
#   where gc_on and pg_on have a default value of 1
#
# parameter dictionary parameters holds all values for simulations
parameters={}
columns={}
# helper lists
# B for breathing rates
# 0 to 620 in increments of 20

# some typical simulation run number assignments:  (is this useful?  maybe delete?)
# let 0 through 63 be pg mediated inhibition with nc[15][]=0 and
# 64 through 127 be  pg mediated inhibition with nc[15][]=1 and
# 128 through 191 be gc mediated inhibition
# with jobs_per_nn becomes
# let 0 through jobs_per_nn-1 be pg mediated inhibition with nc[15][]=0 and
# jobs_per_nn through 2*jobs_per_nn-1 be  pg mediated inhibition with nc[15][]=1 and
# 2*jobs_per_nn through 3*jobs_per_nn-1 be gc mediated inhibition


#for folder_num in range(jobs_per_nn):  # pg mediated inhibition with nc[15][]=0
#  parameters[folder_num]="""
#breathing_period=400
#light_period=398 // 300 is a short run, 398 regular
#breath_peak_rate = %d
#light1_peak_rate = %d
#light2_peak_rate = 0
#breath_half_width=20
#light_half_width=20
#for i=0, n-1 {
#  nc[15][i].weight = 0 // turn off breathing input to pg cells
#}
#toggle_gc_connection() // turns off all gc cell connections
#objref pwm
#pwm=new PWManager()
#pwm.hide(3) // close voltage window for faster run
#do_everything()
#quit()
#""" % both[folder_num]


jobs_per_nn=len(both)
num_of_columns_tested=len(add_columns) # add_columns supplied by create_arrays.p
c_bs_group_len = jobs_per_nn * num_of_columns_tested
parameters=[]
columns=[]
for n_index in range(len(net_type)):
  for c_index in range(num_of_columns_tested):
    for bs_index in range(jobs_per_nn):
      folder_num=n_index*c_bs_group_len+c_index*jobs_per_nn+bs_index
      print "net: "+net_type[n_index]+", run_"+str(folder_num)+", B="+str(both[bs_index][0])+", S="+str(both[bs_index][1])+", cols="+str(add_columns[c_index])
      # after the first one is appended it can be accessed by index
      parameters.append("""
breathing_period=400
light_period=399 // 300 is a short run, 398 regular
breath_peak_rate = %d
light1_peak_rate = %d
light2_peak_rate = 0
breath_half_width=30
light_half_width=30
"""%both[bs_index]) # note that both has tuple pairs of Breath and Stimulus peak_rate in list
      if net_type[n_index] in ['gc_net', 'noinhib_net']:
        parameters[folder_num] = parameters[folder_num] + \
"""// for some reason toggle_pg_connection() was causing an error however
// the below worked
    for i=0, n-1 {
      nc[14][i].weight = 0
      nc[15][i].weight = 0
      nc[16][i].weight = 0
      nc[17][i].weight = 0
      nc[18][i].weight = 0
      nc[19][i].weight = 0
      nc[20][i].weight = 0
      nc[21][i].weight = 0
      nc[22][i].weight = 0
      nc[23][i].weight = 0
      nc[24][i].weight = 0
      nc[25][i].weight = 0
    // xstatebutton automatically sets pg_connection_state=0
    }
// toggle_pg_connection() // turns off all pg cell connections
"""
      if net_type[n_index] in ['pg_net', 'noinhib_net']:
        parameters[folder_num] = parameters[folder_num] + \
"""gc_connection_state=1
toggle_gc_connection() // turns off all gc cell connections
"""
      parameters[folder_num]=parameters[folder_num]+"""objref pwm
pwm=new PWManager()
pwm.hide(3) // close voltage window for faster run
do_everything()
quit()
// net_type %s
""" % net_type[n_index]
      columns.append("""
n = %s // n easier to type than num_of_cols
""" % (add_columns[c_index]))

for folder_num in range(num_of_jobs):
  folder="run_"+str(folder_num)
  fid=open(folder+"/parameters.hoc","w")
  fid.write(parameters[folder_num])
  fid.close()
  fid=open(folder+"/num_of_columns.hoc","w")
  fid.write(columns[folder_num])
  fid.close()

print"*** the following phases are external to this program:"
#************************************************************************
print "phase 3"
#************************************************************************
# phase 3 The parallel job is uploaded to NGS: the mod files are
# compiled and the jobs are started by each run worker executing
# build_net_Shep.hoc and loading the appropriate run_X/parameters.hoc
# and run_X/num_of_columns.hoc scripts. These run_workers automatically
# save the neuron tank data (and optionally spike time data used to
# create subsequent graphs).

print 'start parallel job by uploading this archive to NSG'
#
#
#************************************************************************
print "phase 4"
#************************************************************************
# In phase 4 the matlab tanks are created
#

#************************************************************************
print "phase 5"
#************************************************************************
# In phase 5 the tanks are analyzed to create and save matlab polar plot figures
#

#************************************************************************
print "phase 6"
#************************************************************************
# In phase 6 the raster plots of spike activity are created and saved.
#
