#!/usr/bin/python

"""batch_runs.py will create multiple folders for simple_circuit model
runs with names like run_X where X goes from 1 to the number of jobs.

phase 1 the run_X folders are created and the programs and folders by
recursively copying simple_circuit into them.

phase 2 the num_of_columns.hoc and parameters.hoc files are created in
each run_x folder so that they are all ready to run with their
different assigned values.  A parameters[] list prepares the contents of
each of these files before they are written.

For now the following phases are external to this program:
initiated with
qsub run.pbs
and 
matlab
batch_polar_plots
or
batch_bi_polar_plots

phase 3 the mod files are compiled in each run_x folder and the jobs
are started by running build_net_Shep.hoc in each run_X folder These
jobs automatically save the tank and spike time data used to create
subsequent graphs

phase 4 the matlab tanks are created

phase 5 the tanks are analyzed to create and save matlab polar plot
figures

phase 6 the raster plots of spike activity are created and
saved.
"""


#************************************************************************
print "phase 1"
#************************************************************************
# In phase 1 the run_X folders are created and the programs and folders by recursively copying simple_circuit
# into them.

# the task of create arrays is to create "both" which has both the breathing peak rate
# and the stimulus (light) peak rate in list of (B, S) tuples.  Once this is created
# we can determine the number of jobs by multiplying the length of both by three for
# three different types of network models
execfile("create_arrays.py")
jobs_per_nn=len(both)
num_of_columns_tested=len(add_columns) # add_columns supplied by create_arrays.py
num_of_nn_types=4 # pg, gc, full, noinhib
num_of_jobs=jobs_per_nn*num_of_nn_types*num_of_columns_tested
print "num_of_jobs: "+ str(num_of_jobs)
# absolute_path="/home/tmm46/projects/VerhagenLab/20150611/batch_runs/"
#absolute_path="/home/tmm46/projects/VerhagenLab/20150622/batch_runs/"
#absolute_path="/home/tmm46/projects/VerhagenLab/20150708/batch_runs/"
absolute_path="/home/tmm46/projects/VerhagenLab/20150720/cols_50_400/batch_runs/"

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

# now copy over the simple_circuit folder
# http://stackoverflow.com/questions/12683834/how-to-copy-directory-recursively-in-python-and-overwrite-all
for folder_num in range(num_of_jobs):
  cmd_string="cp -rf ../simple_circuits/* run_"+str(folder_num)
  os.system(cmd_string)

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

num_of_sims = num_of_jobs
num_of_procs = 14 * 4 # 14 nodes at 4 processoers per node
serial_num_of_sims=int(num_of_sims / num_of_procs + 1) # ceiling of division because some procs will 
# be running 1 extra unless exactly an integer result of division
one_sim_time =  (2*60+35) / 4. # it took 2 hours 35 mins to previously run 4 serial jobs
time_to_run = serial_num_of_sims * one_sim_time
# where the 4 serial jobs came from 192 jobs/56 processors = 3.4 jobs/proc
print "this job is predicted to take "+repr(time_to_run)+" mins = "+repr(time_to_run/60)+" hrs"
print "on "+str(num_of_procs)+" processors"

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

for folder_num in range(num_of_columns_tested): #, 2*jobs_per_nn):  # pg mediated inhibition with nc[15][]=1
  parameters[folder_num]="""
breathing_period=400
light_period=399 // 300 is a short run, 398 regular
breath_peak_rate = %d
light1_peak_rate = %d
light2_peak_rate = 0
breath_half_width=30
light_half_width=30
for i=0, n-1 {
  nc[15][i].weight = 1 // turn off breathing input to pg cells
}
toggle_gc_connection() // turns off all gc cell connections
objref pwm
pwm=new PWManager()
pwm.hide(3) // close voltage window for faster run
do_everything()
quit()
// net_type pg_net
""" % both[folder_num%jobs_per_nn]
  columns[folder_num]="""
n = %s // n easier to type than num_of_cols
""" % tuple(add_columns)[folder_num % num_of_columns_tested]

for folder_num in range( num_of_columns_tested, 2* num_of_columns_tested): # gc mediated
  parameters[folder_num]="""
breathing_period=400
light_period=399 // 300 is a short run, 398 regular
breath_peak_rate = %d
light1_peak_rate = %d
light2_peak_rate = 0
breath_half_width=30
light_half_width=30
// for some reason toggle_pg_connection() was causing an error however
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
objref pwm
pwm=new PWManager()
pwm.hide(3) // close voltage window for faster run
do_everything()
quit()
// net_type gc_net
""" % both[folder_num%jobs_per_nn]
  columns[folder_num]="""
n = %s // n easier to type than num_of_cols
""" % tuple(add_columns)[folder_num % num_of_columns_tested]

for folder_num in range(2* num_of_columns_tested, 3* num_of_columns_tested): # gc mediated
  parameters[folder_num]="""
breathing_period=400
light_period=399 // 300 is a short run, 398 regular
breath_peak_rate = %d
light1_peak_rate = %d
light2_peak_rate = 0
breath_half_width=30
light_half_width=30
// for some reason toggle_pg_connection() was causing an error however
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
toggle_gc_connection() // turns off all gc cell connections
objref pwm
pwm=new PWManager()
pwm.hide(3) // close voltage window for faster run
do_everything()
quit()
// net_type noinhib_net
""" % both[folder_num%jobs_per_nn]
  columns[folder_num]="""
n = %s // n easier to type than num_of_cols
""" % tuple(add_columns)[folder_num % num_of_columns_tested]

for folder_num in range(3* num_of_columns_tested, 4* num_of_columns_tested): # gc mediated
  parameters[folder_num]="""
breathing_period=400
light_period=399 // 300 is a short run, 398 regular
breath_peak_rate = %d
light1_peak_rate = %d
light2_peak_rate = 0
breath_half_width=30
light_half_width=30
objref pwm
pwm=new PWManager()
pwm.hide(3) // close voltage window for faster run
do_everything()
quit()
// net_type full_net
""" % both[folder_num%jobs_per_nn]
  columns[folder_num]="""
n = %s // n easier to type than num_of_cols
""" % tuple(add_columns)[folder_num % num_of_columns_tested]

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
# In phase 3 the mod files are compiled in each run_x folder and the
# jobs are started by running build_net_Shep.hoc in each run_X folder
# These jobs automatically save the tank and spike time data used to
# create subsequent graphs
# 
# use simpleque to create a list of nrnivmodl tasks and run them
fid=open(absolute_path+"tasklist","w")
for folder_num in range(num_of_jobs):
  folder=absolute_path+"run_"+str(folder_num)
  fid.write("cd %s; source /home/tmm46/.bash_profile; /home/tmm46/bin/neuron/nrn/x86_64/bin/nrnivmodl; /home/tmm46/bin/neuron/nrn/x86_64/bin/nrngui %s/build_net_Shep.hoc\n" % (folder,folder))
fid.close()

print 'start job by running "qsub run.pbs"'
# 
# Following http://maguro.cs.yale.edu/mediawiki/index.php/SimpleQueue
# I used a command 
# /usr/local/cluster/software/installation/SimpleQueue/sqPBS.py gen 8 tmm46 nrn_task tasklist > run.pbs
# to generate run.pbs however then I edited run.pbs to 1 node and 10 ppn to
# more efficiently use louise (it seems to work for our first case of requiring 10 jobs).
# 
# To run first batch_run.py is run which will create some run_X folders and tasklist.  Then tasklist is run with
# qsub run.pbs
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
