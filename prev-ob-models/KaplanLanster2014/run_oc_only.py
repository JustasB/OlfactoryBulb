import os
import sys
import time
import numpy as np
import simulation_parameters # defines simulation parameters

param_tool = simulation_parameters.parameter_storage()
params = param_tool.params
param_tool.write_parameters_to_file(params["info_file"]) # human readable
param_tool.write_parameters_to_file() # 
param_tool.hoc_export() # write the simulation parameters to a NEURON executable file
param_tool.print_cell_gids()

t_all_0 = time.time()

def add_first_line(pn):
    """"
    This function is required if OC only is run 
    Adds a nrow, ncol first line in the mit_spiketimes_merged_ files 
    """
    fn = params['mit_spiketimes_merged_fn_base'] + '%d.dat' % pn
    fn_out = params['mit_spiketimes_merged_fn_base'] + 'for_neuron_%d.dat' % pn
    print 'Loading', fn
    assert (os.path.exists(fn)), 'File not found: %s\n\tHave you merged the files before?' % fn
    d = np.loadtxt(fn)
    n_row, n_col = d.shape[0], d.shape[1]
    # read the current contents of the file
    f = open(fn)
    text = f.read()
    f.close()
    # open the file again for writing
    f = open(fn_out, 'w')
    f.write('%d %d\n' % (n_row, n_col))
    # write the original contents
    f.write(text)
    f.close()


if __name__ == '__main__':
    for pn in xrange(params['n_patterns']):
        os.system("rm %s*" % (params["pyr_spiketimes_fn_base"]))
        os.system("rm %s*" % (params["pyr_spike_fn_base"]))
        os.system("rm %s*" % (params["pyr_volt_fn_base"]))

#    param_tool.set_gids_to_record([6806])
#    for pn in [0]:
    for pn in xrange(params['n_patterns']):
        t1 = time.time()
        # add a nrow, ncol first line in the mit_spiketimes_merged_ files 
        add_first_line(pn)

    #    neuron_command = "mpirun -np %d $(which nrniv) -mpi -nobanner -nogui \
    #            -c \"x=%d\" -c \"strdef param_file\" -c \"sprint(param_file, \\\"%s\\\")\" start_file_oc_only_recognition_task.hoc > delme%d" \
    #            % (params['n_proc'], pn, params['hoc_file'], pn)

        neuron_command = "mpirun -np %d $(which nrniv) -mpi -nobanner -nogui \
                -c \"x=%d\" -c \"strdef param_file\" -c \"sprint(param_file, \\\"%s\\\")\" start_file_oc_only_recognition_task.hoc" \
                % (params['n_proc'], pn, params['hoc_file'])
    #    neuron_command = "$(which nrniv) -mpi -nobanner -nogui \
    #            -c \"x=%d\" -c \"strdef param_file\" -c \"sprint(param_file, \\\"%s\\\")\" start_file_oc_only_recognition_task.hoc" \
    #            % (pn, params['hoc_file'])

        os.chdir('neuron_files') # this is important to avoide problems with tabchannel files and the functions defined therein
        print 'NEURON simulation starts for pattern %d ...' % pn
        os.system(neuron_command)

        t2 = time.time() - t1
        print "Simulating %d cells for %d ms took %.3f seconds or %.2f minutes" % (params['n_cells_oc'] , params["t_sim"], t2, t2/60.)

        os.chdir('../') # this is important to avoide problems with tabchannel files and the functions defined therein

    t_all_1 = time.time() - t_all_0
    print "Simulating %d patterns %d cells for %d ms took %.3f seconds or %.2f minutes" % (params['n_patterns'], params['n_cells_oc'], params["t_sim"], t_all_1, t_all_1/60.)
