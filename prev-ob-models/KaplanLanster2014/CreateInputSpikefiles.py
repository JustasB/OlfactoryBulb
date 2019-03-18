import os 
import sys
import time
import numpy as np
import MergeSpikefiles
import utils



class CreateInputSpikefiles(object):


    def __init__(self, params, pn_max=None, comm=None):

        self.params = params
        if pn_max == None:
            self.pn_max = self.params['n_patterns']
        else:
            self.pn_max = pn_max
    
        self.Merger = MergeSpikefiles.MergeSpikefiles(params)

        self.comm = comm
        if self.comm == None:
            self.pc_id = 0
            self.n_proc = 1
        else:
            self.pc_id = comm.Get_rank()
            self.n_proc = comm.Get_size()
        self.my_gids = {}
        self.conn_dict_tgt_key = {}


    def set_my_gids(self, cell_type):
        """
        If run on multiple cores, distribute the GIDs among the processes
        """
        gid_min, gid_max = utils.distribute_n(self.params['n_%s' % cell_type], self.n_proc, self.pc_id)
        self.my_gids[cell_type] = (gid_min + self.params['%s_offset' % cell_type], gid_max + self.params['%s_offset' % cell_type])


    def get_spike_times(self, pn, cell_type='mit'):
        """
        This function creates the spike times - output files:
            # 1) tgt_gid    weight     # each connection gets a seperate netcon with its own weight
            # 2) time of spike     tgt_gid     netcon_index
        """

        print "Collecting spiketimes for celltype: %s pattern_nr: %d" % (cell_type, pn)
        fn = self.params['%s_spiketimes_merged_fn_base' % (cell_type)] + str(pn) + '.dat'
        if (os.path.exists(fn) == False):
            print "Merging data for file: ", fn
            self.Merger.merge_spiketimes_files(self.params['%s_spiketimes_fn_base' % (cell_type)], self.params['%s_spiketimes_merged_fn_base' % (cell_type)], pn)
        print "Loading data from:", fn
        d = np.loadtxt(fn)

        info_txt = "\n\tNo spikes were found!\n\tDon't forget that you might need to copy the mit_spiketimes_merged_* files\n \
        from another folder to this folder %s" % (self.params['spiketimes_folder'])
        assert (d.size > 0), info_txt

        spike_times_dict = {}
        for src_gid in xrange(self.params['%s_offset' % cell_type], self.params['%s_offset' % cell_type] + self.params['n_%s' % cell_type]):
            idx = (d[:, 1] == src_gid).nonzero()[0]
            spike_times = d[idx, 0]
            spike_times_dict[src_gid] = spike_times

        return spike_times_dict


    def get_conn_dict(self, tgt_cell_type):
        """
        This function creates a dictionary storing the adjecency list with the target GID as key
        and writes the information connecting the target_gid, the weights and netcons into a file to be
        read by setup_network.hoc
        """

        fn = self.params['conn_list_mit_%s' % tgt_cell_type]
        d = np.loadtxt(fn, skiprows=1)

        output_fn = self.params['mit_%s_tgt_netcon_weight_fn' % tgt_cell_type] + '%d.dat' % self.pc_id
        f = file(output_fn, 'w')
        info_to_write = ''

        self.conn_dict_tgt_key[tgt_cell_type] = {}
        n_lines = 0
#        for tgt_gid in xrange(self.params['%s_offset' % tgt_cell_type], self.params['%s_offset' % tgt_cell_type] + self.params['n_%s' % tgt_cell_type]):
        for tgt_gid in xrange(self.my_gids[tgt_cell_type][0], self.my_gids[tgt_cell_type][1]):
#            print 'Getting sources for %s %d / %d (%.2f percent)' % (tgt_cell_type, tgt_gid - self.params['%s_offset' % tgt_cell_type], self.params['n_%s' % tgt_cell_type], \
#                    100. * (float(tgt_gid - self.params['%s_offset' % tgt_cell_type])/ self.params['n_%s' % tgt_cell_type]))
            idx = (d[:, 1] == tgt_gid).nonzero()[0]
            srcs, weights = d[idx, 0], d[idx, 2]
            conn_info_list = [(srcs[i], weights[i]) for i in xrange(idx.size)]
            

            self.conn_dict_tgt_key[tgt_cell_type][tgt_gid] = conn_info_list
#            if tgt_gid == 0:
#                print 'DEBUG'
#                print 'list', conn_info_list
#            print 'dict', tgt_gid, self.conn_dict_tgt_key[tgt_cell_type][tgt_gid]
            # for each tgt cell write an individual file storing the netcon id and the weight for this connection
            # this information will be retrieved from a function in setup_network.hoc which adds synapses to the tgt cell's soma accordingly
            netcon_id_offset = 0 # could be different, depending on the consistency of the neuron_model template 
            # if there are default netcon objects in the netconlist and if setup_network inserts spikes into this netconlist
            # the netcon_id_offset should be != 0 --> but in general it would be smart if these 'offline spikes'
            # are inserted through a dedicated netconlist into the neuron --> netcon_id_offset = 0
            for i_, src_w in enumerate(conn_info_list):
                info_to_write += '%d\t%d\t%.6e\n' % (tgt_gid, i_ + netcon_id_offset, src_w[1])
                n_lines += 1

        n_columns = 3
        f.write('%d %d\n' % (n_lines, n_columns))
        print 'Writing to:', output_fn
        f.write(info_to_write)
        f.flush()


    def write_spiketimes_for_netcons(self, pn, spike_times_dict, tgt_cell_type='pyr'):
        """
        Arguments:
        """


        output_fn = self.params['mit_%s_spiketimes_tgt_netcon_fn' % tgt_cell_type] + '%d_%d.dat' % (pn, self.pc_id)
        f = file(output_fn, 'w')
        info_to_write = ''

        n_lines = 0
        for tgt_gid in xrange(self.my_gids[tgt_cell_type][0], self.my_gids[tgt_cell_type][1]):

#            print 'debug', tgt_cell_type, tgt_gid, self.conn_dict_tgt_key[tgt_cell_type]
            for i_netcon, src_w in enumerate(self.conn_dict_tgt_key[tgt_cell_type][tgt_gid]):
                src_gid, w = src_w[0], src_w[1]
                # get spike times emitted by src_gid
                spiketimes = spike_times_dict[src_gid]
                for i_spike, spike in enumerate(spiketimes):
                    info_to_write += '%d\t%.2f\t%d\n' % (tgt_gid, spike, i_netcon)
                    n_lines += 1

        n_columns = 3
        f.write('%d %d\n' % (n_lines, n_columns))
        print 'Writing to:', output_fn
        f.write(info_to_write)
        f.flush()


    def merge_files(self, pn, cell_type):

        fn_0 = self.params['mit_%s_tgt_netcon_weight_fn' % cell_type]
        tmp_fn_0 = 'delme_%d' % np.random.randint(10**6)
        cat_cmd = 'cat %s* > %s' % (fn_0, tmp_fn_0)
        mv_cmd = 'mv %s %s' % (tmp_fn_0, self.params['mit_%s_tgt_netcon_weight_fn_merged' % cell_type])
        os.system(cat_cmd)
        os.system(mv_cmd)

        if self.pc_id == 0:
            fn_0 = self.params['mit_%s_spiketimes_tgt_netcon_fn' % cell_type] + str(pn)
            tmp_fn_0 = 'delme_%d' % np.random.randint(10**6)
            cat_cmd = 'cat %s* > %s' % (fn_0, tmp_fn_0)
            mv_cmd = 'mv %s %s' % (tmp_fn_0, self.params['mit_%s_spiketimes_tgt_netcon_fn_merged' % cell_type] + '%d.dat' % pn)
            os.system(cat_cmd)
            os.system(mv_cmd)





if __name__ == '__main__':
    
    if len(sys.argv) > 1:
        param_fn = sys.argv[1]
        if os.path.isdir(param_fn):
            param_fn += '/Parameters/simulation_parameters.json'
        import json
        f = file(param_fn, 'r')
        print 'Loading parameters from', param_fn
        params = json.load(f)

    else:
        import simulation_parameters
        param_tool = simulation_parameters.parameter_storage()
        params = param_tool.params

    try:
        if sys.argv[1].isdigit():
            pn_max = sys.argv[1]
        else:
            pn_max = int(sys.argv[2])
    except:
        print 'Processing all patterns'
        pn_max = params['n_patterns']

#    try:
#        from mpi4py import MPI
#        comm = MPI.COMM_WORLD
#        rank = comm.Get_rank()
#    except:
#        comm, rank = None, 0

    CIS = CreateInputSpikefiles(params, pn_max=pn_max)
    CIS.set_my_gids(cell_type='pyr')
    CIS.set_my_gids(cell_type='rsnp')
    CIS.get_conn_dict(tgt_cell_type = 'pyr') # write the netconid info
    CIS.get_conn_dict(tgt_cell_type = 'rsnp') # write the netconid info
    for pn in xrange(pn_max):
        spike_times_dict = CIS.get_spike_times(pn, cell_type='mit')
        CIS.write_spiketimes_for_netcons(pn, spike_times_dict, tgt_cell_type='pyr')
        CIS.write_spiketimes_for_netcons(pn, spike_times_dict, tgt_cell_type='rsnp')

        CIS.merge_files(pn, 'pyr')
        CIS.merge_files(pn, 'rsnp')
