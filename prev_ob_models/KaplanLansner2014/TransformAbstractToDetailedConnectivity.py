import sys
import re
import numpy as np
import numpy.random as rnd
import random
import os
import pylab
import time
import matplotlib

class GetConnections(object):
    def __init__(self, params, comm=None, rank=0, debug=0):
        '''
        params: simulation parameter dictionary
        comm is the MPI communicator
        rank is the mpi rank of the process
        '''
        self.params = params
        self.comm = comm    # MPI.COMM_WORLD
        self.debug = debug
        if comm != None:
            self.n_proc = self.comm.Get_size()
        else:
            self.n_proc = 1
        self.my_rank = rank # process id
        rnd.seed(self.params['seed_connections'] + 1)


    def draw_connection(self, p, w, noise=0):
        """
        Decide whether a connection is drawn, given the possibility p.
        w : is the mean weight for the type of connection to be drawn.
        noise : is an optional argument and stands for the relative sigma of the normal distributino
        i.e. noise = 0.1 means sigma = 0.1 * w
        """
        if (p > rnd.random()):
            if (noise != 0 and noise > 0):
                weight = rnd.normal(w, noise)
                # check if sign of weight changed
                # if so, return 0
                if (np.sign(weight) != np.sign(w)):
                    return 0
                return weight
            elif (noise < 0): # stupid user, noise should be > 0
                print "WARNING, negative noise given to draw_connection(p, w, noise)!"
                noise *= (-1.0)
                weight = rnd.normal(w, noise)
                if (np.sign(weight) != np.sign(w)):
                    return 0
                return weight
            elif (noise == 0):
                return w
        return 0



    def take_log_weights(self, data):
        """
        if a weight is zero, it stays zero,
        otherwise: take the log
        """

        n_row = data[:, 0].size
        log_data = np.zeros(data.shape)
        for i in xrange(data.shape[0]):
            idx_nonzero = (data[i, :] > 0).nonzero()[0]
            log_data[i, idx_nonzero] = np.log(data[i, idx_nonzero])
        return log_data


    def get_mit_rsnp_connections(self, random_conn=False):
        '''
        This functions creates the OB->OC connections files from the matrix stored in self.params['ob_oc_abstract_weights_fn']
        '''
        print "Drawing MIT - RSNP connections .... "

        abstract_weights_non_negative = np.loadtxt(self.params['ob_oc_abstract_weights_fn'])
        abstract_weights = self.take_log_weights(abstract_weights_non_negative)
        if random_conn:
            rnd.seed(self.params['random_ob_oc_seed'])
            rnd.shuffle(abstract_weights)
            np.savetxt(self.params['ob_oc_abstract_weights_fn'].rsplit('.dat')[0] + '_random.dat', abstract_weights)
        assert (abstract_weights[:, 0].size == self.params['n_mit'])
        assert (abstract_weights[0, :].size == self.params['n_hc'] * self.params['n_mc'])
        # scale the abstract weights into the biophysical range
        w_max_abstract = abstract_weights.min()
        w_mit_rsnp_max = self.params['w_mit_rsnp_max']

        output = ""
        line_cnt = 0
        for tgt_mc in xrange(self.params['n_hc'] * self.params['n_mc']): #column
            for mit in xrange(self.params['n_mit']): # row
                w_in = abstract_weights[mit, tgt_mc]
                if (w_in < 0):
                    tgt_rsnps = self.get_rnd_targets(self.params['n_rsnp_per_mc'], self.params['n_tgt_rsnp_per_mc'])
#                    tgt_rsnps = random.sample(range(self.params['n_rsnp_per_mc']), int(round(self.params['n_tgt_rsnp_per_mc'])))
                    for tgt_rsnp in tgt_rsnps:
                        w_out = (w_in / w_max_abstract) * w_mit_rsnp_max 
                        w_noise = rnd.normal(w_out, w_out * self.params['w_mit_rsnp_sigma_frac'])
#                        w_noise = rnd.normal(w_out, self.params['w_mit_rsnp_max'] * self.params['w_mit_rsnp_sigma_frac'])
                        if (w_noise > self.params['weight_threshold']):
                            src_gid = self.params['mit_offset'] + mit
                            tgt_gid = self.params['rsnp_offset'] + tgt_rsnp + tgt_mc * self.params['n_rsnp_per_mc']
                            output += "%d\t%d\t%.8e\n" % (src_gid, tgt_gid, w_noise)
                            line_cnt += 1

        output_fn = self.params['conn_list_mit_rsnp']
        print 'Saving %d mit-rsnp connections to file %s' % (line_cnt, output_fn)
        first_line = "%d %d\n" % (line_cnt, 3)
        output_file = open(output_fn, 'w')
        output_file.write(first_line)
        output_file.write(output)
        output_file.close()


    def get_mit_pyr_connections(self, random_conn=False):
        '''
        This functions creates the OB->OC connections files from the matrix stored in self.params['ob_oc_abstract_weights_fn']
        '''
        print "Drawing MIT - PYR connections .... "
        abstract_weights_non_negative = np.loadtxt(self.params['ob_oc_abstract_weights_fn'])
        abstract_weights = self.take_log_weights(abstract_weights_non_negative)
        if random_conn:
            rnd.seed(self.params['random_ob_oc_seed'])
            rnd.shuffle(abstract_weights)
            np.savetxt(self.params['ob_oc_abstract_weights_fn'].rsplit('.dat')[0] + '_random.dat', abstract_weights)
        assert (abstract_weights[:, 0].size == self.params['n_mit'])
        assert (abstract_weights[0, :].size == self.params['n_hc'] * self.params['n_mc'])
        # scale the abstract weights into the biophysical range
        w_max_abstract = abstract_weights.max()
        w_min_abstract = abstract_weights.min()
        w_mit_pyr_max = self.params['w_mit_pyr_max']
        w_mit_pyr_matrix = np.zeros((self.params['n_mit'], self.params['n_hc'] * self.params['n_mc']))

        output = ""
        line_cnt = 0
        for tgt_mc in xrange(self.params['n_hc'] * self.params['n_mc']): #column
            for mit in xrange(self.params['n_mit']): # row
                w_in = abstract_weights[mit, tgt_mc]
                if (w_in > 0):
#                    for tgt_pyr in xrange(int(round(self.params['n_tgt_pyr_per_mc']))):
                    tgt_pyrs = self.get_rnd_targets(self.params['n_pyr_per_mc'], self.params['n_tgt_pyr_per_mc'])
                    for tgt_pyr in tgt_pyrs:
                        w_out = (w_in / w_max_abstract) * w_mit_pyr_max
                        w_noise = rnd.normal(w_out, w_out * self.params['w_mit_pyr_sigma_frac'])
                        if (w_noise > self.params['weight_threshold']):
                            src_gid = self.params['mit_offset'] + mit
                            tgt_gid = self.params['pyr_offset'] + tgt_pyr + tgt_mc * self.params['n_pyr_per_mc']
                            output += "%d\t%d\t%.8e\n" % (src_gid, tgt_gid, w_noise)
                            line_cnt += 1

                            w_mit_pyr_matrix[mit, tgt_mc] = w_noise
                            
        output_fn = self.params['conn_list_mit_pyr']
        print 'Saving %d mit-pyr connections to file %s' % (line_cnt, output_fn)
        first_line = "%d %d\n" % (line_cnt, 3)
        output_file = open(output_fn, 'w')
        output_file.write(first_line)
        output_file.write(output)
        output_file.close()

        np.savetxt(self.params['connection_matrix_detailed_ob_oc_dat'], w_mit_pyr_matrix)


    def get_oc_oc_connections(self, random_conn=False):
        """
        This is a serial version. To be implemented
        if comm.n_proc > 1: call the parallel version
        else: call the serial version
        """

        print "Drawing OC - OC connections .... "
        abstract_weights_non_negative = np.loadtxt(self.params['oc_oc_abstract_weights_fn'])
        abstract_weights = self.take_log_weights(abstract_weights_non_negative)
        if random_conn:
            rnd.shuffle(abstract_weights)
            rnd.seed(self.params['random_oc_oc_seed'])
            np.savetxt(self.params['oc_oc_abstract_weights_fn'].rsplit('.dat')[0] + '_random.dat', abstract_weights)

        assert (abstract_weights[:,0].size == self.params['n_hc'] * self.params['n_mc'])
        assert (abstract_weights[0,:].size == self.params['n_hc'] * self.params['n_mc'])
        w_max_abstract = abstract_weights.max()
        w_min_abstract = abstract_weights.min()

        w_pyr_pyr_global_max = self.params['w_pyr_pyr_global_max']
        w_pyr_rsnp_max = self.params['w_pyr_rsnp_max']
        output_pyr_pyr = ""
        line_cnt_pyr_pyr = 0
        output_pyr_rsnp = ""
        line_cnt_pyr_rsnp = 0
        cnt_discarded_conn = 0
        for src_mc in xrange(abstract_weights[:, 0].size):
            for tgt_mc in xrange(abstract_weights[:, 0].size):
                if (src_mc != tgt_mc):
                    w_in = abstract_weights[src_mc, tgt_mc]
                    if (w_in > 0): # draw several pyr -> pyr connections between the two MC
                        src_tgt_dict = {} # src_tgt_dict[src_gid] = [tgt_gid_0, ...] multiple connections between the same source and the same target are forbiddden
                        w_out = (w_in / w_max_abstract) * w_pyr_pyr_global_max
                        src_pyrs = rnd.randint(0, self.params['n_pyr_per_mc'], self.params['n_pyr_pyr_between_2mc'])
                        for src in np.unique(src_pyrs):
                            src_tgt_dict[src] = []
                        for src in src_pyrs:
                            src_pyr = src + src_mc * self.params['n_pyr_per_mc'] + self.params['pyr_offset']
                            tgt_pyr = rnd.randint(0, self.params['n_pyr_per_mc']) + tgt_mc * self.params['n_pyr_per_mc'] + self.params['pyr_offset']
                            src_tgt_dict[src].append(tgt_pyr)

                        # remove multiple instances of the same src-tgt connection
                        for src in src_pyrs:
                            n1 = len(src_tgt_dict[src])
                            src_tgt_dict[src] = np.unique(src_tgt_dict[src]).tolist()
                            cnt_discarded_conn += n1 - len(src_tgt_dict[src])
                            for tgt_pyr in src_tgt_dict[src]:
                                w_noise = self.draw_connection(1.0, w_out, noise=self.params['w_pyr_pyr_global_sigma'])
                                if (w_noise > self.params['weight_threshold']):
                                    output_pyr_pyr += "%d %d %.6e\n" % (src_pyr, tgt_pyr, w_noise)
                                    line_cnt_pyr_pyr += 1

                    elif (w_in < 0):
                        w_out = (w_in / w_min_abstract) * w_pyr_rsnp_max
                        src_pyrs = self.get_rnd_targets(self.params['n_pyr_per_mc'], self.params['n_pyr_rsnp_between_2mc']) # avoid double connections
                        for src in src_pyrs:
                            src_pyr = src + src_mc * self.params['n_pyr_per_mc'] + self.params['pyr_offset'] 
                            tgt_rsnp = rnd.randint(0, self.params['n_rsnp_per_mc']) + tgt_mc * self.params['n_rsnp_per_mc'] + self.params['rsnp_offset']
                            w_noise = self.draw_connection(1.0, w_out, noise=self.params['w_pyr_rsnp_sigma'])
                            if (w_noise > self.params['weight_threshold']):
                                output_pyr_rsnp += "%d %d %.6e\n" % (src_pyr, tgt_rsnp, w_noise)
                                line_cnt_pyr_rsnp += 1

        print 'Number of discarded pyr-pyr connections:', cnt_discarded_conn
        print 'Number of pyr-rsnp connections:', line_cnt_pyr_rsnp
        print 'Number of pyr-pyr connections:', line_cnt_pyr_pyr
        print 'Number of OC-OC connections:', line_cnt_pyr_pyr + line_cnt_pyr_rsnp
        output_fn_pyr_pyr = self.params['conn_list_pyr_pyr']
        output_file_pyr_pyr = open(output_fn_pyr_pyr, 'w')
        output_file_pyr_pyr.write("%d\t%d\n" % (line_cnt_pyr_pyr, 3))
        output_file_pyr_pyr.write(output_pyr_pyr)
        output_file_pyr_pyr.close()

        output_fn_pyr_rsnp = self.params['conn_list_pyr_rsnp']
        output_file_pyr_rsnp = open(output_fn_pyr_rsnp, 'w')
        output_file_pyr_rsnp.write("%d\t%d\n" % (line_cnt_pyr_rsnp, 3))
        output_file_pyr_rsnp.write(output_pyr_rsnp)
        output_file_pyr_rsnp.close()


    def get_pyr_readout_connections(self):
        '''
        '''
        print "Drawing OC - Readout connections .... "
        abstract_weights = np.loadtxt(self.params['oc_readout_abstract_weights_fn'])
        assert (abstract_weights[:, 0].size == self.params['n_hc'] * self.params['n_mc'])
        assert (abstract_weights[0, :].size == self.params['n_readout'])
        # scale the abstract weights into the biophysical range
        w_max_abstract = abstract_weights.max()
        w_min_abstract = abstract_weights.min()
        w_pyr_readout_max = self.params['w_pyr_readout']

        output = ""
        line_cnt = 0
        for tgt_cell in xrange(self.params['n_readout']):
            for src_mc in xrange(self.params['n_hc'] * self.params['n_mc']): # row
                w_in = abstract_weights[src_mc, tgt_cell]
                if (w_in > 0):
                    w_out = (w_in / w_max_abstract) * w_pyr_readout_max
                elif (w_in < 0):
                    w_out = (-1.0) * (w_in / w_min_abstract) * w_pyr_readout_max
                if (abs(w_in > self.params['weight_threshold'])):
                    for src_pyr in xrange(self.params['n_pyr_per_mc']):
                        src_gid = self.params['pyr_offset'] + src_mc * self.params['n_pyr_per_mc'] + src_pyr
                        tgt_gid = self.params['readout_offset'] + tgt_cell
                        output += "%d \t%d\t%.8e\n" % (src_gid, tgt_gid, w_out)
                        line_cnt += 1
        first_line = "%d %d\n" % (line_cnt, 3)
        output_file = open(self.params['conn_list_pyr_readout'], 'w')
        output_file.write(first_line)
        output_file.write(output)
        output_file.close()




    def get_matrix_from_conn_list(self, conn_list_fn, src_type, tgt_type):

        if ((src_type == 'pyr') or (src_type == 'rsnp')):
            n_src = self.params['n_mc'] * self.params['n_hc']
        else:
            n_src = self.params['n_%s' % src_type]
        if ((tgt_type == 'pyr') or (tgt_type == 'rsnp')):
            n_tgt = self.params['n_mc'] * self.params['n_hc']
        else:
            n_tgt = self.params['n_%s' % tgt_type]
        src_offset = self.params['%s_offset' % src_type]
        tgt_offset = self.params['%s_offset' % tgt_type]

        m = np.zeros((n_src, n_tgt))
        d = np.loadtxt(conn_list_fn, skiprows=1)

        if ((tgt_type == 'pyr') or (tgt_type == 'rsnp')):
            n_tgt_per_mc = self.params['n_%s_per_mc' % tgt_type]
        else:
            n_tgt_per_mc = 1
        if ((src_type == 'pyr') or (src_type == 'rsnp')):
            n_src_per_mc = self.params['n_%s_per_mc' % src_type]
        else:
            n_src_per_mc = 1

        for i in xrange(d[:, 0].size):
            src = (d[i, 0] - src_offset) / n_src_per_mc
            tgt = (d[i, 1] - tgt_offset) / n_tgt_per_mc
            m[src, tgt] += d[i, 2]

        return m


    def get_rnd_targets(self, gid_max, n):
        """
        gid_max: upper boundary
        n: number of random integers to be drawn within range [0, gid_max]
        """
        assert (n < gid_max), 'ERROR: Can\'t provide more unique numbers in range (0, gid_max) than gid_max. n is too high!'
        tgts = rnd.randint(0, gid_max, n)
        tgts = np.unique(tgts).tolist()
        while len(tgts) < n:
            # check if rnd_int is already in l 
            tgts.append(rnd.randint(0, gid_max))
            tgts = np.unique(tgts).tolist()
        return tgts

