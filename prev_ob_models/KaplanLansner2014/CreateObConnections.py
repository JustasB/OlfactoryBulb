import numpy
import numpy.random as rnd
import random
import os
import sys

class CreateObConnections(object):
    """
    CreateObConnections creates a connections of the following types:
         ---< means excitatory; ---o means inhibitory
     * ORN ---<  MT
     * ORN ---<  PG
     * MT  ---<  PG: (reciprocal, i.e. the same MT cell that excites a PG cell gets inhibited by that one)
     * PG  ---o  MT: (reciprocal)
     * PG  ---o  MT: (serial (one way)
     * MT  ---<  GRAN: local (within one glomerulus)
     * MT  ---<  GRAN: global (between glomeruli)
     * GRAN---o  MT: local (within one glomerulus)
     * GRAN---o  MT: global (between glomeruli)

    The output files are in the NEURON friendly format:
    n 3  # n = number of connections, 3 stands for number of columns
    src_id  tgt_id  weight

    Example:
    2 3
    1   5   0.3
    2   4   -0.5

    The first line includes the number of rows n and the number of columns 
    This is necessary to make the file easily readable by hoc via matrix.scanf(read_data)
    This means the cell with the global identifier 2 (gid, which must be unique)
    connects to target cell with gid 5 with a weight of 0.3, i.e. excitatory
    Cell ids start with 1.

    """

    def __init__(self, param_dict):
        self.params = param_dict
        self.folder_name = self.params['folder_name']
        self.orn_params_fn_base = self.params['orn_params_fn_base']
        self.mit_params_fn_base = self.params['mit_params_fn_base']
        self.n_or = self.params['n_or']
        self.n_gor = self.params['n_gor']
        self.n_orn_x = self.params['n_orn_x'] # number of ORNs expressing the same OR
        self.n_patterns = self.params['n_patterns']
        self.n_orn_y = self.params['n_orn_y'] # number of different OR
        self.n_orn = self.params['n_orn']
        self.n_mit_x = self.params['n_mit_x'] # number of MT cells per glomerulus
        self.n_mit_y = self.params['n_mit_y'] # number of glomeruli
        self.n_mit = self.params['n_mit']
        self.global_offset = self.params['global_offset']
        self.n_proc = self.params['n_proc']
#        self.orn_mit_exc = {}  # { tgt_gid : numpy.array((src_gid, tgt_gid, weight)) }
#        self.orn_mit_inh = {}  # same for inhibitory connections
        random.seed(self.params['seed'])
        rnd.seed(self.params['seed'])

    def create_orn_mit_weight_matrix(self):
        """
        n_orn_x ORNs connect to n_mit_x.
        The ORNs are ordered by their sensitivity (number of expressed receptors, controlled by their gor value).
        Example for 2 target MT cells:

        ORNs: 1 2 3 4 5 6 7 8 
               \||/    \||/
        MT:      9      10 
        ORNs 1,2,3,4 excite MT cell 9
        ORNs 5,6,7,8 excite MT cell 10
        """

        assert (self.n_mit_x < self.n_orn_x), "ERROR: Having more mitral cells than ORNs does not make sense\n:\tn_mit_x=%d > n_orn_x=%d" % (self.n_mit_x, self.n_orn_x)

        # Each mitral cell gets excitatory input from a set of ORNs
        n_orn_mit = self.params['rel_orn_mit']
        self.conn_mat_orn_mit = [[] for i in xrange(self.n_mit)]

        n_mit_x = self.params['n_mit_x']
        n_mit_y = self.params['n_mit_y']
        n_orn_x = self.params['n_orn_x']
        mit_offset = self.params['mit_offset']
        orn_offset = self.params['orn_offset']
        w_orn_mit = self.params['w_orn_mit_target'] / n_orn_mit
        w_orn_mit_sigma = self.params['w_orn_mit_sigma'] * w_orn_mit # convert percentage -> conductance
        w_orn_mult = self.params['w_orn_mit_mult']
        for glom in xrange(n_mit_y):
            for tgt in xrange(n_mit_x):
                tgt_gid = glom * n_mit_x + tgt + mit_offset
                src_offset = tgt * n_orn_mit + glom * n_orn_x + orn_offset # MT cells receive exc input from discjunct sets of ORNs
                for src in xrange(n_orn_mit):
                    src_gid = src + src_offset
                    orn_id_x = src_gid % n_orn_x
                    # as the output rate of ORNs projecting to the same glom is different, the weight is modified to compensate for this
                    w_mult = w_orn_mult + ((1. -  w_orn_mult) / (n_orn_x - 1)) * orn_id_x
                    w = rnd.normal(w_orn_mit * w_mult, w_orn_mit_sigma)
                    if (w > 0.0):
                        self.conn_mat_orn_mit[tgt_gid - mit_offset].append((src_gid, w))


    def connect_orn_mit(self):
        """
        n_orn_x ORNs connect to n_mit_x.
        The ORNs are ordered by their sensitivity (number of expressed receptors, controlled by their gor value).
        Example for 2 target MT cells:

        ORNs: 1 2 3 4 5 6 7 8 
               \||/    \||/
        MT:      9      10 
        ORNs 1,2,3,4 excite MT cell 9
        ORNs 5,6,7,8 excite MT cell 10
        """

        self.create_orn_mit_weight_matrix()

        # increase weights for the second MT cells in each row
        mit_offset = self.params['mit_offset']
        n_mit_y = self.params['n_mit_y']
        n_mit_x = self.params['n_mit_x']

        # change the ORN -< MT projections for some certain MT cells
        for j in xrange(len(self.params['orn_mit_change_ids'])):
            col = self.params['orn_mit_change_ids'][j]
            tgt_mts_to_increase = [mit_offset + col + i * n_mit_x for i in xrange(n_mit_y)]
            for i in tgt_mts_to_increase:
                self.change_orn_mit_weights(i, self.params['orn_mit_change_factors'][j])

        assert (self.n_mit_x < self.n_orn_x), "ERROR: Having more mitral cells than ORNs does not make sense\n:\tn_mit_x=%d > n_orn_x=%d" % (self.n_mit_x, self.n_orn_x)

        orn_offset = self.params['orn_offset']
        lines_to_write = ""
        conn_cnt = 0
        n_orn = self.params['n_orn']
        n_mit = self.params['n_mit']
        for tgt in xrange(n_mit):
            for src in xrange(len(self.conn_mat_orn_mit[tgt])):
                (src_gid, w) = self.conn_mat_orn_mit[tgt][src]
                if (w > 0.0):
                    conn_cnt += 1
                    tgt_gid = tgt + mit_offset
                    lines_to_write += "%d\t%d\t%.6E\n" % (src_gid, tgt_gid, w)

        # save connections in plain text format
        f_plain = open(self.params['conn_list_orn_mit'], 'w')
        first_line = "%d 3\n" % (conn_cnt)
        f_plain.write(first_line)
        f_plain.write(lines_to_write)
        f_plain.close()


    def create_orn_pg_weight_matrix(self):
        """
        ORNs: 5 6 7 8
               \||/ 
        PG:      a  
                 |  
                 o  
        MT:      9  
        """
        print "Creating connections: ORN ---< PG"

        # there are self.params['rel_pg_mit'] times more periglomerular cells than MT cells
        # and according to Shepherds "Synaptic organization of the brain" approx. 25 % of the PG-MT dendro-dendritic connections are reciprocal
        # i.e. these are treated here as feed-forward inhibitory connections from PG ---o MT
        # There are two populations of PG cells: those receiving direct input from ORNs (TH-ir PG cells according to Toida 2008 p. 210) and those PG cells (CB-ir) lying a bit deeper in the glomerular level which do not receive ORN input
        # The TH-ir population can be split up into to different populations again: one receiving ORN input and providing feedforward inhibition (n_pg_ff_inh),
        # the other population (=n_pg_wta_inh) receives additionally to the ORN input (from different ORNs here) also input from a closeby mitral cell and inhibits _other_ mitral cells via serial synapses  (ref Toida 2008)
        n_pg_ff_inh = int(round(self.params['n_pg_x_serial'] * (1. - self.params['rel_ff_wta_inhibition'])))
        n_pg_wta_inh = int(round(self.params['n_pg_x_serial'] * self.params['rel_ff_wta_inhibition']))
        n_pg_x = self.params['n_pg_x']
        n_pg_y = self.params['n_pg_y']
        n_orn_x = self.params['n_orn_x']
        w_orn_pg = self.params['w_orn_pg_target'] / (self.params['rel_orn_mit'])
        w_orn_pg_sigma = self.params['w_orn_pg_sigma'] * w_orn_pg # convert percentage -> conductance
        w_orn_mult = self.params['w_orn_pg_mult']
        n_orn = self.params['n_orn']
        n_pg = self.params['n_pg']
#        self.conn_mat_orn_pg = numpy.zeros((n_orn, n_pg))
        self.conn_mat_orn_pg = [[] for i in xrange(n_pg)]
        orn_offset = self.params['orn_offset']
        pg_offset = self.params['pg_offset']

        
        for glom in xrange(n_pg_y):
            src_offset = glom * n_orn_x
            tgt_offset = glom * n_pg_x + self.params['pg_offset']
            # 1) Those receiving FF exc from ORNs with lower sensitivity. They get no input from MIT cells
            for tgt in xrange(n_pg_ff_inh):
                # for each of these tgt pg cells decide from which ORN it receives excitation
                # divide the 'serial' population into n_mit_x groups, each group receives excitation from a different set of ORNs (still all ORNs express the same OR)
                group = tgt % self.params['n_mit_x']
                # define the range of ORN ids which excite the given group
                # as the target PG cells shall receive input from less sensitive ORNs than the MT cell to be inhibited by this group of 'serial' PG cells,
                # the set of source ORNs need to be shifted to the left in the ORN array (where the ORNs with smaller gor are)
                id_shift = self.params['orn_inh_shift'] * self.params['rel_orn_mit']
                src_orns = range(int(round(group * self.params['rel_orn_mit'] - id_shift)), int(round((group + 1) * self.params['rel_orn_mit'] - id_shift)))
                # this could be something to try: let all the less sensitive orns inhibit the target MT cell - not only a subgroup
                # this has as a consequence that MT cells coding for lower concentrations would get effectively more inhibitory input (in terms of numbers) than MT cells coding for high concentrations
                # because MT cells coding for high concetrations get input from very insensitive ORNs which have a low gor value
                # Thus, this might lead to different behavior and a different more complex way for scaling the weights might be adequate.
                for orn in src_orns:
                    if (orn >= 0):
                        src_gid = orn + src_offset
                        orn_id_x = src_gid % n_orn_x
                        w_mult = w_orn_mult + ((1. -  w_orn_mult) / (n_orn_x - 1)) * orn_id_x
                        tgt_gid = tgt + tgt_offset
                        w = rnd.normal(w_orn_pg * w_mult, w_orn_pg_sigma)
                        if (w > 0.0):
                            self.conn_mat_orn_pg[tgt_gid - pg_offset].append((src_gid, w))


            # 1) Those receiving FF exc from ORNs with the same sensitivity as the ones that excite the closeby MIT cell. They get additional exc input from this one MIT cell, but these PG cells inhibit other MIT cells.
            # purpose: wta mechanism
            for tgt in xrange(n_pg_ff_inh, n_pg_ff_inh + n_pg_wta_inh):
                # distribute these pg cells evenly among the mit cells
                group = tgt % self.params['n_mit_x']
                src_orns = range(int(round(group * self.params['rel_orn_mit'])), int(round((group + 1) * self.params['rel_orn_mit'])))
                for orn in src_orns:
                    if (orn >= 0):
                        src_gid = orn + src_offset
                        orn_id_x = src_gid % n_orn_x
                        w_mult = w_orn_mult + ((1. -  w_orn_mult) / (n_orn_x - 1)) * orn_id_x
                        tgt_gid = tgt + tgt_offset
                        w = rnd.normal(w_orn_pg * w_mult, w_orn_pg_sigma)
                        if (w > 0.0):
                            self.conn_mat_orn_pg[tgt_gid - pg_offset].append((src_gid, w))




    def connect_orn_pg(self):
        """
        ORNs: 5 6 7 8
               \||/ 
        PG:      a  
                 |  
                 o  
        MT:      9  
        """

        self.create_orn_pg_weight_matrix()
#        for i in self.params['pg_change_target_groups']:
#            self.change_orn_pg_weights(i, self.params['orn_pg_factors'][i])

        n_orn = self.params['n_orn']
        n_pg = self.params['n_pg']
        n_pg_x = self.params['n_pg_x']
        n_pg_y = self.params['n_pg_y']
        n_orn_x = self.params['n_orn_x']
        pg_offset = self.params['pg_offset']

        # certain PG cells may get modified excitatory weights from ORNs to extra-compensate for ORN output rates
        for j in xrange(len(self.params['orn_pg_change_ids'])):
            col = self.params['orn_pg_change_ids'][j]
            tgt_pgs_to_increase = [pg_offset + col + i * n_pg_x for i in xrange(n_pg_y)]
            for i in tgt_pgs_to_increase:
                self.change_orn_pg_weights(i, self.params['orn_pg_change_factors'][j])

        orn_offset = self.params['orn_offset']
        conn_cnt = 0
        lines_to_write = ""

        for tgt in xrange(n_pg):
            for src in xrange(len(self.conn_mat_orn_pg[tgt])):
                (src_gid, w) = self.conn_mat_orn_pg[tgt][src]

#        for src in xrange(n_orn):
#            for tgt in xrange(len(self.conn_mat_orn_pg[src])):
#                (tgt_gid, w) = self.conn_mat_orn_pg[src][tgt]
#                w = self.conn_mat_orn_pg[src][tgt][1]
#            for tgt in xrange(n_pg):
#                w = self.conn_mat_orn_pg[src, tgt]
                if (w > 0.0):
#                    src_gid = src + orn_offset
                    tgt_gid = tgt + pg_offset
                    lines_to_write += "%d\t%d\t%.6E\n" % (src_gid, tgt_gid, w)
                    conn_cnt += 1

        f_plain = open(self.params['conn_list_orn_pg'], 'w')
        first_line = "%d 3\n" % (conn_cnt)
        f_plain.write(first_line)
        f_plain.write(lines_to_write)
        f_plain.close()

    def connect_pg_mit_serial(self):
        """
        There are two types of PG cells that connect serially to MT cells:
        1) Those PG cells that receive direct excitation from ORNs project via serial dendro dendritic inhibitory (DDI) synapses to MT cells providing feedforward inhibition
        2) PG cells that receive both excitation from ORNs and excitation from MIT cell but inhibit via serial synapses other MT cells in the glomerulus. WTA mechanism

        This function also writes serial MIT ---< PG connections into file conn_list_mit_pg_serial.dat
        """
        n_mit_x = self.params['n_mit_x']
        n_mit_y = self.params['n_mit_y']
        n_pg_ff_inh = int(round(self.params['n_pg_x_serial'] * (1. - self.params['rel_ff_wta_inhibition']))) # these cells provide only feedforward inhibition
        n_pg_wta_inh = int(round(self.params['n_pg_x_serial'] * self.params['rel_ff_wta_inhibition']))# rel_ff_wta_inhibition * n_pg_x_serial get excitation from one mitral cell but inhibit all other mitral cells via serial connections
        n_pg_x = self.params['n_pg_x']
        n_pg_y = self.params['n_pg_y']
        n_pg_mit_serial = n_pg_ff_inh + (n_mit_x - 1) * n_pg_wta_inh # total number of incoming serial connections from PG cells to MT cells
        w_pg_mit_serial = self.params['w_pg_mit_serial_target'] / n_pg_mit_serial
        w_pg_mit_serial_sigma = w_pg_mit_serial * self.params['w_pg_mit_serial_sigma'] # convert percentage -> conductance
        w_mit_pg_serial = self.params['w_mit_pg_serial']
        w_mit_pg_serial_sigma = w_mit_pg_serial * self.params['w_mit_pg_serial_sigma'] # convert percentage -> conductance

        lines_to_write_pg_mit = ""
        lines_to_write_mit_pg = ""
        conn_cnt_pg_mit = 0
        conn_cnt_mit_pg = 0
        for glom in xrange(n_pg_y):
            src_offset = glom * n_pg_x + self.params['pg_offset']
            tgt_offset = glom * n_mit_x + self.params['mit_offset']
            for src in xrange(n_pg_ff_inh): # [0 ... n_pg_ff_inh]
                # 1) draw connections from PG cells that receive direct ORN input
                # for each of these src pg cells decide from which ORN group it receives excitation, just as in the function self.connect_orn_pg() above
                # divide the 'serial' population into n_mit_x groups
                group = src % n_mit_x
                # the group is the index of the target MT cell within this glomerulus
                tgt_gid = group + tgt_offset
                src_gid = src + src_offset
                w = rnd.normal(w_pg_mit_serial, w_pg_mit_serial_sigma)
                if (w > 0.0): # inhibitory weights are positive, inhibition is handled via the reversal potential of the synapse in NEURON
                    lines_to_write_pg_mit += "%d\t%d\t%.6E\n" % (src_gid, tgt_gid, w)
                    conn_cnt_pg_mit += 1
                
            # 2) PG cells that receive excitation by ORNs and MT cells inhibit all other MT cells via serial synapses
            for src in xrange(n_pg_ff_inh, self.params['n_pg_x_serial']): # [n_pg_ff_inh ... n_pg_x_serial]
                src_gid = src + src_offset
                tgt_mts = range(0, n_mit_x)
                tgt_mts.remove(src % n_mit_x)
                for tgt in tgt_mts:
                    tgt_gid = tgt + tgt_offset
                    w = rnd.normal(w_pg_mit_serial, w_pg_mit_serial_sigma)
                    if (w > 0.0): # even inhibitory weights are positive, inhibition is handled via the reversal potential of the synapse in NEURON
                        lines_to_write_pg_mit += "%d\t%d\t%.6E\n" % (src_gid, tgt_gid, w)
                        conn_cnt_pg_mit += 1
                # MIT ---< PG
                # the 'src' PG cell receives excitation from the one that it does not inhibit
                w = rnd.normal(w_mit_pg_serial, w_mit_pg_serial_sigma)
                if (w > 0.0): # even inhibitory weights are positive, inhibition is handled via the reversal potential of the synapse in NEURON
                    mit_gid = (src % n_mit_x) + tgt_offset
                    lines_to_write_mit_pg+= "%d\t%d\t%.6E\n" % (mit_gid, src_gid, w)
                    conn_cnt_mit_pg += 1


        # PG ---o MIT
        f_plain = open(self.params['conn_list_pg_mit_serial'], 'w')
        first_line = "%d 3\n" % (conn_cnt_pg_mit)
        f_plain.write(first_line)
        f_plain.write(lines_to_write_pg_mit)
        f_plain.close()

        # MIT ---< PG
        f_plain = open(self.params['conn_list_mit_pg_serial'], 'w')
        first_line = "%d 3\n" % (conn_cnt_mit_pg)
        f_plain.write(first_line)
        f_plain.write(lines_to_write_mit_pg)
        f_plain.close()


    def connect_pg_mit_reciprocal(self):
        """
        Within each glomerulus there is a population of PG cells having reciprocal dendro-dendritic connections with MT cells.
        This population is split up in two:
            1) one having only local reciprocal ddi connections with a single MT cell
            2) cells of the other have additionally ddi connections with all the other mitral cells in this glomerulus
            
        These PG cells also have serial connections to the MT cells to which they do not have DDI connections.
        """
        # src and target is an arbitrary definition here, since it's reciprocal :)
        n_mit_x = self.params['n_mit_x']
        n_mit_y = self.params['n_mit_y']
        n_pg_x = self.params['n_pg_x']
        n_pg_y = self.params['n_pg_y']
        n_pg_rec_local = int(round(self.params['n_pg_x_rec'] * (1. - self.params['rel_reciprocal_intraglom'])))
        n_pg_rec_total = self.params['n_pg_x_rec']
        n_pg_mit_rec = n_pg_rec_total + (n_mit_x - 1) * self.params['rel_reciprocal_intraglom'] * self.params['n_pg_x_rec']
        w_pg_mit_reciprocal = self.params['w_pg_mit_reciprocal_target'] / n_pg_mit_rec
        w_pg_mit_reciprocal_sigma = w_pg_mit_reciprocal * self.params['w_pg_mit_reciprocal_sigma']  # convert percentage -> conductance
        w_mit_pg_reciprocal = self.params['w_mit_pg_reciprocal'] # this is not scaled as only ~ 0.1 percent of the reciprocal PG population receive excitation also from other mitral cells
        w_mit_pg_reciprocal_sigma = w_mit_pg_reciprocal * self.params['w_mit_pg_reciprocal_sigma'] # convert percentage -> conductance
        lines_to_write_pg_mit = ""
        lines_to_write_mit_pg = ""
        conn_cnt_pg_mit = 0
        conn_cnt_mit_pg = 0

        for glom in xrange(n_pg_y):
            src_offset = glom * n_pg_x + self.params['pg_offset']# + round(self.params['alpha_pg_serial'] * self.params['n_pg_x']) # additional offset because of the population of 'serial' PG in that row of the PG array
            tgt_offset = glom * n_mit_x + self.params['mit_offset']
            for src in xrange(self.params['n_pg_x_serial'], self.params['n_pg_x_serial'] + n_pg_rec_local): # these cells have local reciprocal ddi connections
                # source PG cells are assigned to a distinct MT cell 'round-robin'
                tgt_mt = src % n_mit_x
                tgt_gid = tgt_mt + tgt_offset
                src_gid = src + src_offset
                # PG ---o MT
                w_pg_mit = rnd.normal(w_pg_mit_reciprocal, w_pg_mit_reciprocal_sigma)
                if (w_pg_mit > 0.0): # even inhibitory weights are positive, inhibition is handled via the reversal potential of the synapse in NEURON
                    lines_to_write_pg_mit += "%d\t%d\t%.6E\n" % (src_gid, tgt_gid, w_pg_mit)
                    conn_cnt_pg_mit += 1

                # MT ---< PG
                w_mit_pg = rnd.normal(w_mit_pg_reciprocal, w_mit_pg_reciprocal_sigma)
                if (w_mit_pg > 0.0):
                    lines_to_write_mit_pg += "%d\t%d\t%.6E\n" % (tgt_gid, src_gid, w_mit_pg) #tgt and source swapped
                    conn_cnt_mit_pg += 1

            # these cells have BOTH local reciprocal ddi connections and intraglomerular, i.e. with other mitral cells, too
            for src in xrange(self.params['n_pg_x_serial'] + n_pg_rec_local, self.params['n_pg_x']): 
                src_gid = src + src_offset
                tgt_mts = range(0, n_mit_x)
                tgt_mts.remove(src % n_mit_x) # remove the one mt cell with which this cell already has a reciprocal connection to
                for tgt in tgt_mts:
                    tgt_gid = tgt + tgt_offset
                    # PG ---o MT
                    w_pg_mit = rnd.normal(w_pg_mit_reciprocal, w_pg_mit_reciprocal_sigma)
                    if (w_pg_mit > 0.0): # even inhibitory weights are positive, inhibition is handled via the reversal potential of the synapse in NEURON
                        lines_to_write_pg_mit += "%d\t%d\t%.6E\n" % (src_gid, tgt_gid, w_pg_mit)
                        conn_cnt_pg_mit += 1
                    # MT ---< PG
                    w_mit_pg = rnd.normal(w_mit_pg_reciprocal, w_mit_pg_reciprocal_sigma)
                    if (w_mit_pg > 0.0):
                        lines_to_write_mit_pg += "%d\t%d\t%.6E\n" % (tgt_gid, src_gid, w_mit_pg) #tgt and source swapped
                        conn_cnt_mit_pg += 1


        print "conn_cnt_pg_mit_reciprocal = %d\t Number of reciprocal pg-mit connections" % (conn_cnt_pg_mit)
        print "conn_cnt_mit_pg_reciprocal = %d\t Number of reciprocal mit-pg connections" % (conn_cnt_mit_pg)
        f_plain = open(self.params['conn_list_pg_mit_reciprocal'], 'w')
        first_line = "%d 3\n" % (conn_cnt_pg_mit)
        f_plain.write(first_line)
        f_plain.write(lines_to_write_pg_mit)
        f_plain.close()

        f_plain = open(self.params['conn_list_mit_pg_reciprocal'], 'w')
        first_line = "%d 3\n" % (conn_cnt_mit_pg)
        f_plain.write(first_line)
        f_plain.write(lines_to_write_mit_pg)
        f_plain.close()


    def connect_mt_gran_local(self):
        """
        This function draws excitatory connections from MT cells to Gran cells residing in the same glomerular column
        """
        n_mit_x = self.params['n_mit_x']
        n_mit_y = self.params['n_mit_y']
        n_gran_x = self.params['n_gran_x']
        n_gran_y = self.params['n_gran_y']
        n_local_syn = self.params['n_mit_gran_syn_local']
        # scale the target weights to cell-to-cell connections
        w_gran_mit = self.params['w_gran_mit_local_target'] / n_local_syn
        w_gran_mit_sigma = w_gran_mit * self.params['w_gran_mit_local_sigma'] # convert percentage -> conductance
        w_mit_gran = (self.params['w_mit_gran_local_target'] / n_local_syn) * (n_gran_x / float(n_mit_x))
        w_mit_gran_sigma = w_mit_gran * self.params['w_mit_gran_local_sigma'] # convert percentage -> conductance
#        mit_gran_output = []
        lines_to_write_mit_gran = ""
        lines_to_write_gran_mit = ""
        conn_cnt_mit_gran = 0
        conn_cnt_gran_mit = 0
        for glom in xrange(n_mit_y):
            for mit in xrange(n_mit_x): # 'src'
                granule_cells = numpy.unique(rnd.randint(0, n_gran_x, n_local_syn))
                for gran in granule_cells:
                    mit_gid = glom * n_mit_x + self.params['mit_offset'] + mit
                    # before n_gran_y
                    gran_gid = glom * n_gran_x + self.params['gran_offset'] + gran
                    w_exc = rnd.normal(w_mit_gran, w_mit_gran_sigma)
                    if (w_exc > 0.0):
                        lines_to_write_mit_gran += "%d\t%d\t%.6E\n" % (mit_gid, gran_gid, w_exc)
                        conn_cnt_mit_gran += 1
                    w_inh = rnd.normal(w_gran_mit, w_gran_mit_sigma)
                    if (w_inh > 0.0):
                        lines_to_write_gran_mit += "%d\t%d\t%.6E\n" % (gran_gid, mit_gid, w_inh)
                        conn_cnt_gran_mit += 1


        f_mit_gran_plain = open(self.params['conn_list_mit_gran_local'], 'w')
        first_line = "%d 3\n" % (conn_cnt_mit_gran)
        f_mit_gran_plain.write(first_line)
        f_mit_gran_plain.write(lines_to_write_mit_gran)
        f_mit_gran_plain.close()
        f_gran_mit_plain = open(self.params['conn_list_gran_mit_local'], 'w')
        first_line = "%d 3\n" % (conn_cnt_gran_mit)
        f_gran_mit_plain.write(first_line)
        f_gran_mit_plain.write(lines_to_write_gran_mit)
        f_gran_mit_plain.close()


    def connect_mt_gran_global(self):
        """
        This function draws excitatory connections from MT cells to Gran cells residing in a different glomerular column
        """
        n_mit_x = self.params['n_mit_x']
        n_mit_y = self.params['n_mit_y']
        n_gran_x = self.params['n_gran_x']
        n_gran_y = self.params['n_gran_y']
        n_global_syn = self.params['n_mit_gran_syn_global']
        # scale the target weights to cell-to-cell connections
        w_gran_mit = self.params['w_gran_mit_global_target'] / n_global_syn
        w_gran_mit_sigma = w_gran_mit * self.params['w_gran_mit_global_sigma'] # convert percentage -> conductance
        w_mit_gran = (self.params['w_mit_gran_global_target'] / n_global_syn) * (n_gran_x / float(n_mit_x))
        w_mit_gran_sigma = w_mit_gran * self.params['w_mit_gran_global_sigma'] # convert percentage -> conductance
        mit_offset = self.params['mit_offset']
        gran_offset = self.params['gran_offset']
#        mit_gran_output = []
        lines_to_write_mit_gran = ""
        lines_to_write_gran_mit = ""
        conn_cnt_mit_gran = 0
        conn_cnt_gran_mit = 0
        for src_glom in xrange(n_mit_y):
            tgt_gloms = range(n_mit_y)
            tgt_gloms.remove(src_glom)
            for mit in xrange(n_mit_x): # 'src'
                for syn in xrange(n_global_syn):
                    tgt_glom = random.sample(tgt_gloms, 1)[0] # draw a random target glomerulus other than the source glom
                    tgt_gran = rnd.randint(0, n_gran_x)
                    mit_gid = src_glom * n_mit_x + mit_offset + mit
                    gran_gid = tgt_glom * n_gran_x + gran_offset + tgt_gran
                    w_exc = rnd.normal(w_mit_gran, w_mit_gran_sigma)
                    if (w_exc > 0.0):
                        lines_to_write_mit_gran += "%d\t%d\t%.6E\n" % (mit_gid, gran_gid, w_exc)
                        conn_cnt_mit_gran += 1
                    w_inh = rnd.normal(w_gran_mit, w_gran_mit_sigma)
                    if (w_inh > 0.0):
                        lines_to_write_gran_mit += "%d\t%d\t%.6E\n" % (gran_gid, mit_gid, w_inh)
                        conn_cnt_gran_mit += 1

        f_mit_gran_plain = open(self.params['conn_list_mit_gran_global'], 'w')
        first_line = "%d 3\n" % (conn_cnt_mit_gran)
        f_mit_gran_plain.write(first_line)
        f_mit_gran_plain.write(lines_to_write_mit_gran)
        f_mit_gran_plain.close()
        f_gran_mit_plain = open(self.params['conn_list_gran_mit_global'], 'w')
        first_line = "%d 3\n" % (conn_cnt_gran_mit)
        f_gran_mit_plain.write(first_line)
        f_gran_mit_plain.write(lines_to_write_gran_mit)
        f_gran_mit_plain.close()

    def change_orn_mit_weights(self, mit_gid, factor):
        """
        Multiply all weights of ORN - MT connections targeting cell with mit_gid by factor
        """
        tgt_index = mit_gid - self.params['mit_offset']
#        self.conn_mat_orn_mit[:, tgt_index] *= factor
        for src in xrange(len(self.conn_mat_orn_mit[tgt_index])):
            src_gid = self.conn_mat_orn_mit[tgt_index][src][0] 
            w = self.conn_mat_orn_mit[tgt_index][src][1] 
            w_new = w * factor
            self.conn_mat_orn_mit[tgt_index][src] = (src_gid, w_new)

    def change_orn_pg_weights(self, target_group_index, factor):
        """
        Change all ORN-PG weights originating from the group of ORNs that excite the MT cell(s) with the target_group_index
        withing the row.
        """
        n_tgt_pg = int(round(self.params['n_pg_x_serial'] * (1. - self.params['rel_ff_wta_inhibition']))) # these cells provide only feedforward inhibition
        n_pg_x = self.params['n_pg_x']
        n_pg_y = self.params['n_pg_y']
        n_orn_x = self.params['n_orn_x']
        pg_offset = self.params['pg_offset']

        for glom in xrange(n_pg_y): # for each row
            src_offset = glom * n_orn_x + self.params['orn_offset']
            tgt_offset = glom * n_pg_x + self.params['pg_offset']
            for tgt in xrange(n_tgt_pg):
                # for each of these tgt pg cells decide from which ORN it receives excitation
                # divide the 'serial' population into n_mit_x groups, each group receives excitation from a different set of ORNs (still all ORNs express the same OR)
                group = tgt % self.params['n_mit_x']
                if (group == target_group_index):
                    tgt_gid = tgt + tgt_offset
                    conn_list = self.conn_mat_orn_pg[tgt_gid - pg_offset]
                    for i in xrange(len(conn_list)):
                        (src_gid, w) = self.conn_mat_orn_pg[tgt_gid - pg_offset][i]
                        w_new = w * factor
                        self.conn_mat_orn_pg[tgt_gid - pg_offset][i] = (src_gid, w_new)

