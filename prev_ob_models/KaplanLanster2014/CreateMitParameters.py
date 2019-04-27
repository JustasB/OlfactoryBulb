import numpy
import sys
import re

class CreateMitParameters(object):
    def __init__(self, param_dict):
        self.params = param_dict
        self.conn_list_orn_mit = self.params['conn_list_orn_mit']
        self.global_offset = self.params['global_offset']

    def create_parameters(self, test=0):
        """
        if (test == 1):
            it is assumed, that only one pattern is presented
        output:
            gid     oor     c/Kd    glom

            oor and c/Kd are mean values of the oor and c/Kd values from neurons projecting to the MIT
        """
        print "Creating mit parameters..."
        # orn_params contains in first line n_rows n_cols
        for pn in xrange(self.params['n_patterns']):
            if (test==1):
                orn_pf = file(self.params["orn_params_test"], 'r')
            else: 
                orn_pf = file(self.params['orn_params_fn_base'] + "%d.dat" % (pn), 'r')
            firstline = orn_pf.readline()
            values = re.split("\s", firstline)
            n_rows = int(values[0])
            n_cols = int(values[1])
            #print "n_rows, n_cols", n_rows, n_cols
            orn_params = numpy.zeros((n_rows, n_cols)) # (n_orn_x * n_orn_y, n_params)

            # read in the orn_params from file
            row = 0
            for line in orn_pf:
                values = re.split("\s", line)
                for col in xrange(n_cols):
                    orn_params[row][col] = values[col]
                row += 1

            # the conn_list_ file stores connections from orn -> mit cells in the form:
            # src tgt weight
            # this list is sorted according to the tgt cell index
            # for each tgt cell the average parameters for the cells projecting excitatory to the tgt cell
            # are now being calculated
            conn_file = open(self.conn_list_orn_mit, 'r')
            print "Reading data from :\n\t", self.conn_list_orn_mit

            firstline = conn_file.readline()
            #values = re.split("\s", firstline)
            #n_rows_cf = int(values[0])
            #n_cols_cf = int(values[1])
            self.n_mit_params = n_cols + 1 # one more for glomerulus id
            params = numpy.zeros(n_cols-1) # buffer for parameters for one target cell
            mit_params = numpy.zeros((self.params['n_mit'], self.n_mit_params))

            # process the first data line
            first_dataline = conn_file.readline()
            values = re.split("\s", first_dataline)
            src, tgt, weight = int(values[0]), int(values[1]), float(values[2])
            old_tgt = tgt
            tgt_cnt = 0
            src_cnt = 0
            if (weight > 0):
                params += orn_params[0,1:]
                src_cnt = 1 # count the cells projecting excitatory to one tgt 

            for line in conn_file:
                values = re.split("\s", line)
                src, tgt, weight = int(values[0]), int(values[1]), float(values[2])
            #    if ((tgt != old_tgt) and (weight < 0)): # there is another tgt cell or the same tgt cell but with inhibitory connections
                if ((tgt != old_tgt)): # there is another tgt cell
                    # calculate the mean parameters for all source cells

                    mit_params[tgt_cnt,2:] = params / src_cnt
#                    glom_id = (tgt - self.n_orn - self.global_offset) / self.n_mit_x
                    glom_id = (tgt_cnt) / self.params['n_mit_x']
                    mit_params[tgt_cnt,1] = glom_id
                    mit_params[tgt_cnt,0] = old_tgt
            #        if (weight > 0):
                    tgt_cnt += 1 # new mitral cell
                    src_cnt = 0 # reset src_cnt
                    params *= 0 # reset params array
                if (weight > 0.0):
            #        print "src %d \tsrc-glof %d\torn_params.shape %d %d" % (src, src-global_offset, orn_params.shape[0], orn_params.shape[1])
            #        print "orn: %d " % (src-global_offset), orn_params[src-global_offset,:]
                    params += orn_params[src-self.global_offset,1:]
                    src_cnt += 1
                old_tgt = tgt   # update
                
            # process the last parameter package
            mit_params[tgt_cnt,2:] = params / src_cnt
#            glom_id = (tgt - self.n_orn - self.global_offset) / self.n_mit_x
            glom_id = (tgt_cnt) / self.params['n_mit_x']
            mit_params[tgt_cnt,1] = glom_id
            mit_params[tgt_cnt,0] = old_tgt

            # write output to mit parameter file
            # firstline has to be NEURON conform
            output_fn = self.params['mit_params_fn_base'] + "%d.dat" % pn
            print "\tWriting to ", output_fn
            firstline = "%d %d\n" % (self.params['n_mit'], n_cols)
            mit_pf = file(output_fn, 'w')
            mit_pf.write(firstline)
            for mit in xrange(self.params['n_mit']):
                line = "%d\t" % (mit_params[mit,0])
                for param in xrange(1, self.n_mit_params):
                    line += "%.4e\t" % mit_params[mit, param]
#                    mit_pf.write("%.4e\t" % mit_params[mit, param])
                mit_pf.write(line)
                mit_pf.write("\n")
            mit_pf.close()
        return 1
        if (test == 1):
            self.n_pattern = 1
