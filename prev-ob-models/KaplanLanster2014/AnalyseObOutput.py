import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import MergeSpikefiles
#import functions

class AnalyseObOutput(object):
    """
    This classe calculates the normalized MIT cell responses after stimulating OB with EPTH input.
    Usage: pass the parameter dictionary to the class constructor
    The spiking activity of MIT cells within one hypercolumn (HC) (i.e. one glomerulus getting input from 
    a population of ORNs expressing one OR) in response to one pattern presentation is mapped into (0,1) and normalized.
    Normalized means that the sum of the activity for one cell is within 0 and 1 AND that summed activity for all cells within one glomerulus is within 0 and 1.
    E.g. 
    MIT spike output
          \  pattern 
    cell   \    0       1    
    -----------------------
        0  |    0       10 
        1  |    5       10
        2  |    8       2

    Normalized cell activity:
          \  pattern
    cell   \    0   1    
    -----------------------
        0  |    0       1.0
        1  |    0.33    0.67
        2  |    0.8     0.2

    After that a second normalization step is necessary, which guarantees that the sum of activities of cells belonging to one glomerulus is < 1

          \  pattern
    cell   \    0   1    
    -----------------------
        0  |    0       0.536
        1  |    0.294   0.357
        2  |    0.706   0.107

    """

    def __init__(self, param_dict):
#        self.sim_param_fn = sim_param_fn
        self.params = param_dict
        self.mit_spike_output = None # [pattern, cellid] gives number of spikes
        self.Merger = MergeSpikefiles.MergeSpikefiles(self.params)

    def get_output_activity(self):
        """
        This function simply loads the number of spikes fired by each mitral cell and return an array, e.g.
            mit_spikes[pattern, cell_index]
        """
        n_patterns = self.params['n_patterns']
        n_mit = self.params['n_mit']
        mit_offset = self.params['mit_offset']
        self.mit_spikes = np.zeros((n_patterns, n_mit)) # stores the number of spikes

        missing_patterns = []
        # copy the number of spikes for each pattern into an array
        for pattern in xrange(self.params['n_patterns']):
            fn = self.params['mit_spikes_merged_fn_base'] + str(pattern) + '.dat'
            if (os.path.exists(fn) == False):
                self.Merger.merge_nspike_files(self.params['mit_spike_fn_base'], self.params['mit_spikes_merged_fn_base'], pattern)
                print "Merging nspike files for pattern", pattern

            input_fn = self.params['mit_spikes_merged_fn_base'] + str(pattern) + '.dat'
#            print "debug fn", input_fn
            try:
                single_pattern_response = np.loadtxt(input_fn)
                for row in xrange(single_pattern_response[:, 0].size):
                    gid = single_pattern_response[row, 0]
                    cell_index = gid - mit_offset
                    self.mit_spikes[pattern, cell_index] = single_pattern_response[row, 1]
            except:
                missing_patterns.append(pattern)
#            single_pattern_response = functions.fill_up_nspike_data(input_fn, n_mit, mit_offset)
        if (len(missing_patterns) > 0):
            print "Loading nspike files from the following patterns failed!\nEither delete nspikes_*_merged_ files and retry or resimulate missing patterns!", missing_patterns
        
        output_fn = self.params['mit_response']
        print "Saving mit nspike matrix to:", output_fn
        np.savetxt(output_fn, self.mit_spikes)

    def rescale_activity(self):
        """
        Rescale mit_spikes so that the global maximum = 1
        assume that get_output_activity was called before!
        """
        n_max = 0
        for pn in xrange(self.params['n_patterns']):
            n_max = max(n_max, self.mit_spikes[pn, :].max())
        mit_spikes_new = self.mit_spikes.copy() / n_max
        output_fn = self.params['mit_nspikes_rescaled']
        print "AnalyseObOutput output file:", output_fn
        np.savetxt(output_fn, mit_spikes_new)

    def rescale_activity_cellwise(self):
        """
        Rescale mit_spikes so that for all mitral cells holds: summed activity over all patterns = 1
        assume that get_output_activity was called before!
        """

        mit_spikes_new = np.zeros((self.params['n_patterns'], self.params['n_mit']))
        for mit in xrange(self.params['n_mit']):
            spike_sum = self.mit_spikes[:, mit].sum()
            if spike_sum > 0:
                mit_spikes_new[:, mit] = self.mit_spikes[:, mit] / spike_sum
        output_fn = self.params['mit_nspikes_normed_cells']
        print "AnalyseObOutput output file:", output_fn
        np.savetxt(output_fn, mit_spikes_new)

    def rescale_activity_patternwise(self):
        """
        Rescale mit_spikes so that for each pattern holds: summed activity over all mitral cells = 1
        assume that get_output_activity was called before!
        """

        mit_spikes_new = np.zeros((self.params['n_patterns'], self.params['n_mit']))
        for pn in xrange(self.params['n_patterns']):
            mit_spikes_new[pn, :] = self.mit_spikes[pn, :] / self.mit_spikes[pn, :].sum()
        output_fn = self.params['mit_nspikes_normed_patterns']
        print "AnalyseObOutput output file:", output_fn
        np.savetxt(output_fn, mit_spikes_new)


    def rescale_activity_glom_patterns(self):
        """
        Rescale mit_spikes with two normalization steps:
            1) for all patterns: summed activity within a glomerulus = 1 (or 0)
            2) for all patterns: summed activity for all cells = 1
        assume that get_output_activity was called before!
        """

        mit_spikes_new = np.zeros((self.params['n_patterns'], self.params['n_mit']))
        # glomerular normalization first
        for pn in xrange(self.params['n_patterns']):
            for glom in xrange(self.params['n_mit_y']):
                mit_spikes_new[pn, :] = self.mit_spikes[pn, :] / self.mit_spikes[pn, :].sum()
                    # calculate which cells belong to one glomerulus
                id1 = glom * self.params['n_mit_x']
                id2 = (glom + 1) * self.params['n_mit_x']  
                # normalization within one glomerulus
                activity_sum = mit_spikes_new[pn, id1:id2].sum()
                if (activity_sum > 1):
                    mit_spikes_new[pn, id1:id2] /= activity_sum

        # Normalization: the sum of spikes fired by each cell during all patterns is normalized to 1 -> each mitral cell has a normalized activty
        for cell in xrange(self.params['n_mit']):
            spike_sum = mit_spikes_new[:, cell].sum() # sum of spikes for one cell over all patterns
            if (spike_sum != 0):
                mit_spikes_new[:, cell] = mit_spikes_new[:, cell] / spike_sum

        output_fn = self.params['mit_nspikes_normed_glom_cells']
        print "AnalyseObOutput output file:", output_fn
        np.savetxt(output_fn, mit_spikes_new)

    def rescale_activity_cellwise_then_patternwise(self):
        """
        0) Assume that get_output_activity was called before.
        1) Rescale mit_spikes so that for all mitral cells holds: summed activity over all patterns = 1
        2) Rescale mit_spikes so that for each pattern holds: summed activity over all mitral cells = 1
        """

        mit_spikes_new = np.zeros((self.params['n_patterns'], self.params['n_mit']))
        for mit in xrange(self.params['n_mit']):
            spike_sum = self.mit_spikes[:, mit].sum()
            if spike_sum > 0:
                mit_spikes_new[:, mit] = self.mit_spikes[:, mit] / spike_sum
        for pn in xrange(self.params['n_patterns']):
            mit_spikes_new[pn, :] = self.mit_spikes[pn, :] / self.mit_spikes[pn, :].sum()
        output_fn = self.params['mit_nspikes_normed_cells_then_patterns']
        print "AnalyseObOutput output file:", output_fn
        np.savetxt(output_fn, mit_spikes_new)


    def rescale_activity_patternwise_then_cellwise(self):
        """
        0) Assume that get_output_activity was called before.
        1) Rescale mit_spikes so that for each pattern holds: summed activity over all mitral cells = 1
        2) Rescale mit_spikes so that for all mitral cells holds: summed activity over all patterns = 1
        """

        mit_spikes_new = np.zeros((self.params['n_patterns'], self.params['n_mit']))
        for pn in xrange(self.params['n_patterns']):
            mit_spikes_new[pn, :] = self.mit_spikes[pn, :] / self.mit_spikes[pn, :].sum()
        for mit in xrange(self.params['n_mit']):
            spike_sum = self.mit_spikes[:, mit].sum()
            if spike_sum > 0:
                mit_spikes_new[:, mit] = self.mit_spikes[:, mit] / spike_sum
        output_fn = self.params['mit_nspikes_normed_patterns_then_cells']
        print "AnalyseObOutput output file:", output_fn
        np.savetxt(output_fn, mit_spikes_new)




#    def get_output_file(self, normalize=False, remove_silent_mit=False):
    def get_mit_response_normalized(self):
        """
        This function will return the file with the normalized activity in the olfactory bulb
        in the format: pattern (row) vs. mitral cell (column)
            mit1    mit2
        p1
        p2

        For this purpose, 
            1) merge the mitral cell output files: 
                first merge files from different processors
                sort according to cell id
                then merge files for different patterns
            2) Normalization:
                a) the sum of spikes fired by each cell during all patterns is normalized to 1 
                -> each mitral cell has a normalized activty
                b) if the sum of normalized activity of mitral cells within one glomerular unit > 1 -> set it to one

            glomerulus (hypercolumn) equals one
        """
        print "Normalizing OB output patterns..."

        n_patterns = self.params['n_patterns']
        n_mit = self.params['n_mit']
        mit_offset = self.params['mit_offset']
        mit_spikes = np.zeros((n_patterns, n_mit)) # stores the number of spikes

        # copy the number of spikes for each pattern into an array
        for pattern in xrange(self.params['n_patterns']):
            fn = self.params['mit_spikes_merged_fn_base'] + str(pattern) + '.dat'
            if (os.path.exists(fn) == False):
                self.Merger.merge_nspike_files(self.params['mit_spike_fn_base'], self.params['mit_spikes_merged_fn_base'], pattern)
                print "Merging nspike files for pattern", pattern

            assert(os.path.exists(self.params['mit_spikes_merged_fn_base'] + str(pattern) + '.dat')), 'ERROR:\n\tFile does not exist %s' % (self.params['mit_spikes_merged_fn_base'] + str(pattern) + '.dat')
            if self.params['print_debug']:
                print 'AnalyseObOutput: loading:', self.params['mit_spikes_merged_fn_base'] + str(pattern) + '.dat'
            single_pattern_response = np.loadtxt(self.params['mit_spikes_merged_fn_base'] + str(pattern) + '.dat')

            if single_pattern_response.size > 0:
                for row in xrange(single_pattern_response[:, 0].size):
                    gid = single_pattern_response[row, 0]
                    cell_index = gid - mit_offset
                    mit_spikes[pattern, cell_index] = single_pattern_response[row, 1]

        # 2a) Normalization: the sum of spikes fired by each cell during all patterns is normalized to 1 -> each mitral cell has a normalized activty
        normalized_activity = np.zeros((n_patterns, n_mit)) # normalized_activity[:, cell_id].sum() = 1 for all cell_id (or = 0 if no spikes fired)
        for cell in xrange(n_mit):
            spike_sum = mit_spikes[:, cell].sum() # sum of spikes for one cell over all patterns
            if (spike_sum == 0):
                normalized_activity[:, cell] = np.zeros(n_patterns)
            else:
                normalized_activity[:, cell] = mit_spikes[:, cell] / spike_sum

        # 2b) the sum of normalized_activity within one glomerulus must not be > 1

        for pattern in xrange(n_patterns):
            for glom in xrange(self.params['n_mit_y']):
                # calculate which cells belong to one glomerulus
                id1 = glom * self.params['n_mit_x']
                id2 = (glom + 1) * self.params['n_mit_x']  
                # normalization within one glomerulus
                activity_sum = normalized_activity[pattern, id1:id2].sum()
                if (activity_sum > 1):
                    normalized_activity[pattern, id1:id2] /= activity_sum
                # else: do nothing

        wta_activity = np.zeros((n_patterns, n_mit))
        for pattern in xrange(n_patterns):
            for glom in xrange(self.params['n_mit_y']):
                id1 = glom * self.params['n_mit_x']
                id2 = (glom + 1) * self.params['n_mit_x']  
                if mit_spikes[pattern, id1:id2].sum() > 0:
                    wta_activity[pattern, id1 + np.argmax(mit_spikes[pattern, id1:id2])] = 1

        wta_output_fn = self.params['other_folder'] + '/mit_activity_wta.dat'
        print 'Saving mit wta activity to:', wta_output_fn
        np.savetxt(wta_output_fn, wta_activity)

                     
        silent_mit = []
        for mit in xrange(n_mit):
            if (normalized_activity[:, mit].sum() == 0):
                silent_mit.append(mit)

        if (len(silent_mit) > 0):
            silent_cells  = ['%d ' % mit for mit in silent_mit]
            silent_file = open(self.params['silent_mit_fn'], 'w')
            for s in silent_cells:
                silent_file.write(s)
            silent_file.close()

        print 'Silent mitral cells:', silent_mit
        print "AnalyseObOutput output file:", self.params["mit_response_normalized"]
        np.savetxt(self.params["mit_response_normalized"], normalized_activity)

#        spike_column = 1 # column which contains number of spikes
#                spike_sum = single_pattern_response[id1:id2, spike_column].sum()
#                if ((spike_sum <= 1)):
#                    mit_response_data[pattern, id1:id2] = np.zeros(id2-id1)
#                elif (normalize):
#                    mit_response_data[pattern, id1:id2] = single_pattern_response[id1:id2, spike_column] / spike_sum
#                else:
#                    mit_response_data[pattern, id1:id2] = single_pattern_response[id1:id2, spike_column]
#            self.mit_spike_output[pattern] = mit_response_data[pattern]
#        np.savetxt(self.params["mit_response_normalized"], mit_response_data, delimiter=",")
#        if normalize:
#            print "Normalize", self.params["mit_response_normalized"]
#            np.savetxt(self.params["mit_response_normalized"], mit_response_data, delimiter=",")
#        else:
#            print "Not normalized ", self.params["mit_response"]
#            np.savetxt(self.params["mit_response"], mit_response_data, delimiter=",")

    def write_active_cells_to_log_file(self):
        """
        Return the ids of MT cells which have a normalized output activity higher than average + 1 * sigma
        """
        output_fn = self.params["analysis_output_fn"] # = "%s/analysis_output.log" % (self.params["folder_name"])
        output_f = open(output_fn, 'a')  # open file in append mode
        n_patterns = self.params['n_patterns']
        active_mit_cells = [[] for i in xrange(n_patterns)]   # contains all MT cells firing at least one spike
        active_mit_cells_thresholded = [[] for i in xrange(n_patterns)] # contains all MT cells firing more than average n_spikes
        lines = ""

        for p in xrange(n_patterns):
            mean_nspikes = self.mit_spike_output[p].mean()
            std_nspikes = self.mit_spike_output[p].std()
            for mit in xrange(self.params['n_mit']):
                if (self.mit_spike_output[p, mit] > 0):
                    active_mit_cells.append(mit) # + self.params['n_orn'] + self.params['global_offset']
                if (self.mit_spike_output[p, mit] > mean_nspikes):
                    active_mit_cells_thresholded.append(mit) # + self.params['n_orn'] + self.params['global_offset']

         
            lines += "Pattern: " + str(p)
            lines += "\nactive_mit_cells:" 
            lines += str(active_mit_cells)
            lines += "\nn_active_mit_cells:" + str(len(active_mit_cells)) + "\t" + str(len(active_mit_cells) / float(self.params['n_mit']))
            lines += "\nactive_mit_cells_thresholded:"
            lines += str(active_mit_cells_thresholded)
            lines += "\nn_active_mit_cells_thresholded:" + str(len(active_mit_cells_thresholded)) + "\t" + str(len(active_mit_cells_thresholded) / float(self.params['n_mit']))

        print lines
        output_f.write(lines)
        output_f.close()
            # build the correlation matrix, pij
#            pij = np.zeros((n_mit, n_mit))
#            for i in xrange(n_mit):
#                for j in xrange(n_mit):
#                    for pattern in xrange(n_patterns):
#                        pij[i,j] += mit_response_data[i, pattern] * mit_response_data[j,pattern]
#                    pij[i,j] /= n_patterns
