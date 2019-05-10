import os
import numpy
import sys
import simulation_parameters

class MergeSpikefiles(object):

    def __init__(self, params):
        self.params = params


    def merge_nspike_files(self, merge_pattern, sorted_pattern_output, pattern='', sort_idx=1):
        if (pattern == ''): # merge for all patterns
            for pattern in xrange(self.params['n_patterns']):
                rnd_nr1 = numpy.random.randint(0,10**8)
                rnd_nr2 = rnd_nr1 + 1
                fn_out = sorted_pattern_output + str(pattern) + ".dat"
                print 'output_file:', fn_out
                # merge files from different processors
                tmp_file = "tmp_%d" % (rnd_nr2)
                os.system("cat %s%d_* > %s" % (merge_pattern, pattern, tmp_file))
                # sort according to cell id
                os.system("sort -gk %d %s > %s" % (sort_idx, tmp_file, fn_out))
                os.system("rm %s" % (tmp_file))
        else:
            rnd_nr1 = numpy.random.randint(0,10**8)
            rnd_nr2 = rnd_nr1 + 1
            fn_out = sorted_pattern_output + str(pattern) + ".dat"
            print 'output_file:', fn_out
            # merge files from different processors
            tmp_file = "tmp_%d" % (rnd_nr2)
            os.system("cat %s%d_* > %s" % (merge_pattern, pattern, tmp_file))
            # sort according to cell id
            os.system("sort -gk %d %s > %s" % (sort_idx, tmp_file, fn_out))
            os.system("rm %s" % (tmp_file))


    def merge_spiketimes_files(self, merge_pattern, sorted_pattern_output, pattern='', sort_idx=1):
        if (pattern == ''): # merge for all patterns
            for pattern in xrange(self.params['n_patterns']):
                rnd_nr1 = numpy.random.randint(0,10**8)
                rnd_nr2 = numpy.random.randint(0,10**8) + 1
                fn_out = sorted_pattern_output + str(pattern) + ".dat"
                print 'output_file:', fn_out
                # merge files from different processors
                tmp_file = "tmp_%d" % (rnd_nr2)
                os.system("cat %s%d_* > %s" % (merge_pattern, pattern, tmp_file))
                # sort according to cell id
                os.system("sort -gk %d %s > %s" % (sort_idx, tmp_file, fn_out))
                os.system("rm %s" % (tmp_file))
        else:
            rnd_nr1 = numpy.random.randint(0,10**8)
            rnd_nr2 = numpy.random.randint(0,10**8) + 1
            fn_out = sorted_pattern_output + str(pattern) + ".dat"
            print 'output_file:', fn_out
            # merge files from different processors
            tmp_file = "tmp_%d" % (rnd_nr2)
            os.system("cat %s%d_* > %s" % (merge_pattern, pattern, tmp_file))
            # sort according to cell id
            os.system("sort -gk %d %s > %s" % (sort_idx, tmp_file, fn_out))
            os.system("rm %s" % (tmp_file))


    def merge_epth_spiketimes_file(self, pattern=''):
        merge_pattern = self.params["orn_spiketimes_fn_base"]
        sorted_pattern_output = self.params["orn_spiketimes_merged_fn_base"]
        self.merge_spiketimes_files(merge_pattern, sorted_pattern_output, pattern)


    def merge_epth_nspike_files(self, pattern=''):
        merge_pattern = self.params["orn_spike_fn_base"]
        sorted_pattern_output = self.params["orn_spikes_merged_fn_base"]
        self.merge_nspike_files(merge_pattern, sorted_pattern_output, pattern)


    def merge_ob_spiketimes_file(self, pattern='', sort_idx=1):
        merge_pattern = self.params["mit_spiketimes_fn_base"]
        sorted_pattern_output = self.params["mit_spiketimes_merged_fn_base"]
        self.merge_spiketimes_files(merge_pattern, sorted_pattern_output, pattern)

        merge_pattern = self.params["gran_spiketimes_fn_base"]
        sorted_pattern_output = self.params["gran_spiketimes_merged_fn_base"]
        self.merge_spiketimes_files(merge_pattern, sorted_pattern_output, pattern)

        merge_pattern = self.params["pg_spiketimes_fn_base"]
        sorted_pattern_output = self.params["pg_spiketimes_merged_fn_base"]
        self.merge_spiketimes_files(merge_pattern, sorted_pattern_output, pattern)

        if (pattern==''):
            # create spiketimes file from the other three
            for pn in xrange(self.params['n_patterns']):
                mit_fn = self.params['mit_spiketimes_merged_fn_base'] + str(pn) + '.dat'
                gran_fn = self.params['gran_spiketimes_merged_fn_base'] + str(pn) + '.dat'
                pg_fn = self.params['pg_spiketimes_merged_fn_base'] + str(pn) + '.dat'
                ob_fn = self.params['ob_spikes_merged_fn_base'] + str(pn) + '.dat'

                rnd_nr = numpy.random.randint(0,10**8)
                tmp_fn = 'tmp_%d' % (rnd_nr)
                os.system('cat %s %s %s > %s' % (mit_fn, pg_fn, gran_fn, tmp_fn))
                os.system('sort -gk %d %s > %s' % (sort_idx, tmp_fn, ob_fn))
                os.system('rm %s' % tmp_fn)


    def merge_ob_nspike_files(self, pattern='', sort_idx=1):
        merge_pattern = self.params["mit_spike_fn_base"]
        sorted_pattern_output = self.params["mit_spikes_merged_fn_base"]
        self.merge_nspike_files(merge_pattern, sorted_pattern_output, pattern)

        merge_pattern = self.params["gran_spike_fn_base"]
        sorted_pattern_output = self.params["gran_spikes_merged_fn_base"]
        self.merge_nspike_files(merge_pattern, sorted_pattern_output, pattern)

        merge_pattern = self.params["pg_spike_fn_base"]
        sorted_pattern_output = self.params["pg_spikes_merged_fn_base"]
        self.merge_nspike_files(merge_pattern, sorted_pattern_output, pattern)

        if (pattern==''):
            # create spiketimes file from the other three
            for pn in xrange(self.params['n_patterns']):
                mit_fn = self.params['mit_spikes_merged_fn_base'] + str(pn) + '.dat'
                gran_fn = self.params['gran_spikes_merged_fn_base'] + str(pn) + '.dat'
                pg_fn = self.params['pg_spikes_merged_fn_base'] + str(pn) + '.dat'
                ob_fn = self.params['ob_spikes_merged_fn_base'] + str(pn) + '.dat'

                rnd_nr = numpy.random.randint(0,10**8)
                tmp_fn = 'tmp_%d' % (rnd_nr)
                os.system('cat %s %s %s > %s' % (mit_fn, pg_fn, gran_fn, tmp_fn))
                os.system('sort -gk %d %s > %s' % (sort_idx, tmp_fn, ob_fn))
                os.system('rm %s' % tmp_fn)


    def merge_oc_spiketimes_files(self, pattern='', sort_idx=1):
        merge_pattern = self.params["pyr_spiketimes_fn_base"]
        sorted_pattern_output = self.params["pyr_spiketimes_merged_fn_base"]
        self.merge_spiketimes_files(merge_pattern, sorted_pattern_output, pattern)

        merge_pattern = self.params["basket_spiketimes_fn_base"]
        sorted_pattern_output = self.params["basket_spiketimes_merged_fn_base"]
        self.merge_spiketimes_files(merge_pattern, sorted_pattern_output, pattern)

        merge_pattern = self.params["rsnp_spiketimes_fn_base"]
        sorted_pattern_output = self.params["rsnp_spiketimes_merged_fn_base"]
        self.merge_spiketimes_files(merge_pattern, sorted_pattern_output, pattern)

        if (pattern != ''):
            pyr_fn = self.params['pyr_spiketimes_merged_fn_base'] + str(pattern) + '.dat'
            basket_fn = self.params['basket_spiketimes_merged_fn_base'] + str(pattern) + '.dat'
            rsnp_fn = self.params['rsnp_spiketimes_merged_fn_base'] + str(pattern) + '.dat'
            oc_fn = self.params['oc_spiketimes_merged_fn_base'] + str(pattern) + '.dat'

            rnd_nr = numpy.random.randint(0,10**8)
            tmp_fn = 'tmp_%d' % (rnd_nr)
            os.system('cat %s %s %s > %s' % (pyr_fn, rsnp_fn, basket_fn, tmp_fn))
#            print 'cat %s %s %s > %s' % (pyr_fn, rsnp_fn, basket_fn, tmp_fn)
            os.system('sort -gk %d %s > %s' % (sort_idx, tmp_fn, oc_fn))
#            print 'sort -gk 1 %s > %s' % (tmp_fn, oc_fn)
            os.system('rm %s' % tmp_fn)


        else: # merge for all patterns
            # create spiketimes file from the other three
            for pn in xrange(self.params['n_patterns']):
                pyr_fn = self.params['pyr_spiketimes_merged_fn_base'] + str(pn) + '.dat'
                basket_fn = self.params['basket_spiketimes_merged_fn_base'] + str(pn) + '.dat'
                rsnp_fn = self.params['rsnp_spiketimes_merged_fn_base'] + str(pn) + '.dat'
                oc_fn = self.params['oc_spiketimes_merged_fn_base'] + str(pn) + '.dat'

                rnd_nr = numpy.random.randint(0,10**8)
                tmp_fn = 'tmp_%d' % (rnd_nr)
                os.system('cat %s %s %s > %s' % (pyr_fn, rsnp_fn, basket_fn, tmp_fn))
                os.system('sort -gk %d %s > %s' % (sort_idx, tmp_fn, oc_fn))
                os.system('rm %s' % tmp_fn)


    def merge_oc_nspike_files(self, pattern='', sort_idx=1):
        # merge for each cell type individually
        merge_pattern = self.params["pyr_spike_fn_base"]
        sorted_pattern_output = self.params["pyr_spikes_merged_fn_base"]
        self.merge_nspike_files(merge_pattern, sorted_pattern_output, pattern)

        merge_pattern = self.params["basket_spike_fn_base"]
        sorted_pattern_output = self.params["basket_spikes_merged_fn_base"]
        self.merge_nspike_files(merge_pattern, sorted_pattern_output, pattern)

        merge_pattern = self.params["rsnp_spike_fn_base"]
        sorted_pattern_output = self.params["rsnp_spikes_merged_fn_base"]
        self.merge_nspike_files(merge_pattern, sorted_pattern_output, pattern)

        # merge output of all cell_types into one file
        if (pattern != ''):
            pyr_fn = self.params['pyr_spikes_merged_fn_base'] + str(pattern) + '.dat'
            basket_fn = self.params['basket_spikes_merged_fn_base'] + str(pattern) + '.dat'
            rsnp_fn = self.params['rsnp_spikes_merged_fn_base'] + str(pattern) + '.dat'
            oc_fn = self.params['oc_spikes_merged_fn_base'] + str(pattern) + '.dat'

            rnd_nr = numpy.random.randint(0,10**8)
            tmp_fn = 'tmp_%d' % (rnd_nr)
            os.system('cat %s %s %s > %s' % (pyr_fn, rsnp_fn, basket_fn, tmp_fn))
            os.system('sort -gk %d %s > %s' % (sort_idx, tmp_fn, oc_fn))
            os.system('rm %s' % tmp_fn)

        else:
            # create spiketimes file from the other three
            for pn in xrange(self.params['n_patterns']):
                pyr_fn = self.params['pyr_spikes_merged_fn_base'] + str(pn) + '.dat'
                basket_fn = self.params['basket_spikes_merged_fn_base'] + str(pn) + '.dat'
                rsnp_fn = self.params['rsnp_spikes_merged_fn_base'] + str(pn) + '.dat'
                oc_fn = self.params['oc_spikes_merged_fn_base'] + str(pn) + '.dat'

                rnd_nr = numpy.random.randint(0,10**8)
                tmp_fn = 'tmp_%d' % (rnd_nr)
                os.system('cat %s %s %s > %s' % (pyr_fn, rsnp_fn, basket_fn, tmp_fn))
                os.system('sort -gk %d %s > %s' % (sort_idx, tmp_fn, oc_fn))
                os.system('rm %s' % tmp_fn)

    def merge_readout_spiketimes_files(self, pattern=''):

        merge_pattern = self.params["readout_spiketimes_fn_base"]
        sorted_pattern_output = self.params["readout_spiketimes_merged_fn_base"]
        self.merge_spiketimes_files(merge_pattern, sorted_pattern_output, pattern)


    def merge_readout_nspike_files(self, pattern=''):

        merge_pattern = self.params["readout_spike_fn_base"]
        sorted_pattern_output = self.params["readout_spikes_merged_fn_base"]
        self.merge_nspike_files(merge_pattern, sorted_pattern_output, pattern)


if __name__ == '__main__':
    info_txt = \
    """
    Usage:
        python MergeSpikeFiles.py [FOLDER] [CELLTYPE] 
        or
        python MergeSpikeFiles.py [FOLDER] [CELLTYPE] [PATTERN_NUMBER]

    """
    assert (len(sys.argv) > 2), 'ERROR: folder and cell_type not given\n' + info_txt
    folder = sys.argv[1]
    cell_type = sys.argv[2]
    params_fn = os.path.abspath(folder) + '/Parameters/simulation_parameters.json'
    param_tool = simulation_parameters.parameter_storage(params_fn=params_fn)
    params = param_tool.params
    try:
        pn_max = int(sys.argv[3])
    except:
        print 'Merging all patterns'
        pn_max = params['n_patterns']

    if cell_type == 'all':
        cell_types = params['cell_types']
    else:
        cell_types = [cell_type]


    MS = MergeSpikefiles(params)
    for cell_type in cell_types:
        for pattern in xrange(pn_max):
            print 'Merging nspike file for %s pattern %d' % (cell_type, pattern)
            MS.merge_nspike_files(params['%s_spike_fn_base' % cell_type], params['%s_spikes_merged_fn_base' % cell_type], pattern)
            print 'Merging spiketimes file for %s pattern %d' % (cell_type, pattern)
            MS.merge_spiketimes_files(params['%s_spiketimes_fn_base' % cell_type], params['%s_spiketimes_merged_fn_base' % cell_type], pattern)
            
