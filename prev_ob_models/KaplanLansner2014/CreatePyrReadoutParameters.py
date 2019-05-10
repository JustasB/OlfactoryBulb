import numpy as np

class CreatePyrReadoutParameters(object):
    """
    After the run_epth_prelearning.py this class gets called from
    create_ob_oc_connectivity and sets parameters for the OC and Readout cells
    based on the results from the MDS and VQ.
    """

    def __init__(self, param_dict):
        self.params = param_dict


    def write_pyr_parameters(self):
        bias_file = self.params['ob_oc_abstract_bias_fn']
        biases = np.loadtxt(bias_file)
        assert (biases.size == self.params['n_hc'] * self.params['n_mc']), "Number of MC does not match in bias file and parameterfile"
        if (self.params['with_curr_bias']):
            cv = self.get_bias_iamp_conversion_factor()
            i_bias = biases * cv  # [uS / cm^2]
        else:
            cv = self.get_bias_conductance_conversion_factor()
            conductances = biases * cv  # [uS / cm^2]

        output_fn = self.params['pyr_params_file']
        f = open(output_fn, 'w')
        lines_to_write = ""

        if (self.params['with_curr_bias']):
            output_params = i_bias
        else:
            output_params = conductances

        for mc in xrange(self.params['n_hc'] * self.params['n_mc']):
            for pyr in xrange(self.params['n_pyr_per_mc']):
                gid = pyr + mc * self.params['n_pyr_per_mc'] + self.params['pyr_offset']
                g_or_i_bias = output_params[mc]
                lines_to_write += "%d\t%.4e\n" % (gid, g_or_i_bias)

        n_rows = self.params['n_pyr']
        n_cols = 2
        first_line = "%d\t%d\n" % (n_rows, n_cols)
        f.write(first_line)
        f.write(lines_to_write)
        f.close()


    def write_readout_parameters(self):
        bias_file = self.params['oc_readout_abstract_bias_fn']
        biases = np.loadtxt(bias_file)
        assert (biases.size == self.params['n_readout']), "Number of Readout cells does not match in bias file and parameterfile %s" % self.params['oc_readout_abstract_bias_fn']
        if (self.params['with_curr_bias']):
            cv = self.get_bias_iamp_conversion_factor()
            i_bias = biases * cv  # [uS / cm^2]
        else:
            cv = self.get_bias_conductance_conversion_factor()
            conductances = biases * cv  # [uS / cm^2]

        output_fn = self.params['readout_params_file']
        f = open(output_fn, 'w')
        lines_to_write = ""

        if (self.params['with_curr_bias']):
            output_params = i_bias
        else:
            output_params = conductances

        for readout in xrange(self.params['n_readout']):
            gid = readout + self.params['readout_offset']
            g_or_i_bias = output_params[readout]
            lines_to_write += "%d\t%.4e\n" % (gid, g_or_i_bias)

        n_rows = self.params['n_readout']
        n_cols = 2
        first_line = "%d\t%d\n" % (n_rows, n_cols)
        f.write(first_line)
        f.write(lines_to_write)
        f.close()

    def get_bias_conductance_conversion_factor(self):
        fn = self.params['ob_oc_abstract_bias_fn']
        d = np.loadtxt(fn)
        bias_max = d.min() 
        bias_conversion_factor = self.params['g_ka_pyr_max'] / bias_max
        print "bias_conversion_factor: ", bias_conversion_factor
        return bias_conversion_factor


    def get_bias_iamp_conversion_factor(self):
        fn = self.params['ob_oc_abstract_bias_fn']
        d = np.loadtxt(fn)
        bias_max = d.min() 
        bias_conversion_factor = self.params['i_bias_pyr_max'] / bias_max
        print "bias_conversion_factor: ", bias_conversion_factor
        return bias_conversion_factor



    def write_pyr_parameters_currents(self):
        bias_file = self.params['ob_oc_abstract_bias_fn']
        biases = np.loadtxt(bias_file)
        assert (biases.size == self.params['n_hc'] * self.params['n_mc']), "Number of MC does not match in bias file and parameterfile"
        cv = self.get_bias_current_conversion_factor()
        currentamps = biases * cv  # [uS / cm^2]
        output_fn = self.params['pyr_params_file']
        f = open(output_fn, 'w')
        lines_to_write = ""
        for mc in xrange(self.params['n_hc'] * self.params['n_mc']):
            for pyr in xrange(self.params['n_pyr_per_mc']):
                gid = pyr + mc * self.params['n_pyr_per_mc'] + self.params['pyr_offset']
                i_amp = currentamps[mc]
                lines_to_write += "%d\t%.4e\n" % (gid, i_amp)
        n_rows = self.params['n_pyr']
        n_cols = 2
        first_line = "%d\t%d\n" % (n_rows, n_cols)
        f.write(first_line)
        f.write(lines_to_write)
        f.close()


    def write_readout_parameters_currents(self):
        bias_file = self.params['oc_readout_abstract_bias_fn']
        biases = np.loadtxt(bias_file)
        assert (biases.size == self.params['n_readout']), "Number of Readout cells does not match in bias file and parameterfile"
        cv = self.get_bias_current_conversion_factor()
        currentamps = biases * cv  # [uS / cm^2]
        output_fn = self.params['readout_params_file']
        f = open(output_fn, 'w')
        lines_to_write = ""
        for readout in xrange(self.params['n_readout']):
            gid = readout + self.params['readout_offset']
            i_amp = currentamps[readout]
            lines_to_write += "%d\t%.4e\n" % (gid, i_amp)
        n_rows = self.params['n_readout']
        n_cols = 2
        first_line = "%d\t%d\n" % (n_rows, n_cols)
        f.write(first_line)
        f.write(lines_to_write)
        f.close()


    def get_bias_current_conversion_factor(self):
        fn = self.params['ob_oc_abstract_bias_fn']
        d = np.loadtxt(fn)
        bias_max = d.min() 
        bias_conversion_factor = self.params['i_bias_pyr_max'] / bias_max
        print "bias_conversion_factor: ", bias_conversion_factor
        return bias_conversion_factor

