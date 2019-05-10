import numpy as np
import os

class BCPNN(object):
    """
    BCPNN
    1) load a pattern, i.e. the normalized ob activity
    2) init:
        initialize the weights to some uniform value
        initialize biases to 1 / n_patterns
    3) for all patterns:
        calculate the post-synaptic activities: s_j
    4) get new p_j, p_ij values -> recalculate w_ij, b_j
    """

    def __init__(self, n_hc_in, n_mc_in, n_hc_out, n_mc_out, n_patterns, params, change_thresh=1e-10):

        self.params = params
        self.n_hc_in = n_hc_in
        self.n_mc_in = n_mc_in
        self.n_hc_out = n_hc_out
        self.n_mc_out = n_mc_out
        self.n_patterns = n_patterns
        self.gain_parameter = 1.
        self.mc_mc_mask = np.ones((self.n_hc_in * self.n_mc_in, self.n_hc_out * self.n_mc_out))
        n_pre = self.n_hc_in * self.n_mc_in
        n_post = self.n_hc_out * self.n_mc_out
        self.w_ij = np.ones((n_pre, n_post)) * (1./ self.n_patterns) # initialize with ones, in case there is no restricting mc_hc_mask
        self.post_activity = np.zeros((n_patterns, n_post))
        self.no_big_change = False
        self.iteration = 0
        self.p_i = np.zeros(n_pre)
        self.p_j = np.zeros(n_post)
        self.p_ij = np.zeros((n_pre, n_post))
        self.bias = np.ones(n_post) * np.log((1./ self.n_patterns)**2)
        self.change_thresh = change_thresh
        self.pre_loaded = False
        self.post_loaded = False


    def initialize(self):
        """
        This function calculates weights and biases, and probabilities from the loaded input and output activity.
        So, before calling you must call:

            bcpnn.load_input_activity(ob_activity_fn)
            bcpnn.load_output_activity(binary_oc_activation_fn)
            and optionally a mask:
            bcpnn.load_mc_hc_mask(w_ij_mit_hc, silent_units_fn=params['silent_mit_fn'])

        """
#        n_pre = self.n_hc_in * self.n_mc_in
#        n_post = self.n_hc_out * self.n_mc_out
#        self.p_i = np.zeros(n_pre)
#        self.p_j = np.zeros(n_post)
#        self.p_ij = np.zeros((n_pre, n_post))
#        self.bias = np.ones(n_post) * np.log((1./ self.n_patterns)**2)

        # show all patterns once and activate units in the output layer and apply WTA to the post activity
#        for pn in xrange(self.n_patterns):
#            pre_activity = self.input_activity[pn, :]
#            for post in xrange(n_post): # mc index
#                in_j = 0.
#                for pre in xrange(n_pre):
#                    in_j += (self.w_ij[pre, post] * pre_activity[pre])
#                self.post_activity[pn, post] = in_j

#        print "Calculating probabilities: ", self.iteration
#        self.calculate_probabilities()
#        print "Calculating weights and bias: ", self.iteration
#        self.calculate_weights_and_bias()

        debug_fn_1 = self.params['bcpnn_folder'] + "/weights_after_init_wij_mc_hc.dat"
        debug_fn_2 = self.params['bcpnn_folder'] + "/bias_after_init.dat"
        debug_fn_3 = self.params['bcpnn_folder'] + "/p_ij_after_init.dat"
        debug_fn_4 = self.params['bcpnn_folder'] + "/post_activity_after_init.dat"
        debug_fn_5 = self.params['bcpnn_folder'] + "/pi_after_init.dat"
        debug_fn_6 = self.params['bcpnn_folder'] + "/pj_after_init.dat"
        debug_fn_7 = self.params['bcpnn_folder'] + "/input_activity_after_init.dat"

        np.savetxt(debug_fn_1, self.w_ij)
        np.savetxt(debug_fn_2, self.bias)
        np.savetxt(debug_fn_3, self.p_ij)
        np.savetxt(debug_fn_4, self.post_activity)
        np.savetxt(debug_fn_5, self.p_i)
        np.savetxt(debug_fn_6, self.p_j)
        np.savetxt(debug_fn_7, self.input_activity)
#        exit(1)


    def load_input_activity(self, fn_or_array):
        """
        Load a binary activation pattern before the training.
        For OC-OC:
        This binary activation was derived from MIT cell activations projecting to HC and a VQ in this masked MIT-activity space.
        See MDSVQ.create_mitral_response_space(input_fn, activity_fn, output_fn, remove_silent_cells_fn=False) for the code
        """
        if (type(fn_or_array) == type('')):
            print 'BCPNN.load_input_activity from: ', fn_or_array
            self.input_activity = np.loadtxt(fn_or_array)

        elif (type(fn_or_array) == type(np.array([]))):
            self.input_activity = fn_or_array

        self.pre_loaded = True

    def load_output_activity(self, fn_or_array):
        """
        Load a the activation pattern before the training.
        For OB-OC:
        This binary activation was derived from MIT cell activations projecting to HC and a VQ in this masked MIT-activity space.
        See MDSVQ.create_mitral_response_space(input_fn, activity_fn, output_fn, remove_silent_cells_fn=False) for the code
        """
        if (type(fn_or_array) == type('')):
            print 'BCPNN.load_output_activity from: ', fn_or_array
            self.post_activity = np.loadtxt(fn_or_array)
        elif (type(fn_or_array) == type(np.array([]))):
            self.post_activity = fn_or_array

        self.post_loaded = True


    def load_mc_hc_mask(self, mc_hc_conn_fn, silent_units_fn=False):
        """
        Load the result of the VQ, e.g. the result of the VQ in the mitral cell response space.
        Thus, the number of rows might be smaller than the total number of minicolumns in the network (because silent input units have been removed for the VQ).
        """
        mc_hc_mask = np.loadtxt(mc_hc_conn_fn)
        silent_units = []
        if (silent_units_fn):
            silent_units = np.loadtxt(self.params['silent_mit_fn'])
            
        self.mc_mc_mask = np.zeros((self.n_hc_in * self.n_mc_in, self.n_hc_out * self.n_mc_out))
        for src_mc in xrange(self.n_mc_in * self.n_hc_in):
            tgt_hcs = mc_hc_mask[src_mc, :].nonzero()[0]
            for tgt_hc in tgt_hcs:
                mc1 = tgt_hc * self.n_mc_out
                mc2 = (tgt_hc + 1) * self.n_mc_out
                self.mc_mc_mask[src_mc, mc1:mc2] = np.ones(self.n_mc_out)



    def train_network(self, post_act_fn='', weights_fn='', bias_fn='', debug=False):
        assert (self.pre_loaded and self.post_loaded), 'You need to load input and output probabilities first'

        print "Calculating probabilities: ", self.iteration
        self.calculate_probabilities()

        print "Calculating weights and bias: ", self.iteration
        self.calculate_weights_and_bias()

        for pn in xrange(self.n_patterns):
            self.calculate_post_activity(pn)

#        if post_act_fn != '':
#            output_fn = post_act_fn.rsplit('.dat')[0] + '%d.dat' % (self.iteration)
#            print "Saving data for iteration %d:  " % (self.iteration), output_fn
#            np.savetxt(output_fn, self.post_activity, delimiter='\t')
#        if weights_fn != '':
#            output_fn = weights_fn.rsplit('.dat')[0] + '%d.dat' % (self.iteration)
#            print ".... ", output_fn
#            np.savetxt(output_fn, self.w_ij, delimiter='\t')
#        if bias_fn != '':
#            output_fn = bias_fn.rsplit('.dat')[0] + '%d.dat' % (self.iteration)
#            print ".... ", output_fn
#            np.savetxt(output_fn, self.bias, delimiter='\t')

        if debug:
            debug_fn_1 = "test_bcpnn_ob_oc/weights_after_init_wij_mc_hc_%d.dat" % (self.iteration)
            debug_fn_2 = "test_bcpnn_ob_oc/bias_after_init_%d.dat" % (self.iteration)
            debug_fn_3 = "test_bcpnn_ob_oc/p_ij_after_init_%d.dat" % (self.iteration)
            debug_fn_4 = "test_bcpnn_ob_oc/post_activity_after_init_%d.dat" % (self.iteration)
            debug_fn_5 = "test_bcpnn_ob_oc/pi_after_init_%d.dat" % (self.iteration)
            debug_fn_6 = "test_bcpnn_ob_oc/pj_after_init_%d.dat" % (self.iteration)
            debug_fn_7 = "test_bcpnn_ob_oc/input_activity_after_init_%d.dat" % (self.iteration)

            np.savetxt(debug_fn_1, self.w_ij)
            np.savetxt(debug_fn_2, self.bias)
            np.savetxt(debug_fn_3, self.p_ij)
            np.savetxt(debug_fn_4, self.post_activity)
            np.savetxt(debug_fn_5, self.p_i)
            np.savetxt(debug_fn_6, self.p_j)
            np.savetxt(debug_fn_7, self.input_activity)

#            folder = os.path.abspath(record_to_folder) + '/'
#            output_fn = folder + 'activity_after_%d.dat' % (self.iteration)
#            np.savetxt(output_fn, self.post_activity, delimiter='\t')
#            output_fn = folder + 'weights_after_%d.dat' % (self.iteration)
#            np.savetxt(output_fn, self.w_ij, delimiter='\t')
#            output_fn = folder + 'bias_after_%d.dat' % (self.iteration)
#            np.savetxt(output_fn, self.bias, delimiter='\t')
        self.iteration += 1


    def write_to_files(self, post_act_fn, weights_fn, bias_fn):

        output_fn = post_act_fn
        print "Saving data to:  ", output_fn
        np.savetxt(output_fn, self.post_activity, delimiter='\t')
        output_fn = weights_fn
        print "Saving data to:  " , output_fn
        idx_0 = (self.w_ij == 0).nonzero()
        np.savetxt(output_fn, self.w_ij, delimiter='\t')
        output_fn = bias_fn
        print "Saving data to:  " , output_fn
        np.savetxt(output_fn, self.bias, delimiter='\t')


    def calculate_post_activity(self, pn):

        n_pre = self.n_hc_in * self.n_mc_in
        n_post = self.n_hc_out * self.n_mc_out
        pre_activity = self.input_activity[pn, :]
        output_after_softmax = np.zeros(n_post)

        for post in xrange(n_post): # mc index
            in_j = 0.
            s_j = 0. # the support values
            in_j = np.sum(self.w_ij[:, post] * pre_activity * self.mc_mc_mask[:, post])
            if in_j != 0:
                s_j = self.bias[post] + in_j # if weight are already in log domain
#                s_j = self.bias[post] + np.log(in_j)
                self.post_activity[pn, post] = np.exp(s_j)
            else:
                self.post_activity[pn, post] = 0

            # NEW
            for pre_hc_0 in xrange(self.n_hc_in):
                pre_mc_0 = pre_hc_0 * self.n_mc_in
                pre_mc_1 = (pre_hc_0 + 1) * self.n_mc_in
                from_hc = np.sum(self.w_ij[pre_mc_0:pre_mc_1, post] * pre_activity[pre_mc_0:pre_mc_1])
                # what if from_hc == 0 ????


#            print 'debug pn %d post %d in_j %.5e' % (pn, post, in_j.sum())
#            if in_j.sum() != 0:
#                s_j = self.bias[post] + np.log(in_j.sum())
#            else:
#                s_j = self.bias[post]

#            for pre in xrange(n_pre):
#                in_j += (np.exp(self.w_ij[pre, post]) * pre_activity[pre]) * self.mc_mc_mask[pre, post ] # only allow input from masked sources (mc_mc_mask = mask)
#                in_j += np.log(self.w_ij[pre, post] * pre_activity[pre]) * self.mc_mc_mask[pre, post ] # only allow input from masked sources (mc_mc_mask = mask)

#            self.post_activity[pn, post] = s_j # before

        # normalize the post-activity:
        # apply softmax function for each hypercolumn to get the output
        for hc in xrange(self.n_hc_out):
            mc1 = hc * self.n_mc_out
            mc2 = (hc + 1) * self.n_mc_out
            summed_activity = self.post_activity[pn, mc1:mc2].sum()
            if (summed_activity > 1): # normalize
                self.post_activity[pn, mc1:mc2] /= summed_activity
#            if (summed_activity > 0): # this can lead to inf in np.exp(summed_activity)
#                print 'DEbug', np.exp(summed_activity)
#                self.post_activity[pn, mc1:mc2] = np.exp(self.post_activity[pn, mc1:mc2]) / np.exp(summed_activity)
#                    self.post_activity[pn, mc] /= summed_activity




    def calculate_probabilities(self):
        n_pre = self.n_hc_in * self.n_mc_in
        n_post = self.n_hc_out * self.n_mc_out

        # p_j
        print 'debug', self.n_hc_out, self.n_mc_out
        for post in xrange(n_post):
            self.p_j[post] = self.post_activity[:, post].sum() * (1./ self.n_patterns)
        # p_i
        for pre in xrange(n_pre):
            self.p_i[pre] = self.input_activity[:, pre].sum() * (1./ self.n_patterns)

        # p_ij
        for post in xrange(n_post):
            for pre in xrange(n_pre):
                self.p_ij[pre, post] = 0.
                self.p_ij[pre, post] = np.sum(self.input_activity[:, pre] * self.post_activity[:, post]) / float(self.n_patterns)
#                self.p_ij[pre, post] /= float(self.n_patterns)
#                for pn in xrange(self.n_patterns):

        assert (self.p_i.max() <= 1.0)
        assert (self.p_i.min() >= 0.0)
        assert (self.p_j.max() <= 1.0), "P_j max = %.2f argmax = %d" % (self.p_j.max(), self.p_j.argmax())
        assert (self.p_j.min() >= 0.0)
        assert (self.p_ij.max() <= 1.0)
        assert (self.p_ij.min() >= 0.0)


    def calculate_weights_and_bias(self):

        n_pre = self.n_hc_in * self.n_mc_in
        n_post = self.n_hc_out * self.n_mc_out
        self.weights_old = self.w_ij.copy()
        print "Calculating weights and bias ... "

        for post in xrange(n_post):
            # bias
            if (self.p_j[post] == 0.):
                self.bias[post] = np.log((1. / self.n_patterns)**2)
            else:
                self.bias[post] = np.log(self.p_j[post])

            # weights
            hc = post / self.n_mc_out
            for pre in xrange(n_pre):
                if ((self.p_i[pre] == 0.) or (self.p_j[post] == 0.)):
                    self.w_ij[pre, post] = 0.
                elif (self.p_ij[pre, post] == 0.):
                    self.w_ij[pre, post] = (1. / self.n_patterns) * self.mc_mc_mask[pre, post]
                else:
                    self.w_ij[pre, post] = np.log(self.p_ij[pre, post] / (self.p_i[pre] * self.p_j[post])) * self.mc_mc_mask[pre, post]
#                    self.w_ij[pre, post] = (self.p_ij[pre, post] / (self.p_i[pre] * self.p_j[post])) * self.mc_mc_mask[pre, post]



    def silence_mit(self, silent_units_fn):
        """
        Connections involving input units with ids stored in silent_units_fn, are set to zero.
        """
        silent_units = np.loadtxt(silent_units_fn)
        for i in xrange(len(silent_units)):
            self.w_ij[silent_units[i], :] = np.zeros(self.w_ij[0, :].size)

        # apply the mc_hc_mask to all weights
        for i in xrange(self.n_mc_in * self.n_hc_in): # src
            for tgt_mc in xrange(self.n_hc_out * self.n_mc_out):
                if (self.mc_mc_mask[i, tgt_mc] == 0):
                    if (self.w_ij[i, tgt_mc] != 0):
                        print "ERROR: Something went wrong with the input-output connection mask, check load_mc_hc_mask..."
                        self.w_ij[i, j] = 0



         
    def testing(self, input_fn_or_array, weights_fn_or_array, bias_fn_or_array, output_fn=None):

        # load from a filename or array
        if (type(input_fn_or_array) == type('')):
            try: # for OB activity
                input_activity = np.loadtxt(input_fn_or_array)
            except:  # for OC activity
                input_activity = np.loadtxt(input_fn_or_array, delimiter=',')
        elif (type(input_fn_or_array) == type(np.array([]))):
            input_activity = input_fn_or_array

        if (type(weights_fn_or_array) == type('')):
            try: # for OB activity
                weights = np.loadtxt(weights_fn_or_array)
            except:  # for OC activity
                weights = np.loadtxt(weights_fn_or_array, delimiter=',')
        elif (type(weights_fn_or_array) == type(np.array([]))):
            weights = weights_fn_or_array

        if (type(bias_fn_or_array) == type('')):
            try: # for OB activity
                bias_test = np.loadtxt(bias_fn_or_array)
            except:  # for OC activity
                bias_test = np.loadtxt(bias_fn_or_array, delimiter=',')
        elif (type(bias_fn_or_array) == type(np.array([]))):
            bias_test = bias_fn_or_array
        
        # calculate the post activity with the test input
        n_pre = self.n_hc_in * self.n_mc_in
        n_post = self.n_hc_out * self.n_mc_out
#        output_after_softmax = np.zeros((self.n_patterns, n_post))
        post_activity_test = np.zeros((self.n_patterns, n_post))
        wta_activity = np.zeros((self.n_patterns, n_post))

        for pn in xrange(self.n_patterns):

            for post in xrange(n_post): # mc index
                in_j = 0.
                s_j = 0. # the support values
                in_j = np.sum(weights[:, post] * input_activity[pn, :] * self.mc_mc_mask[:, post ])
#                in_j = (self.w_ij[:, post] * input_activity[pn, :] * self.mc_mc_mask[:, post ]).sum() # only allow input from masked sources (mc_mc_mask = mask)
#                for pre in xrange(n_pre):
#                    in_j += (np.exp(self.w_ij[pre, post]) * input_activity[pn, pre]) * self.mc_mc_mask[pre, post ] # only allow input from masked sources (mc_mc_mask = mask)
#                in_j = ((np.exp(self.w_ij[:, post]) * input_activity[pn, :]) * self.mc_mc_mask[:, post ]).sum() # only allow input from masked sources (mc_mc_mask = mask)
#                    in_j += (self.w_ij[pre, post] * input_activity[pn, pre]) * self.mc_mc_mask[pre, post ] # only allow input from masked sources (mc_mc_mask = mask)

#                s_j = np.exp(self.bias[post]) + in_j
#                print 'debug in_j %d = %.6f, bias=%.2f' % (post, in_j, bias_test[post]), 
                s_j = bias_test[post] + np.log(in_j)
                post_activity_test[pn, post] = np.exp(s_j)
#                print ' post_activity_test[%d, %d] = %.6f' % (pn, post, post_activity_test[pn, post])
#                post_activity_test[pn, post] = s_j

            post_activity_test[pn, :] /= post_activity_test[pn, :].sum()
            wta_activity[pn, np.argmax(post_activity_test[pn, :])] = 1

            # apply softmax function for each hypercolumn to get the output
#            for hc in xrange(self.n_hc_out):
#                mc1 = hc * self.n_mc_out
#                mc2 = (hc + 1) * self.n_mc_out
#                summed_activity = post_activity_test[pn, mc1:mc2].sum()
#            output_after_softmax[pn, :] = post_activity_test[pn, :] / post
#            wta_activity[pn, np.argmax(post_activity_test[pn, :])] = 1
#                if (summed_activity > 1):
#                    output_after_softmax[pn, mc1:mc2] = post_activity_test[pn, mc1:mc2] / summed_activity
             
    
        if (output_fn != None):
#            softmax_output_fn = output_fn.rsplit('.dat')[0] + '_softmax.dat'
#            np.savetxt(softmax_output_fn, output_after_softmax)
            np.savetxt(output_fn, post_activity_test)

            wta_output_fn = output_fn.rsplit('.dat')[0] + '_wta.dat'
            np.savetxt(wta_output_fn, wta_activity)
            print "Saving activity during testing to:\n", output_fn, '\n', wta_output_fn

        list_of_winners = np.ones(self.n_patterns) # list of most active readout cells
        list_of_winners *= -1
        correct_patterns = []
        for pn in xrange(self.n_patterns):
            most_active_unit_test = post_activity_test[pn, :].argmax()
            most_active_unit_train = self.post_activity[pn, :].argmax()
            if (most_active_unit_test != most_active_unit_train):
                pass
#                print "Pattern: %d Activity_test[man_test=%d]=%.2e, Activity_train[man_test=%d]=%.2e, activity_test[man_train=%d]=%.2e activity_train[man_train=%d]=%.2e" \
#                        % (pn, most_active_unit_test, post_activity_test[pn, most_active_unit_test], most_active_unit_test, self.post_activity[pn, most_active_unit_test], \
#                        most_active_unit_train, post_activity_test[pn, most_active_unit_train], most_active_unit_train, self.post_activity[pn, most_active_unit_train])

            else:
                correct_patterns.append(pn)
            list_of_winners[pn] = most_active_unit_test

        winner_readouts = np.unique(list_of_winners)
        wrong_patterns = set(range(self.n_patterns)).difference(set(winner_readouts))
        print "List of winner readouts:", np.unique(list_of_winners)
        print "Number of different winner readouts:", np.unique(list_of_winners).size
        print "Wrong patterns", wrong_patterns
        print "Correct patterns:", len(correct_patterns), correct_patterns
