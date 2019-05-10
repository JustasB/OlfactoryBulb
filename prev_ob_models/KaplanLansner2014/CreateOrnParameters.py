import numpy as np
import numpy.random
import matplotlib.mlab as mlab
import random

class CreateOrnParameters(object):

    def __init__(self, param_dict):

        self.params = param_dict
        self.orn_params_fn_base = self.params['orn_params_fn_base']
        self.n_or = self.params['n_or']
        self.n_gor = self.params['n_gor']
        self.n_orn_x = self.params['n_orn_x']
        self.n_patterns = self.params['n_patterns']
        self.n_orn_y = self.params['n_orn_y']
        self.n_orn = self.params['n_orn']
        self.num_params = 9 + 1 # + 1 for the id

        self.param_list = []
        self.min_values = []
        self.max_values = []

        self.gor_values = np.zeros(self.n_orn_x).tolist()

        self.vgor = np.zeros(self.n_orn)
        self.voor = np.zeros(self.n_orn)
        self.vc_Kd= np.zeros(self.n_orn)
        self.vgna = np.zeros(self.n_orn)
        self.vgk = np.zeros(self.n_orn)
        self.vgkcag = np.zeros(self.n_orn)
        self.vgcal = np.zeros(self.n_orn)
        self.vgleak= np.zeros(self.n_orn)
        self.vtau_cadec = np.zeros(self.n_orn) # time constant for calcium removal

        self.set_gor_values()
        self.current_param_values = []
        # a list of lists containing for each row the following parameters:
                #"gna","gk","gkcag","gcal","gleak_orn", "tau_cadec"]

    def create_pattern_completion_test(self, activation_matrix_fn):
        """
        Load the given activation matrix and sample only half of the activation elements 
        to create the simpler test patterns.
        """
        activation_matrix_complex = np.loadtxt(activation_matrix_fn)
        activation_matrix = np.zeros((self.params['n_patterns'], self.params['n_or']))

        # take n_active_OR

        for pn in xrange(self.params['n_patterns']):
            all_activated_ORs = activation_matrix_complex[pn, :].nonzero()[0]
            n_OR_to_stay_active = int(round(all_activated_ORs.size * self.params['frac_ORs_incomplete_patterns']))
            remaining_activated_ORs = random.sample(all_activated_ORs, n_OR_to_stay_active)
            for OR in remaining_activated_ORs:
                activation_matrix[pn, OR] = activation_matrix_complex[pn, OR]

        output_fn = self.params['params_folder'] + '/activation_matrix_for_completion_test.dat'
        print "Activation matrix fn:", output_fn
        np.savetxt(output_fn, activation_matrix)
        for pn in xrange(self.params["n_patterns"]):
            # calculate the number of receptors activated by the current odorant or pattern
            for OR in xrange(self.params['n_or']): # equal to a row in the orn array
                # all ORNs within one row share the same Kd and concentration input value (because same OR is expressed by the these ORNs)
                activation = activation_matrix[pn, OR]# * self.params['scale_factor_activation_matrix']
                for x in xrange(self.params['n_orn_x']):
                    orn_id = OR * self.params['n_orn_x'] + x
                    self.voor[orn_id] = activation #* conc
            self.set_gor_values()
            self.set_gna_values()
            self.set_gk_values()
            self.set_gkcag_values()
            self.set_gcal_values()
            self.set_gleak_values()
            self.set_tau_cadec_values()
            output_fn = self.orn_params_fn_base + "%d.dat" % (pn)
            print "Writin orn parameters for pattern \t... %d / %d to file: %s" % (pn, self.n_patterns, output_fn)
            self.write_params_to_file(output_fn)


    def create_pattern_completion_training(self, n_active_OR):
        """
        Create a number of 'complex' patterns, in which 
        are activated.
        """
       
        activation_matrix = self.create_activation_matrix_for_pattern_completion_training(n_active_OR)
        for pn in xrange(self.params["n_patterns"]):
            # calculate the number of receptors activated by the current odorant or pattern
            for OR in xrange(self.params['n_or']): # equal to a row in the orn array
                # all ORNs within one row share the same Kd and concentration input value (because same OR is expressed by the these ORNs)
                activation = activation_matrix[pn, OR]# * self.params['scale_factor_activation_matrix']
                for x in xrange(self.params['n_orn_x']):
                    orn_id = OR * self.params['n_orn_x'] + x
                    self.voor[orn_id] = activation #* conc

            self.set_gor_values()
            self.set_gna_values()
            self.set_gk_values()
            self.set_gkcag_values()
            self.set_gcal_values()
            self.set_gleak_values()
            self.set_tau_cadec_values()
            output_fn = self.orn_params_fn_base + "%d.dat" % (pn)
            print "Writin orn parameters for pattern \t... %d / %d to file: %s" % (pn, self.n_patterns, output_fn)
            self.write_params_to_file(output_fn)
        return 1


    def create_activation_matrix_for_pattern_completion_training(self, n_active_OR):
        """
        Create random patterns with a fixed number of activated OR per pattern
        --> for pattern completion :
        n_active_OR = int(round(params['n_or'] * params['frac_max_active_OR']))
        --> for rivalry:
        n_active_OR = int(round(params['n_or'] * params['frac_min_active_OR']))
        """
        activation_matrix = np.zeros((self.params['n_patterns'], self.params['n_or']))
        p1, p2, p3 = self.set_odorant_distribution_params()

        np.random.seed(self.params['seed_activation_matrix'])
        random.seed(self.params['seed_activation_matrix'])
        distances = np.zeros((self.params['n_patterns'], self.params['n_or']))
        for pn in xrange(self.params['n_patterns']):
            activated_ORs = random.sample(range(self.params['n_or']), n_active_OR)
            while np.unique(activated_ORs).size != n_active_OR:
                activated_ORs = random.sample(range(self.params['n_or']), n_active_OR)
            for OR in activated_ORs:
                # draw a random distance from the fitted distance distribution
                dist = self.odorant_odor_distance_distribution((p1, p2, p3))
                distances[pn, OR] = dist
                affinity = np.exp(-(dist)**2 / self.params['distance_affinity_transformation_parameter_exp'])
                if affinity > 1.:
                    affinity = 1.
                if affinity < 0.:
                    affinity = 0.
                activation_matrix[pn, OR] = affinity

        if self.params['OR_activation_normalization']:
            for pn in xrange(self.params['n_patterns']):
                # this normalization increases the likelihood of having ORs with a high affinity to 
                # each given pattern --> receptors have specialized to odorants (?)
                activation_matrix[pn, :] /= activation_matrix[pn, :].max() 
                if not activation_matrix[pn, :].sum() == 0:
                    activation_matrix[pn, :] /= activation_matrix[pn, :].sum() 
                else:
                    print '\n\tWARNING: activation_matrix.sum for pattern %d == 0\nWill now quit, because that makes no sense' % (pn)
                    exit(1)
        print 'Activation matrix sum:', activation_matrix.sum()
        print "Activation matrix fn:", self.params['activation_matrix_fn']
        np.savetxt(self.params['activation_matrix_fn'], activation_matrix)
        return activation_matrix


    def create_pattern_rivalry_test(self, activation_matrix_fn):
        """
        Load the given activation matrix from the 'training' run and build composite patterns built from the training patterns.

        """
        simple_activation_matrix = np.loadtxt(activation_matrix_fn)
        activation_matrix = np.zeros((self.params['n_patterns'], self.params['n_or']))

        dgn = self.params['OR_affinity_noise']
        for pn in xrange(self.params['n_patterns'] - 1):
            activation_matrix[pn, :] = simple_activation_matrix[pn, :] + simple_activation_matrix[pn + 1, :]
            for OR in xrange(self.params['n_or']):
                activation_matrix[pn, OR] += np.random.uniform(-dgn, dgn)
                if (activation_matrix[pn, OR] > 1):
                    activation_matrix[pn, OR] = 1
                elif (activation_matrix[pn, OR] < 0):
                    activation_matrix[pn, OR] = 0.

        activation_matrix[self.params['n_patterns'] - 1, :] = simple_activation_matrix[self.params['n_patterns'] - 1, :] + simple_activation_matrix[self.params['n_patterns'] - 2, :]

        output_fn = self.params['params_folder'] + '/activation_matrix_for_rivalrytest.dat'
        print "Activation matrix fn:", output_fn
        np.savetxt(output_fn, activation_matrix)
        for pn in xrange(self.params["n_patterns"]):
            # calculate the number of receptors activated by the current odorant or pattern
            for OR in xrange(self.params['n_or']): # equal to a row in the orn array
                # all ORNs within one row share the same Kd and concentration input value (because same OR is expressed by the these ORNs)
                activation = activation_matrix[pn, OR]# * self.params['scale_factor_activation_matrix']
                for x in xrange(self.params['n_orn_x']):
                    orn_id = OR * self.params['n_orn_x'] + x
                    self.voor[orn_id] = activation #* conc
            self.set_gor_values()
            self.set_gna_values()
            self.set_gk_values()
            self.set_gkcag_values()
            self.set_gcal_values()
            self.set_gleak_values()
            self.set_tau_cadec_values()
            output_fn = self.orn_params_fn_base + "%d.dat" % (pn)
            print "Writin orn parameters for pattern \t... %d / %d to file: %s" % (pn, self.n_patterns, output_fn)
            self.write_params_to_file(output_fn)



    def create_pattern_rivalry_morphing_test(self, activation_matrix_fn):
        """
        Load the given activation matrix from the 'training' run and build composite patterns built from the training patterns.
        """
        simple_activation_matrix = np.loadtxt(activation_matrix_fn)
        n_simple_patterns = simple_activation_matrix.shape[0]
        activation_matrix = np.zeros((self.params['n_patterns'], self.params['n_or']))

        pattern_idx = 0
        random.seed(self.params['seed_activation_matrix'])
        for pn in xrange(self.params['n_patterns_test_rivalry']): # ~ 10 patterns
            pn0 = pn 
            pn1 = (pn + 1) % n_simple_patterns
#            pn0 = (pn * 2) % n_simple_patterns
#            pn1 = (pn * 2 + 1) % n_simple_patterns

            activated_in_simple_0 = simple_activation_matrix[pn0, :].nonzero()[0]
            activated_in_simple_1 = simple_activation_matrix[pn1, :].nonzero()[0]
            activated_from_pn0 = set([]) # start as if the pattern is currently empty and add ORs for each stage
            activated_from_pn1 = set(activated_in_simple_1) # start as if it is the original pattern
            ORs_active_in_simple_pattern_not_yet_activated_pn0 = list(activated_in_simple_0)

            for i_, stage in enumerate(self.params['rivalry_morph_stages']):
                tgt_n_from_pn0 = int(round(stage * activated_in_simple_0.size))
                tgt_n_from_pn1 = int(round((1. - stage) * activated_in_simple_1.size))
                n_to_be_added_to_pn0 =  tgt_n_from_pn0 - len(activated_from_pn0)
                while (len(activated_from_pn0) < tgt_n_from_pn0):
                    new_OR = random.choice(ORs_active_in_simple_pattern_not_yet_activated_pn0)
                    activated_from_pn0.update([new_OR])
                assert (tgt_n_from_pn0 == len(activated_from_pn0))
                
                n_to_be_removed_from_pn1 = len(activated_from_pn1) - tgt_n_from_pn1
                for j_ in xrange(n_to_be_removed_from_pn1):
#                    activated_from_pn1.remove(random.choice(activated_from_pn1))
                    activated_from_pn1.pop()

                print 'Combining patterns %d (%d) & %d (%d) (stage %.1f) take_over_0 (%d)' % (pn0, activated_in_simple_0.size, pn1, activated_in_simple_1.size, stage, len(activated_from_pn0)), activated_from_pn0,
                print '\t take_over_1 (%d)' % (len(activated_from_pn1)), activated_from_pn1
                activation_matrix[pattern_idx, list(activated_from_pn0)] = simple_activation_matrix[pn0, list(activated_from_pn0)]
                activation_matrix[pattern_idx, list(activated_from_pn1)] += simple_activation_matrix[pn1, list(activated_from_pn1)]

                for OR in xrange(self.params['n_or']):
                    if (activation_matrix[pattern_idx, OR] > 1):
                        activation_matrix[pattern_idx, OR] = 1
                    elif (activation_matrix[pattern_idx, OR] < 0):
                        activation_matrix[pattern_idx, OR] = 0.

                pattern_idx += 1

        output_fn = self.params['params_folder'] + '/activation_matrix_for_rivalrytest.dat'
        print "Activation matrix fn:", output_fn
        np.savetxt(output_fn, activation_matrix)
        for pn in xrange(self.params["n_patterns"]):
            # calculate the number of receptors activated by the current odorant or pattern
            for OR in xrange(self.params['n_or']): # equal to a row in the orn array
                # all ORNs within one row share the same Kd and concentration input value (because same OR is expressed by the these ORNs)
                activation = activation_matrix[pn, OR]# * self.params['scale_factor_activation_matrix']
                for x in xrange(self.params['n_orn_x']):
                    orn_id = OR * self.params['n_orn_x'] + x
                    self.voor[orn_id] = activation #* conc
            self.set_gor_values()
            self.set_gna_values()
            self.set_gk_values()
            self.set_gkcag_values()
            self.set_gcal_values()
            self.set_gleak_values()
            self.set_tau_cadec_values()
            output_fn = self.orn_params_fn_base + "%d.dat" % (pn)
            print "Writin orn parameters for pattern \t... %d / %d to file: %s" % (pn, self.n_patterns, output_fn)
            self.write_params_to_file(output_fn)





    def create_conc_invariance_patterns(self, given_activation_matrix=None):
        """
        orthogonal patterns means a single odorant is presented to the system
        """
        if not given_activation_matrix:
            activation_matrix = self.create_single_odorant_activation_matrix()
        else:
            activation_matrix = np.loadtxt(given_activation_matrix)

        activation_matrix = self.modify_patterns_for_conc_inv(activation_matrix)

        for pn in xrange(self.params["n_patterns"]):
            # calculate the number of receptors activated by the current odorant or pattern
            for OR in xrange(self.params['n_or']): # equal to a row in the orn array
                # all ORNs within one row share the same Kd and concentration input value (because same OR is expressed by the these ORNs)
                activation = activation_matrix[pn, OR]# * self.params['scale_factor_activation_matrix']
                for x in xrange(self.params['n_orn_x']):
                    orn_id = OR * self.params['n_orn_x'] + x
                    self.voor[orn_id] = activation #* conc

            self.set_gor_values()
            self.set_gna_values()
            self.set_gk_values()
            self.set_gkcag_values()
            self.set_gcal_values()
            self.set_gleak_values()
            self.set_tau_cadec_values()
            output_fn = self.orn_params_fn_base + "%d.dat" % (pn)
            print "Writin orn parameters for pattern \t... %d / %d to file: %s" % (pn, self.n_patterns, output_fn)
            self.write_params_to_file(output_fn)
        return 1




    def create_single_odorant_patterns(self, given_activation_matrix=None):
        """
        orthogonal patterns means a single odorant is presented to the system
        """

        if not given_activation_matrix:
            activation_matrix = self.create_single_odorant_activation_matrix()
        else:
            activation_matrix = np.loadtxt(given_activation_matrix)

        if self.params['OR_affinity_noise'] != 0.0:
            print 'Adding noise to the activation matrix ...'
            activation_matrix = self.add_pattern_noise(activation_matrix)

        for pn in xrange(self.params["n_patterns"]):
            # calculate the number of receptors activated by the current odorant or pattern
            for OR in xrange(self.params['n_or']): # equal to a row in the orn array
                # all ORNs within one row share the same Kd and concentration input value (because same OR is expressed by the these ORNs)
                activation = activation_matrix[pn, OR]# * self.params['scale_factor_activation_matrix']
                for x in xrange(self.params['n_orn_x']):
                    orn_id = OR * self.params['n_orn_x'] + x
                    self.voor[orn_id] = activation #* conc

            self.set_gor_values()
            self.set_gna_values()
            self.set_gk_values()
            self.set_gkcag_values()
            self.set_gcal_values()
            self.set_gleak_values()
            self.set_tau_cadec_values()
            output_fn = self.orn_params_fn_base + "%d.dat" % (pn)
            print "Writin orn parameters for pattern \t... %d / %d to file: %s" % (pn, self.n_patterns, output_fn)
            self.write_params_to_file(output_fn)
        return 1


    def create_single_odorant_activation_matrix(self):
        """
        create_single_odorant_activation_matrix 
        """
        activation_matrix = np.zeros((self.params['n_patterns'], self.params['n_or']))
        p1, p2, p3 = self.set_odorant_distribution_params()
        np.random.seed(self.params['seed_activation_matrix'])
        random.seed(self.params['seed_activation_matrix'])
        distances = np.zeros((self.params['n_patterns'], self.params['n_or']))
        n_min_active_OR = int(round(self.params['n_or'] * self.params['frac_min_active_OR']))
        n_max_active_OR = int(round(self.params['n_or'] * self.params['frac_max_active_OR']))
        n_above_thresh = np.zeros(self.params['n_patterns']) # number of glomeruli expected to get an activation larger than thresh
        thresh = 0.01
        for pn in xrange(self.params['n_patterns']):
            n_active_OR = np.random.randint(n_min_active_OR, n_max_active_OR)
            activated_ORs = random.sample(range(self.params['n_or']), n_active_OR)
            while np.unique(activated_ORs).size != n_active_OR:
                activated_ORs = random.sample(range(self.params['n_or']), n_active_OR)
            for OR in activated_ORs:
                dist = self.odorant_odor_distance_distribution((p1, p2, p3))
                distances[pn, OR] = dist
                affinity = np.exp(-(dist)**2 / self.params['distance_affinity_transformation_parameter_exp'])
                if affinity > 1.:
                    affinity = 1.
                if affinity < 0.:
                    affinity = 0.
                activation_matrix[pn, OR] = affinity
                n_above_thresh[pn] = (activation_matrix[pn ,:] > thresh).nonzero()[0].size

        print 'ActivationMatrix: per pattern number of activated glom: mean %.2f +- %.2f, \t%.2f +- %.2f percent' % (n_above_thresh.mean(), n_above_thresh.std(), \
                (n_above_thresh.mean() / float(self.params['n_or'])) * 100., (n_above_thresh.std() / float(self.params['n_or'])) * 100.)

        if self.params['OR_activation_normalization']:
            for pn in xrange(self.params['n_patterns']):
                # this normalization increases the likelihood of having ORs with a high affinity to 
                # each given pattern --> receptors have specialized to odorants (?)
#                if activation_matrix[pn, :].sum() == 0:
#                    print '\n\tWARNING: activation_matrix.sum for pattern %d == 0' % (pn)
                activation_matrix[pn, :] /= activation_matrix[pn, :].max() 
                if not activation_matrix[pn, :].sum() == 0:
                    activation_matrix[pn, :] /= activation_matrix[pn, :].sum() 
#                pass

#            for OR in xrange(self.params['n_or']):
#                activation_matrix[:, OR] /= activation_matrix[:, OR].sum()



        print 'Activation matrix sum:', activation_matrix.sum()
        print "Activation matrix fn:", self.params['activation_matrix_fn']
        np.savetxt(self.params['activation_matrix_fn'], activation_matrix)
        return activation_matrix



    def add_pattern_noise(self, activation_matrix):
        """
        This function loads an existing OR - odorant - affinity matrix 
        and adds noise to each pattern:
            each value in the matrix gets +- rand(0, degree_of_noise)
        """
        dgn = self.params['OR_affinity_noise']
        np.random.seed(self.params['OR_pattern_noise_seed'])
        for pn in xrange(self.params['n_patterns']):
            for OR in xrange(self.params['n_or']):
                activation_matrix[pn, OR] += np.random.uniform(-dgn, dgn)
                if (activation_matrix[pn, OR] > 1):
                    activation_matrix[pn, OR] = 1
                elif (activation_matrix[pn, OR] < 0):
                    activation_matrix[pn, OR] = 0.

        if self.params['OR_activation_normalization']:
#            for pn in xrange(self.params['n_patterns']):
#                activation_matrix[pn, :] /= activation_matrix[pn, :].max()
#                activation_matrix[pn, :] /= activation_matrix[pn, :].sum()
            for OR in xrange(self.params['n_or']):
                activation_matrix[:, OR] /= activation_matrix[:, OR].sum()

        np.savetxt(self.params['activation_matrix_fn_with_noise'], activation_matrix)
        return activation_matrix


    def modify_patterns_for_conc_inv(self, activation_matrix):
        """
        Change activation_matrix by to model different concentrations.
        activation_matrix values are modified additively
        """

        new_activation_matrix = np.zeros((self.params['n_patterns'], self.params['n_or']))
#        np.random.seed(self.params['seed_activation_matrix'])
#        selected_patterns = []
#        while len(selected_patterns) != self.params['n_patterns_test_conc_inv']:
#            selected_patterns = np.unique(np.random.randint(0, activation_matrix.shape[0], self.params['n_patterns_test_conc_inv']))
        selected_patterns = range(0, self.params['n_patterns_test_conc_inv'])

        conc_modifiers = np.linspace(-self.params['conc_inv_modifier'], self.params['conc_inv_modifier'], self.params['n_conc_check'])
        print 'Using the following patterns to check for concentration invariance:', selected_patterns
        print 'Using the following modifiers:', conc_modifiers

        # select n old patterns to modify
        pattern_index = 0
        for i_pn, pn in enumerate(selected_patterns):
            for i_conc in xrange(self.params['n_conc_check']):
                for OR in xrange(self.params['n_or']): 
                    activation = activation_matrix[pn, OR]# * self.params['scale_factor_activation_matrix']
                    if activation != 0:
                        new_activation = activation + conc_modifiers[i_conc]
                        if new_activation < 0:
                            new_activation = 0
                        elif new_activation > 1.:
                            new_activation = 1
                        new_activation_matrix[pattern_index, OR] = new_activation
                pattern_index += 1
        np.savetxt(self.params['activation_matrix_fn_conc_inv'], new_activation_matrix)
        print 'Saving modified activation matrix to:', self.params['activation_matrix_fn_conc_inv']
        return new_activation_matrix


    def set_odorant_distribution_params(self):
        """
        Distances are drawn from a tri-modal normal distributions
        This is done by choosing one of the three possible distributions 
        according to its probability p_i. 
        The probability with which one of the normal distributions is chosen 
        is determined by the respective integral (with respect to the other normal distributions).
        """
        p = self.params['odorant_receptor_distance_distribution_parameters']
        dist_range = self.params['odorant_receptor_distance_range']
        x = np.linspace(dist_range[0], dist_range[1], 1000)
        A1 = np.array(p[0] * mlab.normpdf(x, p[1], p[2])).sum() # integral of Gauss 1
        A2 = np.array(p[3] * mlab.normpdf(x, p[4], p[5])).sum() # integral of Gauss 2
        A3 = np.array(p[6] * mlab.normpdf(x, p[7], p[8])).sum() # integral of Gauss 3
        p1 = A1 / (A1 + A2 + A3)
        p2 = A2 / (A1 + A2 + A3)
        p3 = A3 / (A1 + A2 + A3)
        return p1, p2, p3


    def odorant_odor_distance_distribution(self, which_gauss):
        """
        Returns a distance between a virtual OR and a virtual odorant.
        This can in principle be any distribution.
        Here, we use a tri-modal Gaussian distribution which has been fitted to
        the real distance distribution gained by k-means clustering (the centroids being the ORs).
        For details see script: average_OR_affinity_distributions.py 
        Keyword arguments:
        which_gauss -- tuple of size three containing the probabilities with which each gauss distributions is picked
                        (p1, p2, p3) = which_gauss
        """
        (p1, p2, p3) = which_gauss
        p = self.params['odorant_receptor_distance_distribution_parameters']
        dist_range = self.params['odorant_receptor_distance_range']
        which_gauss = np.random.uniform(0, 1.)
        if which_gauss < p1:
            return np.random.normal(p[1], p[2])
        elif (which_gauss < p2 + p1):
            return np.random.normal(p[4], p[5])
        elif (which_gauss < p3 + p2 + p1):
            return np.random.normal(p[7], p[8])


    def create_parameters_for_concentration_sweep(self):
        """
        orthogonal patterns means a single odorant is presented to the system
        """
        self.create_activation_matrix_for_concentration_sweep()

        if self.params['OR_affinity_noise'] != 0.0:
            self.add_pattern_noise()

        for pn in xrange(self.params["n_patterns"]):
            # calculate the number of receptors activated by the current odorant or pattern
            for OR in xrange(self.params['n_or']): # equal to a row in the orn array
                # all ORNs within one row share the same Kd and concentration input value (because same OR is expressed by the these ORNs)
                activation = self.activation_matrix[pn, OR]# * self.params['scale_factor_activation_matrix']
#                if (activation > 1.0):
#                    activation = 1.0
                for x in xrange(self.params['n_orn_x']):
                    orn_id = OR * self.params['n_orn_x'] + x
                    self.voor[orn_id] = activation #* conc

            self.set_gor_values()
            self.set_gna_values()
            self.set_gk_values()
            self.set_gkcag_values()
            self.set_gcal_values()
            self.set_gleak_values()
            self.set_tau_cadec_values()
            output_fn = self.orn_params_fn_base + "%d.dat" % (pn)
#            print "Writing orn parameters for pattern \t... %d / %d to file: %s" % (pn, self.n_patterns, output_fn)
            self.write_params_to_file(output_fn)
        return 1
        # stop new






    def create_params_for_response_curve(self):
        """
        This function creates the parameters needed to
        measure a response curve or concentration sweep
        (one odorant with maximum affinity and different concentrations
        """
        output_fn = self.orn_params_fn_base + "0.dat"
        self.set_oor_value_for_conc_sweep()
        self.set_gor_values()
        self.set_gna_values()
        self.set_gk_values()
        self.set_gkcag_values()
        self.set_gcal_values()
        self.set_gleak_values()
        self.set_tau_cadec_values()
        self.write_params_to_file(output_fn)
        return 1


    def set_oor_value_for_conc_sweep(self):
        # c / Kd values are equally distributed on a log10 scale
        # c_Kd = 10**exp_min .. 10**exp_max
        exp_min = -2.0
#        exp_min = -1.5
        exp_max = 2.0
        z = (exp_max - exp_min) / float(self.n_or)
        c_Kd = []
        c_Kd.append(10**exp_min)
        for i in xrange(self.n_orn_y - 1):
            c_Kd.append(c_Kd[-1] * 10**z)

        # all ORNs within one row share the same Kd and concentration input value (because same OR is expressed by the these ORNs) and concentration increases with the row
        for y in xrange(self.n_orn_y):
            for x in xrange(self.n_orn_x): # increasing gor
                orn_id = y * self.n_orn_x + x
                self.vc_Kd[orn_id] = c_Kd[y]
                self.voor[orn_id] = 1 - 1. / (1 + c_Kd[y])


    def set_gor_values(self):

        def gor_func(x, a1, a2, exp):
            """this function calculates the gor for  a given x """
            return a1 * x**exp + a2
        # --------- CREATE ORN PARAMETERS -------------
        # gor values are distributed between a min and max value
        # gor(x) = a * x^gor_exp + b
        # gor(x_min) = b = gor_min
        # gor(x_max) = gor_max
        # parameters for gor_func, i.e. how gor is distributed among the n_gor different ORNs expressing the same OR
        x_min = 0
        x_max = self.n_orn_x

        gor_exp = self.params["gor_exp"]
        b = self.params["gor_min"]
        a = (self.params["gor_max"] - self.params["gor_min"]) / ((x_max - 1)**self.params["gor_exp"] - x_min**self.params["gor_exp"])
        for i in xrange(self.n_orn_y):  # y-axis
            for j in xrange(self.n_orn_x):
                orn_id = i * self.n_orn_x + j
                gor_value = gor_func(j, a, b, gor_exp)
                self.vgor[orn_id] = gor_value
                self.gor_values[j] = gor_value
        return self.gor_values


    def set_gna_values(self):
        for i in xrange(self.n_orn_y):  # y-axis
            for j in xrange(self.n_orn_x):
                orn_id = i * self.n_orn_x + j
                # map j into the interval [gor_min, gor_max] which is the x value for the function gna(gor)
                self.vgna[orn_id] = self.params["gna"]


    def set_gk_values(self):
        for i in xrange(self.n_orn_y):  # y-axis
            for j in xrange(self.n_orn_x):
                orn_id = i * self.n_orn_x + j
                # map j into the interval [gor_min, gor_max] which is the x value for the function gk(gor)
                self.vgk[orn_id] = self.params["gk"]


    def set_gkcag_values(self):
        # before this function is called set_gor_values must be called
        gkcag_values = self.linear_transformation(self.gor_values, self.params['gkcag_params'][0], self.params['gkcag_params'][1])

        for i in xrange(self.n_orn_y):  # y-axis
            for j in xrange(self.n_orn_x):
                orn_id = i * self.n_orn_x + j
                self.vgkcag[orn_id] = gkcag_values[j]

                # map j into the interval [gor_min, gor_max] which is the x value for the function gkcag(gor)
#                x = self.vgor[orn_id]
#                gkcag = p[0] * np.log10(x) + p[1]
#                self.vgkcag[orn_id] = gkcag


    def set_gcal_values(self):
        # before this function is called set_gor_values must be called
        gcal_values = self.linear_transformation(self.gor_values, self.params['gcal_params'][0], self.params['gcal_params'][1])
        for i in xrange(self.n_orn_y):  # y-axis
            for j in xrange(self.n_orn_x):
                orn_id = i * self.n_orn_x + j
                self.vgcal[orn_id] = gcal_values[j]
                # map j into the interval [gor_min, gor_max] which is the x value for the function gcal(gor)
#                x = self.vgor[orn_id]
#                gcal = p[0] * np.log10(x) + p[1]
#                self.vgcal[orn_id] = gcal


    def set_gleak_values(self):
        # before this function is called set_gor_values must be called
        gleak_values = self.linear_transformation(self.gor_values, self.params['gleak_params'][0], self.params['gleak_params'][1])


        for i in xrange(self.n_orn_y):  # y-axis
            for j in xrange(self.n_orn_x):
                orn_id = i * self.n_orn_x + j
                self.vgleak[orn_id] = gleak_values[j]

                # map j into the interval [gor_min, gor_max] which is the x value for the function gleak(gor)
#                x = self.vgor[orn_id]
#                gleak = p[0] * np.exp(p[1] * x ** 2) + p[2]
#                self.vgleak[orn_id] = gleak


    def set_tau_cadec_values(self):
        # before this function is called set_gor_values must be called
        for i in xrange(self.n_orn_y):  # y-axis
            for j in xrange(self.n_orn_x):
                orn_id = i * self.n_orn_x + j
                # map j into the interval [gor_min, gor_max] which is the x value for the function tau_cadec(gor)
                self.vtau_cadec[orn_id] = self.params["tau_cadec"]


    def write_params_to_file(self, output_fn):
        """
        Write the ORN parameters into a NEURON readable file
        """
#        print "writing orn params to ", output_fn
        orn_pf = file(output_fn, 'w')
        num_params = 9 + 1 # + 1 for the id
        first_line = "%d %d\n"  % (self.n_orn, num_params) # neuron readability
        orn_pf.write(first_line)

        # write these values to file
        for row in xrange(self.params["n_orn_y"]):
            for col  in xrange(self.params["n_orn_x"]):
                orn_id = col + row * self.params["n_orn_x"]
                line = "%d\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n" % \
                        (orn_id, self.voor[orn_id], self.vc_Kd[orn_id], self.vgor[orn_id], \
                            self.vgna[orn_id], self.vgk[orn_id], self.vgkcag[orn_id], self.vgcal[orn_id],\
                            self.vgleak[orn_id], self.vtau_cadec[orn_id])
                orn_pf.write(line)
        orn_pf.close()


    def write_current_param_values_to_file(self, output_fn=None):
        """
        This function is required for hand tuning the orn parameters.
        """
        if output_fn == None:
            output_fn = self.params['orn_params_fn_base'] + '0.dat'
        print "writing orn params to ", output_fn
        orn_pf = file(output_fn, 'w')
        num_params = 9 + 1 # + 1 for the id
        first_line = "%d %d\n"  % (self.n_orn, num_params) # neuron readability
        orn_pf.write(first_line)

        # write these values to file
        for row in xrange(self.params["n_orn_y"]):
            for col  in xrange(self.params["n_orn_x"]):
                orn_id = col + row * self.params["n_orn_x"]
                line = "%d\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n" % (orn_id, \
                        self.voor[orn_id], self.vc_Kd[orn_id], self.vgor[orn_id], \
                            self.current_param_values[col][0],\
                            self.current_param_values[col][1],\
                            self.current_param_values[col][2],\
                            self.current_param_values[col][3],\
                            self.current_param_values[col][4],\
#                            self.vgleak[orn_id],
                            self.current_param_values[col][5])
                orn_pf.write(line)
        orn_pf.flush()
        orn_pf.close()


    def linear_transformation(self, x, y_min, y_max):
        """
        x : the range to be transformed
        y_min, y_max : lower and upper boundaries for the range into which x
        is transformed to
        Returns y = f(x), f(x) = m * x + b
        """
        x_min = np.min(x)
        x_max = np.max(x)
        if x_min == x_max:
            x_max = x_min * 1.0001
        return (y_min + (y_max - y_min) / (x_max - x_min) * (x - x_min))


