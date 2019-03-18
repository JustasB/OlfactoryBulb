import re
import os
import sys
import numpy as np
import pylab
import simulation_parameters
from FigureCreator import plot_params
import json


class AnalyseResponseCurves(object):
    """
    Calculates for a set of response curves the use defined fitness
    """

    def __init__(self, params):
        self.params = params


    def sweep_orn_response_fitness(self):
        match_files = 'orn_response_curve_(\d+).dat'
#        match_files = 'orn_response_curve_%d.dat'
        file_dict = self.find_files(match_files)
        n_matches = len(file_dict.keys())
        assert (n_matches > 0), 'No files found matching the pattern %s in folder %s\nCorrect folder in simulation_parameters.py?' % (match_files, self.params['other_folder'])
        fitness_values = np.zeros(n_matches)
        print 'n_matches', n_matches
        for i in xrange(n_matches):
            fn = file_dict[i]
            d = np.loadtxt(fn)
            fnv = self.get_orn_response_fitness(d)
            print 'Analysing', fn, fnv
            fitness_values[i] = fnv

        n_best = 50
        self.best_sims = fitness_values.argsort()[-n_best:]
        print 'Best simulations and fitness values'
        print 'sim_ids:', self.best_sims
        print 'fitness_values:', fitness_values[self.best_sims]
        output_fn = 'best_sim_ids.log'
        f = file(output_fn, 'w')
        print 'Printing best sim IDs to:', output_fn
        json.dump(list(self.best_sims), f, indent=0)


    def show_best_sim_results(self):
        fig_fns = []
        display_command = 'ristretto '
        for i_ in xrange(len(self.best_sims)):
            fig_fn = self.params['figure_folder'] + '/hand_tuned_%d.png' % self.best_sims[i_]
            fig_fns.append(fig_fn)
            display_command += '%s ' % fig_fn
            print fig_fn, 
        print ' '
        output_fn = self.params['other_folder'] + '/good_figure_filenames.json' 
        output_f = file(output_fn, 'w')

#        print 'debug', type(fig_fns)
        print 'Writing good figure filenames to:', output_fn
        os.system(display_command)
        json.dump(fig_fns, output_f, indent=0)



    def get_orn_response_fitness(self, xy_data):
        """
        Returns the fitness value for a set of response curves.

        Keyword arguments:
        xy_data -- numpy.ndarray
        xy_data[:, 0] = concentration values
        xy_data[:, i] = output-rate of curve i
        """
        fitness = 0.
        
        x_data = xy_data[:, 0]
        y_data = xy_data[:, 1:]
        n_x = x_data.size
        n_curves = xy_data[0, 1:].size # the number of curves to be analysed

        # define some fixed characteristica to be applied to the set of curves
        # for each curve there are some individual criteria, i.e. how the curve should look like
        y_min, y_max = .05, .9 # percentage of the curve's global y-min (max) which marks the x_min, x_max points
        x_mins, x_maxs = (-1) * np.ones(n_curves), (-1) * np.ones(n_curves) # negative to check if they've already been set
        y_max_global = np.zeros(n_curves)

        punishment = 1
        # get the curve characteristics
        for curve_ in xrange(n_curves):
            y_max_curve = y_max * y_data[:, curve_].max()
            y_min_curve = y_min * y_max_curve
            y_max_global[curve_] = y_data[:, curve_].max()
            for i_ in xrange(n_x):
                x, y = x_data[i_], y_data[i_, curve_]
                if y > y_min_curve and (x_mins[curve_] == -1.):
                    x_mins[curve_] = x
                if y > y_max_curve and (x_maxs[curve_] == -1.):
                    x_maxs[curve_] = x

        # compare the curves within the set
        for curve_ in xrange(1, n_curves):
            # compare global maxima
            dy = y_max_global[curve_] - y_max_global[curve_-1]
            if dy < 0: # quite bad
#                fitness += dy
                fitness -= punishment
            else:
                fitness += punishment
            # compare the point when curves cross the y_min threshold
            dx_mins = x_mins[curve_] - x_mins[curve_-1]
            if dx_mins < 0:
                fitness -= punishment
            else:
                fitness += punishment
            # compare the point when curves cross the y_max threshold
            dx_maxs = x_maxs[curve_] - x_maxs[curve_-1]
            if dx_maxs < 0:
                fitness -= punishment
            else:
                fitness += punishment

            # make sure that the curves don't cross
            for i_ in xrange(n_x):
                y = y_data[i_, curve_]
                y_prv = y_data[i_, curve_-1]
                if y < y_prv:
                    fitness -= 100 * punishment
                else:
                    fitness += 10 * punishment

        # ensure monotony
        for curve_ in xrange(n_curves):
            for i_ in xrange(1, n_x):
                if y_data[i_, curve_] < y_data[i_-1, curve_]:
                    fitness -= 10 * punishment
                else:
                    fitness += punishment

        # compare the highest and smallest ymax values of all curves
        dy_max_global = np.max(y_max_global) - np.min(y_max_global)
        fitness -= dy_max_global / (np.min(y_max_global) + 1e-6)
        # punish large differences in maximum output rates
        if dy_max_global > 60:
            fitness -= 10 * punishment
        else:
            fitness += 10 * punishment

        # punish too small output rates
        if np.min(y_max_global) < 5:
            fitness -= 10 * punishment

        # punish too high max output rates
        if np.max(y_max_global) > 160:
            fitness -= 10 * punishment

        # compare the points when curves cross their y_min
        dx_min_global = x_mins[0] - x_mins[-1]
        fitness += dx_min_global * 10.

        if np.max(y_max_global) < 10:
            fitness -= 10e6 * punishment
        else:
            fitness += 10 * punishment

        return fitness




    def find_files(self, pattern, folder=None):
        """
        Returns a dictionary with all absolute file names 
        which match a certain pattern in a given folder.
        The key of the file dictionary is the simulation ID 
        and the value is the absolute path.
        """

        if folder == None:
            folder = self.params['other_folder']
        file_dict = {}
        for f in os.listdir(folder):
            m = re.match(pattern, f)
            if m:
                sim_id = int(m.groups()[0])
                file_dict[sim_id] = os.path.abspath(folder + '/' + f)
        return file_dict


    def create_new_parameters_from_selected_sets(self, list_of_sim_ids=None, param_log_fn_base=None, n_new_sets=None):
        """
        list_of_sim_ids: list of integer values with the simulation ids that contain good parameter sets
        param_log_fn_base: the file which contains the sim_id and the parameters to all the simulations.
        """

        if param_log_fn_base == None:
            param_log_fn_base = 'orn_sweep_params.log'

        old_params = np.loadtxt(param_log_fn_base)
        n_params = old_params.shape[1] - 1
        if list_of_sim_ids == None:
            list_of_sim_ids = self.best_sims
        assert (len(list_of_sim_ids) > 1), 'The length of suggested simulations should be longer than one'
        if n_new_sets == None:
            n_new_sets = len(list_of_sim_ids)
        new_params = np.zeros((n_new_sets, old_params.shape[1]))

#        np.random.seed(0)
        rnd_mod = .05
        for i_ in xrange(n_new_sets):
            old_set_id = np.random.randint(0, len(self.best_sims))
            original_set = old_params[old_set_id, :]
            rnd_factors = rnd_mod * 2. * np.random.random(n_params) + (1 - rnd_mod)
            new_set = rnd_factors * original_set[1:]
            new_params[i_, 1:] = new_set



        new_params[:, 0] = np.arange(n_new_sets)
        output_fn = 'new_orn_params.dat'
        print 'Saving new parameters to', output_fn
        np.savetxt(output_fn, new_params)


    def combine_good_parameter_sets_to_new(self, n_new_sets=100, list_of_sim_ids=None, param_log_fn_base=None, output_fn=None, sim_id_offset=0):
        if param_log_fn_base == None:
            param_log_fn_base = 'orn_sweep_params.log'

        old_params = np.loadtxt(param_log_fn_base)
        n_params = old_params.shape[1] - 1
        if list_of_sim_ids == None:
            list_of_sim_ids = self.best_sims
        assert (len(list_of_sim_ids) > 1), 'The length of suggested simulations should be longer than one'
        if n_new_sets == None:
            n_new_sets = len(list_of_sim_ids)
        new_params = np.zeros((n_new_sets, old_params.shape[1]))

        rnd_mod = .1
        for i_ in xrange(n_new_sets):
            old_set_id_0 = np.random.randint(0, len(self.best_sims))
            original_set_0 = old_params[old_set_id_0, 1:]
            old_set_id_1 = np.random.randint(0, len(self.best_sims))
            original_set_1 = old_params[old_set_id_1, 1:]
            new_set = .5 * (original_set_0 + original_set_1)
            new_params[i_, 1:] = new_set

        new_params[:, 0] = np.arange(n_new_sets) + sim_id_offset
        if output_fn == None:
            output_fn = 'new_combined_orn_params.dat'
        print 'Saving new parameters to', output_fn
        np.savetxt(output_fn, new_params)


if __name__ == '__main__':

    param_tool = simulation_parameters.parameter_storage()
    params = param_tool.params

    ARC = AnalyseResponseCurves(params) 
    ARC.sweep_orn_response_fitness()
    n_random_sets = 1000
    ARC.create_new_parameters_from_selected_sets(n_new_sets=n_random_sets, param_log_fn_base='orn_sweep_params_2.log')
    ARC.combine_good_parameter_sets_to_new(n_new_sets=1000, param_log_fn_base='orn_sweep_params_2.log', sim_id_offset=n_random_sets)
    ARC.show_best_sim_results()
