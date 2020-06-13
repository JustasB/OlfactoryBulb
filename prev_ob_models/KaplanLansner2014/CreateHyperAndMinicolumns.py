import os
import sys
import numpy
import numpy.random as rnd
from scipy.spatial.distance import euclidean


class CreateHyperAndMinicolumns(object):

    def __init__(self, param_dict):
        self.params = param_dict
        self.folder_name = self.params['folder_name']
        self.n_hc = self.params['n_hc']
        self.n_mc = self.params['n_mc']
        self.n_pyr = self.n_hc * self.n_mc * self.params['n_pyr_per_mc']
        self.n_pyr_per_mc = self.params['n_pyr_per_mc']
        self.n_rsnp_per_mc = self.params['n_rsnp_per_mc']
        self.n_rsnp = self.params['n_rsnp']
        self.n_basket_per_hc = self.params['n_basket_per_hc']
        self.n_basket = self.params['n_basket']
        self.n_rsnp = self.params['n_rsnp']
        self.p_rsnp_pyr = self.params["p_rsnp_pyr"]
        self.w_rsnp_pyr = self.params["w_rsnp_pyr"]
        self.p_pyr_pyr_local = self.params["p_pyr_pyr_local"]
        self.w_pyr_pyr_local = self.params["w_pyr_pyr_local"]
        self.p_pyr_pyr_global = self.params["p_pyr_pyr_global"]
        self.p_pyr_basket = self.params["p_pyr_basket"]
        self.w_pyr_basket = self.params["w_pyr_basket"]
        self.p_basket_pyr = self.params["p_basket_pyr"]
        self.w_basket_pyr = self.params["w_basket_pyr"]
        self.p_pyr_rsnp = self.params["p_pyr_rsnp"]
        self.p_basket_basket = self.params["p_basket_basket"]
        self.w_basket_basket = self.params["w_basket_basket"]
        rnd.seed(self.params['seed_connections'])
        self.output_fn = self.params['conn_list_layer23']

    def get_closest_basket_cell_ids(self, mc):
        """
        mc = minicolumn index of the pyramidal cell
        returns the a list of gids of the n_tgt_basket_per_mc closest basket cells to which the pyr cells in the given minicolum shall connect to
        list of gids is without the offset
        """
        n_mc = self.params['n_mc'] * self.params['n_hc']
        pyr_grid_x = int(numpy.sqrt(self.params['n_mc']))
        if (numpy.sqrt(self.params['n_mc'] / int(numpy.sqrt(self.params['n_mc'])) != 1.0)):
            # if self.params['n_mc'] can not be put on a rectangular grid
            pyr_grid_y = int(numpy.sqrt(self.params['n_mc'])) + 1
        else:
            pyr_grid_y = int(numpy.sqrt(self.params['n_mc']))

        # x, y positions of the minicolumn in the HC grid
        pyr_pos_x = mc % pyr_grid_x
        pyr_pos_y = mc / pyr_grid_x

        basket_grid_x = int(numpy.sqrt(self.params['n_basket_per_hc']))
        if ((numpy.sqrt(self.params['n_basket_per_hc']) / int(numpy.sqrt(self.params['n_basket_per_hc']))) != 1.0):
            # if self.params['n_basket_per_hc'] can not be put on a rectangular grid
            basket_grid_y = int(round(numpy.sqrt(self.params['n_basket_per_hc'])) + 1)
        else:
            basket_grid_y = int(numpy.sqrt(self.params['n_basket_per_hc']))

        distances_to_basket_cells = numpy.zeros(self.params['n_basket_per_hc'])
        for basket_id in xrange(self.params['n_basket_per_hc']):
            x = basket_id % basket_grid_x
            y = basket_id / basket_grid_x
#        for y in xrange(basket_grid_y):
#            for x in xrange(basket_grid_x):
            scaled_basket_pos_x = (float(x) / basket_grid_x) * pyr_grid_x
            scaled_basket_pos_y = (float(y) / basket_grid_y) * pyr_grid_y
#            basket_id = y * basket_grid_x + x

#            print "basket_id", basket_id, basket_grid_y, basket_grid_x, (numpy.sqrt(self.params['n_basket_per_hc']) / int(numpy.sqrt(self.params['n_basket_per_hc'])) != 1.0)
#            print "if:", (numpy.sqrt(self.params['n_basket_per_hc']) / int(numpy.sqrt(self.params['n_basket_per_hc']))) != 1.0
#            print "check:", numpy.sqrt(self.params['n_basket_per_hc']) / int(numpy.sqrt(self.params['n_basket_per_hc']))

#            print "basket_grid_y", numpy.sqrt(self.params['n_basket_per_hc'] / int(numpy.sqrt(self.params['n_basket_per_hc']))), (numpy.sqrt(self.params['n_basket_per_hc'] / int(numpy.sqrt(self.params['n_basket_per_hc']))) != 1.0)
#            print "debug basket_id %d   x %f    y %f  distance to mc %d: %f" % (basket_id, scaled_basket_pos_x, scaled_basket_pos_y, mc, distances_to_basket_cells[basket_id])
#            print "debug pyr_pos x %d    y %d" % (pyr_pos_x, pyr_pos_y)
#            print "d", distances_to_basket_cells.size, self.params['n_basket_per_hc']
            distances_to_basket_cells[basket_id] = euclidean([pyr_pos_x, pyr_pos_y], [scaled_basket_pos_x, scaled_basket_pos_y])

        sorted_ids = distances_to_basket_cells.argsort()
        closest_basket_cells = numpy.take(sorted_ids, range(self.params['n_tgt_basket_per_mc']))
        closest_basket_cells.sort()
        return closest_basket_cells


    def create_connections(self):
        my_min_pyr = 0
        my_min_rsnp = 0
        my_min_basket = 0
        my_max_pyr = self.n_hc * self.n_mc * self.n_pyr_per_mc
        my_max_rsnp = self.n_hc * self.n_mc * self.n_rsnp_per_mc
        my_max_basket = self.n_hc * self.n_basket_per_hc

        lines = ""
        n_connections = 0 # needed for NEURON output file

        # within one MC:
        # 1) pyr - pyr: 25%, 1.2
        # 2) rsnp - pyr: 70%, -0.8
        # within one HC:
        # 3) pyr - basket: 70%,, 1.9
        # 4) basket - basket: 70%,  -2.2

        # from one HC to another:
        # 5) pyr - rsnp: 30%, 0.2 if in different patterns
        # 6) pyr - pyr: 30%, 0.4 if in the same patterns

        # iterate over source pyramidal cells
        for pyr in xrange(my_min_pyr, my_max_pyr):
            # 1) - 4) within same HC
            # 1 PYR -> PYR: same HC, same MC
            # coordinates of source pyr
            mc_pyr = (pyr / self.n_pyr_per_mc) % self.n_mc # mc_pyr = 0 .. n_mc
            hc_pyr = pyr / (self.n_pyr_per_mc * self.n_mc)
            # min/max ids of possible target pyr cells
            min_pyr = mc_pyr * self.n_pyr_per_mc + hc_pyr * self.n_mc * self.n_pyr_per_mc
            max_pyr = (mc_pyr + 1) * self.n_pyr_per_mc + hc_pyr * self.n_mc * self.n_pyr_per_mc
            
            for tgt_pyr in xrange(min_pyr, max_pyr):
                if (pyr != tgt_pyr):
                    w = self.draw_connection(self.p_pyr_pyr_local, self.w_pyr_pyr_local, noise=self.params["w_pyr_pyr_local_sigma"])
                    if ((w != 0) and (w > self.params['weight_threshold'])):
                        lines += "%d\t%d\t%.6e\n" % (pyr + self.params['pyr_offset'], tgt_pyr + self.params['pyr_offset'], w)
                        n_connections += 1

            # 2 RSNP -> PYR: within same minicolumn
            min_rsnp = mc_pyr * self.n_rsnp_per_mc
            max_rsnp = (mc_pyr + 1) * self.n_rsnp_per_mc
#            for rsnp in xrange(min_rsnp, max_rsnp): # global loop for inter hypercolumnar connections
            for rsnp in xrange(0, self.n_rsnp_per_mc): # global loop for inter hypercolumnar connections
#                mc_rsnp = rsnp / self.n_rsnp_per_mc
#                hc_rsnp = rsnp / (self.n_rsnp_per_mc * self.n_mc)
                w = self.draw_connection(self.p_rsnp_pyr, self.w_rsnp_pyr, noise=self.params["w_rsnp_pyr_sigma"])
                if ((w != 0) and (w > self.params['weight_threshold'])):
                    src_id = rsnp + self.params['rsnp_offset'] + mc_pyr * self.n_rsnp_per_mc + hc_pyr * self.n_mc * self.n_rsnp_per_mc
                    lines += "%d\t%d\t%.6e\n" % (src_id, pyr + self.params['pyr_offset'] , (-1.0) * w)
                    n_connections += 1


            # pyr within one minicolumn connect to the 8 'closest' basket cells
            min_basket = self.n_hc * self.params['n_basket_per_hc'] # within this hypercolumn
            max_basket = (self.n_hc + 1) * self.params['n_basket_per_hc']
#            for basket in xrange(min_basket, max_basket):
            # 3) BASKET -> PYR
            basket_gid_offset = self.params['basket_offset'] + hc_pyr * self.params['n_basket_per_hc']
            for basket in xrange(0, self.params['n_basket_per_hc']):
                w = self.draw_connection(self.p_basket_pyr, self.w_basket_pyr, noise=self.params["w_basket_pyr_sigma"])
                if ((w != 0) and (w > self.params['weight_threshold'])):
                    lines += "%d\t%d\t%.6e\n" % (basket + basket_gid_offset , pyr+ self.params['pyr_offset'] , (-1.0) * w)
                    n_connections += 1

            # 4) PYR -> BASKET
            basket_gids_without_offset = self.get_closest_basket_cell_ids(mc_pyr)
            for basket in basket_gids_without_offset:
                w = self.draw_connection(self.p_pyr_basket, self.w_pyr_basket, noise=self.params["w_pyr_basket_sigma"])
                if ((w != 0) and (w > self.params['weight_threshold'])):
                    lines += "%d\t%d\t%.6e\n" % (pyr + self.params['pyr_offset'] , basket + basket_gid_offset, w)
                    n_connections += 1

        # BASKET -> BASKET
        for hc in xrange(self.n_hc):
            for src in xrange(0, self.params['n_basket_per_hc']):
                for tgt in xrange(0, self.params['n_basket_per_hc']):
                    if (src != tgt):
                        w = self.draw_connection(self.p_basket_basket, self.w_basket_basket, noise=self.params["w_basket_basket_sigma"])
                        if ((w != 0) and (w > self.params['weight_threshold'])):
                            gid_offset = self.params['basket_offset'] + hc * self.params['n_basket_per_hc']
                            lines += "%d\t%d\t%.6e\n" % (src + gid_offset, tgt + gid_offset, (-1.0) * w)
                            n_connections += 1

        first_line = "%d %d\n" % (n_connections, 3)
        print "Writing connections to file", self.output_fn
        self.output_f = file(self.output_fn, 'w')
        self.output_f.write(first_line)
        self.output_f.write(lines)
        self.output_f.close()


    def create_orthogonal_connections(self):
        my_min_pyr = 0
        my_min_rsnp = 0
        my_min_basket = 0
        my_max_pyr = self.n_hc * self.n_mc * self.n_pyr_per_mc
        my_max_rsnp = self.n_hc * self.n_mc * self.n_rsnp_per_mc
        my_max_basket = self.n_hc * self.n_basket_per_hc

        lines = ""
        n_connections = 0 # needed for NEURON output file

        # within one MC:
        # 1) pyr - pyr: 25%, 1.2
        # 2) rsnp - pyr: 70%, -0.8
        # within one HC:
        # 3) pyr - basket: 70%,, 1.9
        # 4) basket - basket: 70%,  -2.2

        # from one HC to another:
        # 5) pyr - rsnp: 30%, 0.2 if in different patterns
        # 6) pyr - pyr: 30%, 0.4 if in the same patterns

        # iterate over source pyramidal cells
        for pyr in xrange(my_min_pyr, my_max_pyr):
            # 1) - 4) within same HC
            # 1 PYR -> PYR: same HC, same MC
            # coordinates of source pyr
            mc_pyr = (pyr / self.n_pyr_per_mc) % self.n_mc
            hc_pyr = pyr / (self.n_pyr_per_mc * self.n_mc)
            assert (hc_pyr< self.params['n_hc'])
            # min/max ids of possible target pyr cells
            min_pyr = mc_pyr * self.n_pyr_per_mc + hc_pyr * self.n_mc * self.n_pyr_per_mc
            max_pyr = (mc_pyr + 1) * self.n_pyr_per_mc + hc_pyr * self.n_mc * self.n_pyr_per_mc
            
            for tgt_pyr in xrange(min_pyr, max_pyr):
                tgt_hc = tgt_pyr / (self.n_pyr_per_mc * self.n_mc)
                if (pyr != tgt_pyr) and (tgt_hc == hc_pyr):
                    w = self.draw_connection(self.p_pyr_pyr_local, self.w_pyr_pyr_local, noise=self.params["w_pyr_pyr_local_sigma"])
                    if (w != 0):
                        lines += "%d\t%d\t%.6e\n" % (pyr + self.params['pyr_offset'], tgt_pyr + self.params['pyr_offset'], w)
                        n_connections += 1

            # 2 RSNP -> PYR: within same minicolumn
            min_rsnp = mc_pyr * self.n_rsnp_per_mc
            max_rsnp = (mc_pyr + 1) * self.n_rsnp_per_mc
#            for rsnp in xrange(min_rsnp, max_rsnp): # global loop for inter hypercolumnar connections
            for rsnp in xrange(0, self.n_rsnp_per_mc): # global loop for inter hypercolumnar connections
#                mc_rsnp = rsnp / self.n_rsnp_per_mc
#                hc_rsnp = rsnp / (self.n_rsnp_per_mc * self.n_mc)
                w = self.draw_connection(self.p_rsnp_pyr, self.w_rsnp_pyr, noise=self.params["w_rsnp_pyr_sigma"])
                if (w != 0):
                    gid_offset = self.params['rsnp_offset'] + hc_pyr * self.n_mc * self.n_rsnp_per_mc + mc_pyr * self.n_rsnp_per_mc
                    lines += "%d\t%d\t%.6e\n" % (rsnp + gid_offset, pyr + self.params['pyr_offset'] , (-1.0) * w)
                    n_connections += 1

            min_basket = self.n_hc * self.n_basket_per_hc
            max_basket = (self.n_hc + 1) * self.n_basket_per_hc
#            for basket in xrange(min_basket, max_basket):
            for basket in xrange(0, self.n_basket_per_hc):
                # 3) BASKET -> PYR
                w = self.draw_connection(self.params['p_basket_pyr'], self.w_basket_pyr, noise=self.params["w_basket_pyr_sigma"])
                gid_offset = self.params['basket_offset'] + hc_pyr * self.n_basket_per_hc
                if (w != 0):
                    lines += "%d\t%d\t%.6e\n" % (basket + gid_offset, pyr + self.params['pyr_offset'] , (-1.0) * w)
                    n_connections += 1

            # 4) PYR -> BASKET
            basket_gids_without_offset = self.get_closest_basket_cell_ids(mc_pyr)
            basket_gid_offset = self.params['basket_offset'] + hc_pyr * self.n_basket_per_hc
            for basket in basket_gids_without_offset:
                w = self.draw_connection(self.p_pyr_basket, self.w_pyr_basket, noise=self.params["w_pyr_basket_sigma"])
                if ((w != 0) and (w > self.params['weight_threshold'])):
                    lines += "%d\t%d\t%.6e\n" % (pyr + self.params['pyr_offset'] , basket + basket_gid_offset, w)
                    n_connections += 1


            # all connections targeting another hypercolumn
            for tgt_hc in xrange(self.n_hc):
                if (tgt_hc != hc_pyr): # only if target cell is in another HC
                    for tgt_mc in xrange(self.n_mc):
                        same_pattern = (mc_pyr == tgt_mc) # compare mc coordinate
                        if not same_pattern:
                        # 5) PYR - RSNP: 30%, 0.2 if in different patterns
                            # calculate min and max ids of target pyramidal cells
#                            min_rsnp = tgt_hc * self.n_mc * self.n_rsnp_per_mc + tgt_mc * self.n_rsnp_per_mc
#                            max_rsnp = tgt_hc * self.n_mc * self.n_rsnp_per_mc + (tgt_mc + 1) * self.n_rsnp_per_mc
#                            for tgt_rsnp in xrange(min_rsnp, max_rsnp):
                            for tgt_rsnp in xrange(0, self.n_rsnp_per_mc):
                                w = self.draw_connection(self.p_pyr_rsnp, self.params['w_pyr_rsnp_max'], noise=self.params["w_pyr_rsnp_sigma"])
                                if (w != 0):
                                    gid_offset = self.params['rsnp_offset'] + tgt_hc * self.n_mc * self.n_rsnp_per_mc + tgt_mc * self.n_rsnp_per_mc
                                    lines += "%d\t%d\t%.6e\n" % (pyr + self.params['pyr_offset'], tgt_rsnp + gid_offset, w)
                                    n_connections += 1


                        else: # same_pattern == True
                        # 6) PYR - PYR : 30%, 0.4 if in the same patterns
                            min_pyr = tgt_hc * self.n_mc * self.n_pyr_per_mc + tgt_mc * self.n_pyr_per_mc
                            max_pyr = tgt_hc * self.n_mc * self.n_pyr_per_mc + (tgt_mc + 1) * self.n_pyr_per_mc
                            for tgt_pyr in xrange(min_pyr, max_pyr):
                                w = self.draw_connection(self.p_pyr_pyr_global, self.params['w_pyr_pyr_global_max'], noise=self.params["w_pyr_pyr_global_sigma"])
                                if (w != 0):
                                    lines += "%d\t%d\t%.6e\n" % (pyr + self.params['pyr_offset'], tgt_pyr + self.params['pyr_offset'], w)
                                    n_connections += 1




        # BASKET -> BASKET
        for hc in xrange(self.n_hc):
            for src in xrange(0, self.n_basket_per_hc):
                for tgt in xrange(0, self.n_basket_per_hc):
                    if (src != tgt):
                        w = self.draw_connection(self.p_basket_basket, self.w_basket_basket, noise=self.params["w_basket_basket_sigma"])
                        if (w != 0):
                            gid_offset = self.params['basket_offset'] + hc * self.n_basket_per_hc
                            lines += "%d\t%d\t%.6e\n" % (src + gid_offset, tgt + gid_offset, (-1.0) * w)
                            n_connections += 1

        first_line = "%d %d\n" % (n_connections, 3)
#        first_line = "%d %d\n" % (3, n_connections)
        print "Writing connections to file", self.output_fn
        self.output_f = file(self.output_fn, 'w')
        self.output_f.write(first_line)
        self.output_f.write(lines)
        self.output_f.close()

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
                if (numpy.sign(weight) != numpy.sign(w)):
                    return 0
                return weight
            elif (noise < 0): # stupid user, noise should be > 0
                print "WARNING, negative noise given to functions.draw_connection(p, w, noise)!"
                noise *= (-1.0)
                weight = rnd.normal(w, noise)
                if (numpy.sign(weight) != numpy.sign(w)):
                    return 0
                return weight
            elif (noise == 0):
                return w
        return 0


