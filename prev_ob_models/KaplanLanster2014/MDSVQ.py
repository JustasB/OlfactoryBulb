import numpy as np
import os 
import sys
import time
# classes for setting up connectivity and the individual cell parameters
import MergeSpikefiles
import Orange
import scipy.cluster.vq as scvq
import pylab
#from scipy.cluster.vq import vq, whiten, kmeans, kmeans2
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
from scipy.spatial.distance import euclidean as euclidean_dist

class MDSVQ(object):
    """
    """

    def __init__(self, params):
        self.params = params


    def mds(self, input_fn_or_array, mds_output_fn, n_steps=5000, thresh=1e-4, cell_type='mit'):
        if (type(input_fn_or_array) == type('')):
            input_fn = input_fn_or_array
            print "MDSVQ.calculate_distances_from_normalized_output from :", input_fn
            try: 
                input_array = np.loadtxt(input_fn)
            except:
                input_array = np.loadtxt(input_fn, delimiter=',')
        else:
            assert (type(input_fn_or_array) == type(np.array([]))), "ERROR: Wrong type %s [must be a file name or np array]" % (type(input_fn_or_array))
            input_array = input_fn_or_array

#        if cell_type == 'mit':
#            input_array = self.remove_silent_mitral_cells(input_array)
        distances = self.calculate_distances_from_normalized_output(input_array, cell_type='mit')

        dist_mat = Orange.core.SymMatrix(distances)

        # without initial guess
        n_dim = self.params['n_dim_mds']
        mds = Orange.projection.mds.MDS(dist_mat, dim=n_dim)
        # Optimization loop; calculate the stress only after each 10 optimization steps:
        for i in xrange(n_steps):
            old_stress = mds.avg_stress
            for j in xrange(20):
#                mds.SMACOFstep()
                mds.smacof_step()
            mds.calc_stress(Orange.projection.mds.KruskalStress)
            print "i: %d stress: %.3e" % (i, mds.avg_stress)
            if ((old_stress * thresh) > np.abs(old_stress - mds.avg_stress)):
                print "MDS stops after iteration %d because average_stress doesn't change more than %.3e percent" % (i, thresh * 100.)
                break

        # with initial guess
        #mds = Orange.projection.mds.MDS(dist_mat, dim=n_dim, points=init_guess_list)

        #print "run mds"
        #n_runs = 10000
        #mds.run(n_runs)
        # calculate the stress value
        print "stress", mds.avg_stress
        #, mds.calc_stress(Orange.projection.mds.KruskalStress)

        #print mds.points
        points = []
        for i in xrange(len(mds.points)):
            x = mds.points[i][0]
            y = mds.points[i][1]
            z = mds.points[i][2]
            points.append((x, y, z))
        #    pylab.scatter(x, y)#, c=colors[i % (len(colors))])

        mds_output = np.array(points)
        print "saving mds output to:", mds_output_fn
        np.savetxt(mds_output_fn, mds_output)
        return mds_output


#    def remove_silent_mitral_cells(self, input_array):
#        if not os.path.exists(self.params['silent_mit_fn']): don't do anything
#            return input_array
#        else: how to avoid indexing inconsistencies in the input_array used for the mds ?!?



    def vq_ob_oc(self, points_fn, overlap=0, thresh=1e-6):

        print "loading mds from:", points_fn
        points = np.loadtxt(points_fn)

        n_mit = self.params['n_mit']
        n_hc = self.params['n_hc']
        n_mc = self.params['n_mc']
        #n_clusters = n_hc * n_mc
        n_clusters = n_hc

        d = scvq.whiten(points) # gives unit variance, necessary before using scipy's kmeans 
        centroids, distortions = scvq.kmeans(d, n_clusters, thresh=thresh)
        #centroids, distortions = scvq.kmeans(d, n_clusters)
        code, dist = scvq.vq(d, centroids)
        #print "code:", code
        dist_mat = self.get_distances(points, centroids) # distance matrix between mitral cells and hypercolumns represented by centroids

        membership_array = np.zeros((n_mit, overlap+1))
        membership_dict = {} # membership_dict['hypercolumn'] = list of source mitral cells
        # membership_dict gives MIT - HC connectivity

        for i in xrange(n_hc):
            membership_dict['%d' % i] = []
        # get the overlap + 1 closest centroids for each 
        if overlap > 0:
            for mit in xrange(n_mit):
                # get the centroid indices sorted by distance 
                indices_sorted_by_dist = dist_mat[mit, :].argsort()
                closest_centroids = indices_sorted_by_dist[0: overlap + 1]
                membership_array[mit, :] = closest_centroids
                for centroid in closest_centroids:
                    membership_dict['%d' % centroid].append(mit)
#                membership_dict['%d' % mit] = closest_centroids.tolist()

        else:
#            for i in xrange(n_hc):
            for mit in xrange(n_mit):
                hc = code[mit]
                membership_dict['%d' % hc].append(mit)
                membership_array[mit, 0] = code[mit]

        mit_hc_conn_fn = self.params['mit_hc_conn_fn']
        np.savetxt(mit_hc_conn_fn, membership_array, fmt='%d')

        binary_conn_mat = np.zeros((n_mit, n_hc * n_mc)) # binary connectivity matrix
        for cluster in membership_dict.keys():
            #        for cluster in xrange(n_clusters):
            # get the members belonging to this cluster
            members = membership_dict[cluster]
            member_coords = []
            for i in xrange(len(members)):
                member_coords.append(points[members[i], :])
#                print "debug mits", members[i], member_coords[members[i], :]
            # cluster the space spanned by mitral cells projecting to this hypercolumn
            member_coords = np.array(member_coords)
            mc_coords, mc_distortions = scvq.kmeans(member_coords, n_mc)

            # get the coordinates and code book, i.e. which mit projects to which minicolumn in the hypercolumn
            mc_code, mc_dist = scvq.vq(member_coords, mc_coords)

            mc_offset = int(cluster) * n_mc
            for i in xrange(len(members)):
                mit = members[i]
                binary_conn_mat[mit, mc_code[i] + mc_offset] = 1
        output_fn = self.params['abstract_binary_conn_mat_fn'] 
        print "output_fn:", output_fn
        np.savetxt(output_fn, binary_conn_mat)



    def vq(self, input_fn_or_array, output_fn, n_clusters_or_guess, overlap=0, thresh=1e-3, show=False, remove_silent_cells_fn=False):
        """
        input_fn_or_array: space to be clustered
        output_fn: target file storing the output_mask the code book as array
        n_clusters_or_guess: 
            is either the target number of clusters to be put in the input space
            or 
            it is a guess of initial centroids (usefull when clustering the mitral cell activity space among the minicolumns, because then kmeans crashes ....)
        remove_silent_cells_fn: filename storing the indices of points (in the MDS space) which are to be ignored (e.g. because did not have any activity so invalid MI-MDS placement)
            this needs to be done, because otherwise scvq.whiten(points) will give NaN for silent cells (since their output is always zero, so variance is NaN) and make the kmeans and VQ crash

        Writes the output mask (projection matrix) to output_fn
        """

        # Process the input parameters
        if (type(input_fn_or_array) == type('')):
            input_fn = input_fn_or_array
            print "VQ: loading coordinates from:", input_fn
            points = np.loadtxt(input_fn)
        else: # assert it's in array
            assert (type(input_fn_or_array) == type(np.array([]))), "ERROR: Wrong type %s [must be a file name or np array]" % (type(input_fn_or_array))
            points = input_fn_or_array

        n_points = points[:, 0].size
        d = scvq.whiten(points) # gives unit variance, necessary before using scipy's kmeans 

        # avoid empty clusters, i.e. no target unit should be without a source
        n_empty_hc = 1
        kmeans_trial = 0
        while (n_empty_hc != 0) and (kmeans_trial < 1000):
            print 'Running kmeans ...'
            if (type(n_clusters_or_guess) == type(np.array([]))):
                n_clusters = n_clusters_or_guess[:, 0].size
                guessed_centroids = n_clusters_or_guess
                centroids, distortions = scvq.kmeans2(d, guessed_centroids, minit='matrix')
            else:
                assert (type(n_clusters_or_guess) == type(0)), "MDQVQ.vq got wrong type for n_clusters_or_guess: %s\nNeed either coordinates or an integer" % (type(n_clusters_or_guess))
                n_clusters = n_clusters_or_guess
                centroids, distortions = scvq.kmeans(d, n_clusters, thresh=thresh)

            print 'Running vq...', kmeans_trial
            code, dist = scvq.vq(d, centroids)
            
#            print "code:", code
#            print "dist:", dist
            n_empty_hc = 0
            for hc in xrange(self.params['n_hc']):
                # take into account those units, that did not show activity (silent_mits) but still exist in the MI-space
                count_projections_to_tgt = 0
                units_projecting_to_hc = (hc == code).nonzero()[0]
                if os.path.exists(remove_silent_cells_fn):
                    silent_mits = np.loadtxt(remove_silent_cells_fn).tolist()
#                    print 'debug units_projecting_to_hc', units_projecting_to_hc
#                    print 'debug silent_mits', silent_mits, type(silent_mits), len(silent_mits)
                    if type(silent_mits) == type(.1): # it's only one silent mitral cell
                        silent_mits = [silent_mits]
                    n_units_without_silent = len(set(units_projecting_to_hc).difference(silent_mits))
                else:
                    n_units_without_silent = units_projecting_to_hc.size
                print 'n_units projecting to hc without silent ones %d (without overlap): %d' % (hc, n_units_without_silent)

                if n_units_without_silent == 0:
                    n_empty_hc += 1
            kmeans_trial += 1

        # add the MTs of the next 
        print 'MDSVQ: kmeans_trials to get non-empty targets %d ' % (kmeans_trial)
        if overlap > 0:
            # for each point get the 'overlap' closest centroids
            nn = self.get_nearest_neighbors(points, centroids, overlap)
            # nn[point, :] = indices of the centroids being neares to point
            output_mask = np.zeros((n_points, n_clusters))
            for i in xrange(n_points):
                tgt_indices = nn[i, :]
                output_mask[i, code[i]] = 1 # set the connection to the centroid assigned to the source point
                output_mask[i, tgt_indices] = 1 # set the connection between the source point and the nearest neighbour centroids
        else: 
            output_mask = np.zeros((n_points, n_clusters))
            for i in xrange(n_points):
                output_mask[i, code[i]] = 1

        for hc in xrange(self.params['n_hc']):
            n_src_to_hc = output_mask[:, hc].nonzero()[0].size
            print 'n_units projecting to hc %d (with overlap, not considering if MT is silent or not): %d' % (hc, n_src_to_hc)
            assert (n_src_to_hc > 0)


        # check for each centroid (e.g. target HC), that it gets at least one vector (source)
        for c_ in xrange(output_mask[0, :].size):
#            assert (output_mask[:, c_].sum() > 0), 'ERROR: MDSVQ.vq resulted in an empty centroid: %d' % c_
            if (output_mask[:, c_].sum() == 0): 
                print 'ERROR: MDSVQ.vq resulted in an empty centroid: %d' % c_
        print "VQ results: ", output_fn
        np.savetxt(output_fn, output_mask)

        # plotting
        fig = pylab.figure()
        ax = Axes3D(fig)

        colored = True
        # map centroid numbers to different colors
        if (colored):
            min_4d = np.min(code)
            max_4d = np.max(code)
            range_4d = float(max_4d - min_4d)
            norm = matplotlib.mpl.colors.Normalize(vmin=min_4d, vmax=max_4d)
            m = matplotlib.cm.ScalarMappable(norm=norm, cmap=matplotlib.cm.jet)
            m.set_array(np.arange(min_4d, max_4d, 0.01))
            colors = m.to_rgba(code)
            cax = ax.scatter(d[:,0], d[:,1], d[:,2], c=colors, marker='o', linewidth='5', edgecolor=colors)

        np.savetxt(output_fn, output_mask, delimiter='\t')
        ax.set_xlabel('X Label')
        ax.set_ylabel('Y Label')
        ax.set_zlabel('Z Label')
        output_fig = self.params['figure_folder'] + '/ob_oc_mds_vq.png'
        print 'Saving figure to:', output_fig
        pylab.savefig(output_fig)
        if show:
            pylab.show()
        return code, dist, centroids


    def get_nearest_neighbors(self, points, centroids, overlap):
        """
        for each point find the overlap closest centroids
        """

        n_points = points[:, 0].size
        n_centroids = centroids[:, 0].size
        distances = np.zeros((n_points, n_centroids))
        dist_mat = self.get_distances(points, centroids) # get distance matrix between points and centroids
        nearest_neighbors_indices = np.zeros((n_points, overlap), dtype=np.int)

        for p in xrange(n_points):
            indices_sorted_by_dist = dist_mat[p, :].argsort()
            nearest_neighbors_indices[p, :] = indices_sorted_by_dist[:overlap]

        return nearest_neighbors_indices


    def calculate_distances_from_normalized_output(self, normed_activity, cell_type='mit'):
        """
        Calculates the mutual information based on the normed_activity
        """

        n_units = normed_activity[0, :].size # number of columns in the activity file
        n_patterns = self.params['n_patterns']

        prob_i = np.zeros(n_units)
        prob_ij = np.zeros((n_units, n_units))

        for i in xrange(n_units):
            prob_i[i] = normed_activity[:, i].sum() / n_patterns

        for i in xrange(n_units):
            for j in xrange(n_units):
                joint_prob = 0.
                for pn in xrange(n_patterns):
                    joint_prob += normed_activity[pn, i] * normed_activity[pn, j]
                prob_ij[i, j] = joint_prob / n_patterns

        print "build the distance matrix, calculate the mutual information between units .... "
        # build the distance matrix
        # calculate the mutual information between units
        mi = np.zeros((n_units, n_units))
        joint_entropy = np.zeros((n_units, n_units))
        distances = np.zeros((n_units, n_units))
        for i in xrange(n_units):
            for j in xrange(n_units):
                product = prob_i[i] * prob_i[j]
                if ((product != 0) and (prob_ij[i, j] != 0.)):
                    mi[i, j] = prob_ij[i, j] * np.log(prob_ij[i, j] / product)
                else:
                    mi[i, j] = 0.0

                if (prob_ij[i, j] != 0.):
                    joint_entropy[i, j] = -1.0 * prob_ij[i, j] * np.log(prob_ij[i, j])
                else:
                    joint_entropy[i, j] = 0.

                if (joint_entropy[i, j] != 0.):
                    distances[i, j] = 1 - mi[i, j] / joint_entropy[i, j]
                else:
                    distances[i, j] = 1.

        mi_output_fn = self.params['mutual_information_fn_base'] + cell_type + '.dat'
        je_output_fn = self.params['joint_entropy_fn_base'] + cell_type + '.dat'
        dist_output_fn = self.params['information_distance_fn_base'] + cell_type + '.dat'

        print "saving to", mi_output_fn, je_output_fn, dist_output_fn
        np.savetxt(mi_output_fn, mi)
        np.savetxt(je_output_fn, joint_entropy)
        np.savetxt(dist_output_fn, distances)

        return distances

    def create_mitral_response_space(self, input_fn, activity_fn, output_fn, remove_silent_cells_fn=False, optional_mds=False, mit_mc_kmeans_trial=0):
        """
            input_fn contains the mask of mitral cells (rows) projecting to the respective hypercolumn (column), i.e. the result of the first vector quantization in the mutual information space
            # 4) For each hypercolumn, create a new space spanned by the mitral cells projecting to the hc
            #   Each mitral cell represents one dimension and each pattern represents one vector or point in that space.
            #   The value of one component in such a vector is equal to the normalized activation of the respective mitral cell.
            #   The n_patterns vectors are clustered by VQ among the minicolumns in the hypercolumn.

            Should create the output file self.params['binary_oc_activation_fn'] = '%s/binary_oc_activation.dat' % (self.params['other_folder']) 
        """
        mask = np.loadtxt(input_fn)
        try:
            ob_activity = np.loadtxt(activity_fn)
        except:
            ob_activity = np.loadtxt(activity_fn, delimiter=',')
        n_mc = self.params['n_mc']
        n_hc = mask[0, :].size
        n_mit = mask[:, 0].size
        n_patterns = self.params['n_patterns']
        n_mc_oc = self.params['n_hc'] * self.params['n_mc']
        assert (n_hc == self.params['n_hc']), "Mismatch between parameters in network_parameters_Cray and files used here, e.g.: %s" % (input_fn)
        output_array = np.zeros((n_patterns, n_mc_oc))

        if (remove_silent_cells_fn):
            silent_mits = np.loadtxt(self.params['silent_mit_fn'])
            if silent_mits.size == 1:
                silent_mits = silent_mits.tolist() # this returns a float for some reason ...
                silent_mits = [silent_mits]
        else:
            silent_mits = []

        n_src_mit = np.zeros(n_hc)
        for hc in xrange(n_hc):
            src_mit = np.nonzero(mask[:, hc])[0].tolist()
            to_remove = set(silent_mits).intersection(set(src_mit))
            print "\nRunning VQ for HC: ", hc, " n_mit_src:", len(src_mit)
            print "Silent mitral cells need to be ignored:", to_remove
            for silent_mit in to_remove:
                src_mit.remove(silent_mit)

            activity_space = np.zeros((n_patterns, len(src_mit)))
            n_src_mit[hc] = len(src_mit)
            for pattern in xrange(n_patterns):
                # mask only the mitral cells projecting to the hc (src_mits)
                activity_space[pattern, :] = ob_activity[pattern, :].take(src_mit)

            if (optional_mds):
                mit_coords_fn = self.params['mit_response_space_fn_base']
                self.mds(activity_space, mit_coords_fn, cell_type='mc')

            # now run VQ in that activity space
            vq_output_fn = self.params['folder_name'] + '/' + 'tmp_vq_hc_%d.dat' % hc

            # instead of giving n_mc as the number of target clusters to find in the activity space
            # one can give a guess of initial centroids to scipy.cluster.vq.kmeans 
            guessed_centroids = np.random.rand(n_mc, len(src_mit))

            d = scvq.whiten(activity_space) # gives unit variance, necessary before using scipy's kmeans 

            # NOTE: it makes a big difference if you take random initial guessed centroids or just take k=n_mc
#            centroids, distortions = scvq.kmeans(d, guessed_centroids, thresh=1e-8)
#            centroids, distortions = scvq.kmeans(d, n_mc, thresh=1e-8)
#            centroids, distortions = scvq.kmeans2(d, guessed_centroids, minit='matrix')
            # Distortion is defined as the sum of the squared differences between the observations and the corresponding centroid.
            centroids, distortions = scvq.kmeans2(d, n_mc, minit='points')
            codes, dist = scvq.vq(d, centroids)

            dist_to_write = '\t%.2f\t%.2f\t%.2f\n' % (dist.mean(), dist.std(), dist.sum())
            dist_fn = self.params['mit_mc_vq_distortion_fn'] + '%d.dat' % mit_mc_kmeans_trial
            distortion_file = open(dist_fn, 'a')
            distortion_file.write(dist_to_write)

            if (optional_mds):
                mit_coords_fn = self.params['mit_response_space_fn_base']
                centroids_coords = self.params['mit_response_space_centroids_fn_base']

            assert (len(codes) == n_patterns), "ERROR: Something went wrong when clustering the activity space for hc %d, maybe some MC got no input?" % (hc)
            for pn in xrange(len(codes)):
#                print "codes[%d] = %d" % (pn, codes[pn])
                mc = hc * n_mc + codes[pn]
                output_array[pn, mc] = 1


        print 'Debug MDSVQ: n_src_mit = %.2f +- %.2f' % (n_src_mit.mean(), n_src_mit.std())
        print "MDSVQ.create_mitral_response_space writes to:", output_fn, '\n\t', dist_fn
        np.savetxt(output_fn, output_array)



    def get_distances(self, x, y):
        """
        x: array of list of points
        y: array of list of points (same dimensionality as x)
        """
        n_x = len(x)
        n_y = len(y)
        dist_matrix = np.zeros((n_x, n_y))
        dist_function = euclidean_dist
        for i in xrange(n_x):
            for j in xrange(n_y):
                dist_matrix[i, j] = euclidean_dist(x[i], y[j])
        return dist_matrix
