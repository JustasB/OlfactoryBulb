"""
This script puts a variable number of clusters in the 32-dimensional
data given in  Haddad 2008 "A metric for odorant comparison".
In there 447 odorants are described by 32 physico-chemical descriptors.
The cluster centers (centroid) are interpreted as virtual odorant receptors (ORs) and the 
distance between ORs and the centroids is stored in files.

Output:
    Each resulting distance matrix of one clustering trial is stored in a seperate file:
        dist_mat_fn = "%s/distance_matrix_OR%d_%d.dat" % (folder, n_clusters, run)

The script average_OR_affinity_distributions.py processes the results.
"""
import os
import numpy as np
import time
import scipy
import numpy.random as rnd
from scipy.cluster.vq import vq, whiten, kmeans, kmeans2
from scipy.spatial.distance import euclidean as euclidean_dist


def cluster_odor_space(num_clusters):
    global_mean = np.zeros(n_dim) # the "average odorant", needed for the F-Test
    for dim in xrange(n_dim):
        global_mean[dim] = d[:, dim].mean()

    # here it can happen, that kmeans could not place the desired number of centroids in the space
    centroids, distortions = kmeans(d, num_clusters)#, iter=100)#, thresh=1e-6)

    retrial = 0
    while (centroids.shape != (num_clusters, n_dim)):
        print "Run again trial...", retrial
        centroids, distortions = kmeans(d, num_clusters)#, iter=100)#, thresh=1e-6)
        retrial += 1
        if (retrial > 50):
            return # it's not possible to put the desired number of centroids in the space --> end of script
    code, dist = vq(d, centroids)

    # calculate the F value: F = between-group-variability / within-group-variability
    # between-group-variability = sum over i n_i * (sample_mean_for_group_i - global_mean)^2 / (K -1)
    bgv = 0.
    sample_mean = [np.zeros(n_dim) for i in xrange(num_clusters)] # a list of odorants representing the mean of that cluster
    for cluster in xrange(num_clusters):
        n_members = code.tolist().count(cluster) # how many observations belong to cluster
        # calculate the sample mean of the current cluster
        # i.e. average all odorants belonging to this cluster
        sample_size = 0
        for odorant in xrange(n_odorants):
            if (code[odorant] == cluster):
                sample_mean[cluster] += d[odorant,:]
                sample_size += 1
        sample_mean[cluster] /= sample_size
        bgv += n_members * euclidean_dist(sample_mean[cluster], global_mean) ** 2/ (num_clusters - 1)

    # within-group-variability
    wgv = 0.
    # for all clusters sum the difference between sample belonging to that cluster and the sample_mean
    for cluster in xrange(num_clusters):
        for odorant in xrange(n_odorants):
            if (code[odorant] == cluster):
                wgv += euclidean_dist(d[odorant], sample_mean[cluster]) ** 2 / (n_odorants - num_clusters)

    # total variance
    variance = 0.
    mean_odor = np.zeros(n_dim)
    for odorant in xrange(n_odorants):
        mean_odor += d[odorant, :]
    mean_odor /= n_odorants
    for odorant in xrange(n_odorants):
        variance += (euclidean_dist(d[odorant, :], mean_odor))**2
    variance /= n_odorants

    dist_matrix = np.zeros((n_odorants, num_clusters))
    for i in xrange(n_odorants):
        for c in xrange(num_clusters):
            dist_matrix[i, c] = euclidean_dist(d[i, :], centroids[c, :])

    # build affinity matrix --- this is not used in the end
    affinities = np.zeros((n_odorants, num_clusters))
    for i in xrange(n_odorants):
        for c in xrange(num_clusters):
            if (dist_matrix[i, c] == 0):
                affinities[i, c] = 1.
            else:
                affinities[i, c] = 1. / dist_matrix[i, c]


    return variance, distortions, bgv, wgv, bgv/wgv, code, dist_matrix, affinities


if __name__ == '__main__':

    t1 = time.time()
    scipy.random.seed(1)
    rnd.seed(0)

    # load data, source: http://www.nature.com/nmeth/journal/v5/n5/extref/nmeth.1197-S3.xls
    fn_raw_data = 'Haddad_data/odorant_data_first_row_is_weights.csv' 
    raw_data = np.loadtxt(fn_raw_data, delimiter=',')
    w = raw_data[0, :]
    d = raw_data[1:, :]
    n_dim = d[0, :].size

    # weigh the descriptors with the weights given in the first row
    for i in xrange(n_dim):
        d[:, i] *= w[i]
        print i, w[i]

    d = whiten(d) # necessary before using scipy's kmeans 

    n_odorants = d[:,0].size

    n_cluster_trials = 100
    folder = "OR_placement/"

    # create folders and subfolders
    if not os.path.exists(folder):
        print 'Creating directory:', folder
        os.system('mkdir %s' % (folder))
    dist_folder = '%sDistanceMatrices/' % (folder)
    aff_folder = '%sAffinityMatrices/' % (folder)
    code_folder = '%sCodeBooks/' % (folder)
    os.system('mkdir %s %s %s' % (dist_folder, aff_folder, code_folder))

    output_fn = "%s/OR_placement_evaluation.txt" % folder
    f_out = open(output_fn, 'w')
    output = ""
    output += "# n_clusters variance distortion between-cluster-variability(bgv) within-cluster-variability(wcv) F=bgv/wgv\n"
    f_out.write(output)

    for n_clusters in xrange(2, 70):
        for run in xrange(n_cluster_trials):
            print "n_clusters", n_clusters, "run: ", run
            variance, distortion, bgv, wgv, F, code, dist_matrix, affinity_matrix = cluster_odor_space(n_clusters)
            output = "%d\t%.6e\t%6e\t%.6e\t%.6e\t%.6e\n" % (n_clusters, variance, distortion, bgv, wgv, F)
            f_out.write(output)
            f_out.flush()

            code_fn = "%s/code_matrix_nOR%d_%d.dat" % (code_folder, n_clusters, run)
            dist_mat_fn = "%s/distance_matrix_OR%d_%d.dat" % (dist_folder, n_clusters, run)
            aff_mat_fn = "%s/affinity_matrix_OR%d_%d.dat" % (aff_folder, n_clusters, run)
            np.savetxt(aff_mat_fn, affinity_matrix)
            # optional:
#            np.savetxt(code_fn, code)
            np.savetxt(dist_mat_fn, dist_matrix)

    f_out.close()

    t2 = time.time()
    t = t2 - t1
    print "Time %.1f sec, %.1f min" % (t, t/60.)

