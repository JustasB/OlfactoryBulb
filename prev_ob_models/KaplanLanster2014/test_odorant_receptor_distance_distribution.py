import numpy as np
import matplotlib.mlab as mlab
import pylab

def distance_generator():

#    p = [1.631787e+02, 6.670855e+00, 1.977871e+00, \
#         1.909487e+01, 1.110809e+01, 3.353855e+00, \
#         4.188897e+00, 4.088460e+01, 4.966478e-01]

#    p = [162.310869565, 6.67080434783, 1.98630434783,\
#          19.8056521739, 10.8089130435, 3.32682608696, \
#           4.4382173913, 40.8932608696, 0.456293478261] # these values are taken from clustering the odorant space with 40 ORs

    # fit params obtained from the averaged (normalized by nOR) distance distribution:
    p = [ 162.796725363 ,  6.67708323673 ,  1.98973030466 ,  19.3388197446 ,  10.853527839 ,  3.3441114626 ,  4.43148000439 ,  40.897695747 ,  0.461489675411]
    w1 = p[0]
    mu1 = p[1]
    sigma1 = p[2]
    w2 = p[3]
    mu2 = p[4]
    sigma2 = p[5]
    w3 = p[6]
    mu3 = p[7]
    sigma3 = p[8]
    which_gauss = np.random.uniform(0, 1.)
    if which_gauss < p1:
#        print 'G1 ', 
        return max(0, np.random.normal(mu1, sigma1))
    elif (which_gauss < p2 + p1):
#        print 'G2 ', 
        return max(0, np.random.normal(mu2, sigma2))
    elif (which_gauss < p3 + p2 + p1):
#        print 'G3 ', 
        return max(0, np.random.normal(mu3, sigma3))
    else:
        print '\nERROR\n'


def get_prob_for_distributions():
    """
    Based on the integral of the three normal distributions,
    the likelihood from which of the three distributions a distance is to be drawn 
    is calculated here.
    Returns the three probabilities for the three distributions.
    """
    p = [1.631787e+02, 6.670855e+00, 1.977871e+00, \
         1.909487e+01, 1.110809e+01, 3.353855e+00, \
         4.188897e+00, 4.088460e+01, 4.966478e-01]
    w1 = p[0]
    mu1 = p[1]
    sigma1 = p[2]
    w2 = p[3]
    mu2 = p[4]
    sigma2 = p[5]
    w3 = p[6]
    mu3 = p[7]
    sigma3 = p[8]

    dist_range = (0, 4.330310991999920844e+01)
    x = np.linspace(dist_range[0], dist_range[1], 1000)
    A1 = np.array(w1 * mlab.normpdf(x, mu1, sigma1)).sum()
    A2 = np.array(w2 * mlab.normpdf(x, mu2, sigma2)).sum()
    A3 = np.array(w3 * mlab.normpdf(x, mu3, sigma3)).sum()
    print 'debug A', A1, A2, A3
    p1 = A1 / (A1 + A2 + A3)
    p2 = A2 / (A1 + A2 + A3)
    p3 = A3 / (A1 + A2 + A3)
    return p1, p2, p3


def transform_dist_to_affinity_exp(dist, alpha):
    return np.exp(-(dist)**2 / alpha)

#    return np.exp(-dist * alpha)

def transform_dist_to_affinity(dist):
    return 1. / (1. + dist)


p1, p2, p3 = get_prob_for_distributions()
print 'Debug p1 p2 p3', p1, p2, p3, p1 + p2 + p3
np.random.seed(0)

n_or = 50
n_pattern = 1000
n_samples = n_or * n_pattern
dist_samples = np.zeros(n_samples)
for i_ in xrange(n_samples):
    dist_samples[i_] = distance_generator()


pylab.rcParams.update({'figure.subplot.hspace' : .4})
fig = pylab.figure()
ax1 = pylab.subplot(211)
ax2 = pylab.subplot(212)

count_dist, bins_dist = np.histogram(dist_samples, bins=200)
count_norm = count_dist / float(n_samples)
#expected_value_dist = np.sum(count_norm * bins_dist[:-1])
expected_value_dist = dist_samples.mean()
ax1.bar(bins_dist[:-1], count_norm, width=bins_dist[1]-bins_dist[0])
ax1.set_title('Distance distribution')
ax1.set_xlabel('Distances')
ax2.set_ylabel('Normalized count')
#print 'debug counts', count_dist, n_samples
#print 'debug normalized counts', count_norm
#print 'bins_dist * count_norm', bins_dist[:-1] * count_norm, (bins_dist[:-1] * count_norm).sum()
print 'Expected value dist:', expected_value_dist

alpha = 1 * expected_value_dist**2
affinities = transform_dist_to_affinity_exp(dist_samples, alpha)
#affinities = transform_dist_to_affinity(dist_samples)

count_aff, bins_aff = np.histogram(affinities, bins=200)
count_norm = count_aff / float(n_samples)
#expected_value_aff = np.sum(count_norm * bins_aff[:-1])
expected_value_aff = affinities.mean()
ax2.bar(bins_aff[:-1], count_norm, width=bins_aff[1]-bins_aff[0])
ax2.set_title('Affinity distribution')
ax2.set_xlabel('Affinities')
ax2.set_ylabel('Normalized count')
#print 'debug counts', count_aff, n_samples
#print 'debug normalized counts', count_norm
#print 'bins_aff * count_norm', bins_aff[:-1] * count_norm, (bins_aff[:-1] * count_norm).sum()
print 'Expected value affinities:', expected_value_aff


pylab.show()

