"""
1) cluster_odorant_space.py
2) average_OR_affinity_distributions.py

--> expect the averaged distance distributions to be in one folder 
e.g. OR_placement/DistanceDistributions/

"""

import pylab
import numpy as np
import sys
import average_OR_affinity_distributions as AVG
from scipy.optimize import leastsq
import matplotlib.mlab as mlab

def transform_to_affinity(x, sigma):
    y = np.exp(- x**2 / sigma**2)
    return y

def distance_generator(p, p1, p2, p3):
    """
    p -- fit params for the three normal distributions
    # fit params obtained from the averaged (normalized by nOR) distance distribution:
    p = [ 162.796725363 ,  6.67708323673 ,  1.98973030466 ,  19.3388197446 ,  10.853527839 ,  3.3441114626 ,  4.43148000439 ,  40.897695747 ,  0.461489675411]

    p1, p2, p3 -- the probabilities with which the respective normal distribution is chosen
    """

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
        return max(0, np.random.normal(mu1, sigma1))
    elif (which_gauss < p2 + p1):
        return max(0, np.random.normal(mu2, sigma2))
    elif (which_gauss < p3 + p2 + p1):
        return max(0, np.random.normal(mu3, sigma3))
    else:
        print '\nERROR\n'


def get_prob_for_distributions(p):
    """
    Based on the integral of the three normal distributions,
    the likelihood from which of the three distributions a distance is to be drawn 
    is calculated here.
    Returns the three probabilities for the three distributions.
    """
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
    p1 = A1 / (A1 + A2 + A3)
    p2 = A2 / (A1 + A2 + A3)
    p3 = A3 / (A1 + A2 + A3)
    return p1, p2, p3



data_folder = 'OR_placement/DistanceDistributions/'
fn_base = 'mean_distance_distribution_and_fit_params_' 

n_ORs = [20, 66]

fn = data_folder + fn_base + str(n_ORs[0]) + '.txt'
d = np.loadtxt(fn)
n_bins = d[:, 0].size
all_distance_distributions = np.zeros((n_bins, n_ORs[1] - n_ORs[0]))
bins = d[:, 0].tolist() # x, axis = the distances
bincenters = 0.5*(np.array(bins[1:]) + np.array(bins[:-1]))
x = bincenters

for i_, n_OR in enumerate(range(n_ORs[0], n_ORs[1])):
    fn = data_folder + fn_base + str(n_OR) + '.txt'
    d = np.loadtxt(fn)
    all_distance_distributions[:, i_] = d[:, 1]

fn = 'OR_placement/all_distance_distributions_ybins_xORs.dat'
print 'Saving all distance distributions into one file:', fn
np.savetxt(fn, all_distance_distributions)

mean_distance_distribution = np.zeros((n_bins, 2))
mean_affinity = np.zeros(n_bins)
for i_bin in xrange(n_bins):
    mean_distance_distribution[i_bin, 0] = all_distance_distributions[i_bin, :].mean()
    mean_distance_distribution[i_bin, 1] = all_distance_distributions[i_bin, :].std()

distance_probabilities = mean_distance_distribution[:, 0] / np.sum(mean_distance_distribution[:, 0])
expected_distance = np.sum(bins * distance_probabilities)
print 'Expected distance value:', expected_distance
mean_affinity = transform_to_affinity(all_distance_distributions[:, 0], expected_distance)

guess_params_trimodalgauss = [165., 7., 2., 20., 11., 3., 3., 40., .5]
opt_params = leastsq(AVG.residuals_trimodal_gauss, guess_params_trimodalgauss, args=(mean_distance_distribution[:-1, 0], x))[0]
#opt_params = leastsq(AVG.residuals_trimodal_gauss, guess_params_trimodalgauss, args=(mean_count[x0:x1], x))[0]
opt_fit = AVG.peval_trimodal_gauss(x, opt_params)
reduced_chi_square = AVG.get_reduced_chi_square(mean_distance_distribution[:-1, 0], opt_fit, mean_distance_distribution[:-1, 1], len(opt_params))
print 'Reduced chi square:', reduced_chi_square

plot_params = {'backend': 'png',
              'axes.labelsize': 20,
              'axes.titlesize': 20,
              'text.fontsize': 20,
              'xtick.labelsize': 16,
              'ytick.labelsize': 16,
              'legend.pad': 0.2,     # empty space around the legend box
              'legend.fontsize': 14,
               'lines.markersize': 0,
               'lines.linewidth': 3,
              'font.size': 12,
              'path.simplify': False}
#              'figure.subplot.left':.10,
#              'figure.subplot.bottom':.13,
#              'figure.subplot.right':.94,
#              'figure.subplot.top':.94}
#              'figure.figsize': get_fig_size(800)}

pylab.rcParams.update(plot_params)
fig = pylab.figure()
ax = fig.add_subplot(111)
ax.bar(bins, mean_distance_distribution[:, 0], yerr=mean_distance_distribution[:, 1], width=bins[1]-bins[0])
label = 'Tri-modal normal distribition\n \
        p1=(%.1f, %.1f, %.1f),\n \
        p2=(%.1f, %.1f, %.1f),\n \
        p3=(%.1f, %.1f, %.1f)' % \
        (opt_params[0], opt_params[1], opt_params[2], \
        opt_params[3], opt_params[4], opt_params[5], \
        opt_params[6], opt_params[7], opt_params[8])
print 'label', label
print 'fit params:', '[',
for i_ in xrange(9):
    print opt_params[i_], ', ',
print ']'
opt_params[0], opt_params[1], opt_params[2], \
        opt_params[3], opt_params[4], opt_params[5], \
        opt_params[6], opt_params[7], opt_params[8]
lw = 5
ax.plot(x, opt_fit, 'r--', label=label, lw=lw)
ax.plot(x, AVG.peval_gauss(x, [abs(opt_params[0]), abs(opt_params[1]), abs(opt_params[2])]), ls='--', c='#00FFFF', lw=lw)
ax.plot(x, AVG.peval_gauss(x, [abs(opt_params[3]), abs(opt_params[4]), abs(opt_params[5])]), 'y--', lw=lw)
ax.plot(x, AVG.peval_gauss(x, [abs(opt_params[6]), abs(opt_params[7]), abs(opt_params[8])]), 'k--', lw=lw)
ax.set_ylim((0, ax.get_ylim()[1]))

xlabel = 'Distance between vOR and odorant'
ax.set_xlabel(xlabel)
ylabel = 'Normalized number of occurrence'
ax.set_ylabel(ylabel)

output_fig = data_folder + 'averaged_distance_distribution.png'
print 'Saving fig to:', output_fig
fig.savefig(output_fig, dpi=200)
output_fig = data_folder + 'averaged_distance_distribution.pdf'
print 'Saving fig to:', output_fig
fig.savefig(output_fig, dpi=200)

# based on the parameters for the normal distributions, one can
# obtain with the sample probability for each of the three
p1, p2, p3 = get_prob_for_distributions(opt_params[:9])
# sample from the fitted distribution
n_or = 500
n_pattern = 1000
n_samples = n_or * n_pattern
dist_samples = np.zeros(n_samples)
affinities = np.zeros(n_samples)
for i_ in xrange(n_samples):
    dist_samples[i_] = distance_generator(opt_params[:9], p1, p2, p3)
    affinities[i_] = transform_to_affinity(dist_samples[i_], expected_distance)

count_aff, bins_aff = np.histogram(affinities, bins=100)
count_norm = count_aff / float(n_samples)
#expected_value_aff = np.sum(count_norm * bins_aff[:-1])
expected_value_aff = affinities.mean()


plot_params = {'backend': 'png',
              'axes.labelsize': 20,
              'axes.titlesize': 20,
              'text.fontsize': 20,
              'xtick.labelsize': 16,
              'ytick.labelsize': 16,
              'legend.pad': 0.2,     # empty space around the legend box
              'legend.fontsize': 14,
               'lines.markersize': 0,
               'lines.linewidth': 3,
              'font.size': 12,
              'path.simplify': False,
              'figure.subplot.left':.14}
#              'figure.subplot.bottom':.13,
#              'figure.subplot.right':.94,
#              'figure.subplot.top':.94}
#              'figure.figsize': get_fig_size(800)}
pylab.rcParams.update(plot_params)
fig2 = pylab.figure()
ax2 = fig2.add_subplot(111)
ax2.bar(bins_aff[:-1], count_norm, width=bins_aff[1]-bins_aff[0])
ax2.set_xlabel('Affinity between odorant and vOR')
ax2.set_ylabel('Probability')
print 'Expected value affinities:', expected_value_aff
output_fig = data_folder + 'affinity_distribution.png'
print 'Saving fig to:', output_fig
fig2.savefig(output_fig, dpi=200)
output_fig = data_folder + 'affinity_distribution.pdf'
print 'Saving fig to:', output_fig
fig2.savefig(output_fig, dpi=200)

pylab.show()


