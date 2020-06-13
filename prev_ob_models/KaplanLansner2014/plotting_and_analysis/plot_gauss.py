
import numpy as np
import matplotlib.mlab as mlab
import pylab


p = [1.631787e+02, 6.670855e+00, 1.977871e+00, \
     1.909487e+01, 1.110809e+01, 3.353855e+00, \
     4.188897e+00, 4.088460e+01, 4.966478e-01]

def trimodal_gauss(x, p):
    w1 = p[0]
    mu1 = p[1]
    sigma1 = p[2]
    w2 = p[3]
    mu2 = p[4]
    sigma2 = p[5]
    w3 = p[6]
    mu3 = p[7]
    sigma3 = p[8]

    return w1 * mlab.normpdf(x, mu1, sigma1) + w2 * mlab.normpdf(x, mu2, sigma2) + w3 * mlab.normpdf(x, mu3, sigma3)

def gauss(x, mu, sigma):
    return 1. / (sigma * np.sqrt(2 * np.pi)) * np.exp( - (x - mu)**2 / (2. * sigma ** 2))


bins = np.linspace(0, 4.330310991999920844e+01, 100)
#G1 = p[0] * gauss(bins, p[1], p[2])
G1 = p[0] * mlab.normpdf(bins, p[1], p[2])
print 'G1', G1
G2 = p[3] * gauss(bins, p[4], p[5])
print 'G2', G2
G21 = p[3] * mlab.normpdf(bins, p[4], p[5])
print 'G2 normpdf', G21 - G2
G3 = p[6] * gauss(bins, p[7], p[8])
print 'G3', G3
pylab.plot(bins, G1, c='k', lw=3)
pylab.plot(bins, G2, c='b', lw=3)
pylab.plot(bins, G3, c='g', lw=3)

pylab.xlim((bins[0], bins[-1]))
#pylab.ylim((0, 1))
pylab.show()

