# -*- coding: utf-8 -*-
"""
Created on Sun Oct 17 22:10:35 2010

@author: -
"""
import numpy
from enthought.mayavi import mlab

x, y = numpy.mgrid[0:3:1,0:3:1]
s = mlab.surf(x, y, numpy.asarray(x*0.1, 'd'))
#
for i in range(10):
    s.mlab_source.scalars = numpy.asarray(x*0.1*(i+1), 'd')
#
# Produce some nice data.
n_mer, n_long = 6, 11
pi = numpy.pi
dphi = pi/1000.0
phi = numpy.arange(0.0, 2*pi + 0.5*dphi, dphi, 'd')
mu = phi*n_mer
x = numpy.cos(mu)*(1+numpy.cos(n_long*mu/n_mer)*0.5)
y = numpy.sin(mu)*(1+numpy.cos(n_long*mu/n_mer)*0.5)
z = numpy.sin(n_long*mu/n_mer)*0.5

# View it.
l = mlab.plot3d(x, y, z, numpy.sin(mu), tube_radius=0.025, colormap='Spectral')

# Now animate the data.
ms = l.mlab_source
for i in range(10):
    x = numpy.cos(mu)*(1+numpy.cos(n_long*mu/n_mer +
                                    numpy.pi*(i+1)/5.)*0.5)
    scalars = numpy.sin(mu + numpy.pi*(i+1)/5)
    ms.set(x=x, scalars=scalars)
