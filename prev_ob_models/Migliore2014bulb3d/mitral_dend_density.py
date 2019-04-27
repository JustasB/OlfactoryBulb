import params
from neuron import h
from neuron import numpy
import mkmitral

pc = h.ParallelContext()
nhost = int(pc.nhost())
rank = int(pc.id())

matrank = (41, 61)
domain = ((-500, 1500, 50), (-500, 2500, 50))

def f(ii):
  i = int(ii)
  density = numpy.zeros(matrank)
  for gid in range(i , params.Nmitral, nhost):
    m = mkmitral.mkmitral(gid)
    for sec in m.secdens:
      accumulate_density(sec, density, domain)
    print gid
  return density

def accumulate_density(sec, density, domain):
  sec.push()
  for i in range(int(h.n3d())):
    x,y = (h.x3d(i), h.y3d(i))
    r = (round(x, domain[0]),round(y, domain[1]))
    if not False in r:
      density[r] += 1
  h.pop_section()

def round(x, d):
  # integer toward 0 min, max, inc where min is 0, if not in interval
  #return False
  if (x < d[0] or x > d[1]):
    return False
  return int((x - d[0])/d[2])


def compute():
  for i in range(nhost):
    pc.submit(f, i)
  den = numpy.zeros(matrank)
  while(pc.working()):
    den += pc.pyret()  
  return den

if __name__ == '__main__':
  pc.runworker()
  density = compute()
  pc.done()
  print "density max = ", density.max()
  density = density * (20/density.max())

  import pickle
  pickle.dump(density, open('density.dat', 'w'))

#following works on linux if using openmpi
from mayavi.mlab import barchart,show
barchart(density)
show()

