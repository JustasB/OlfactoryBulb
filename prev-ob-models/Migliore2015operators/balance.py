from util import *
from lpt import lpt, statistics
from all2all import all2all

def load_bal(cx, npart):
  ''' cx is list of (complexity, gid) pairs on this process
      return: an LPT balanced list of gids that should belong to this process
  '''
  elapse = h.startsw()
  #send to rank 0
  r = all2all({0:cx})
  # make a list of all the (cx, gid)
  s = {}
  if rank == 0:
    c = []
    for i in r.values():
      c += i
    del r
    #distribute by LPT
    parts = lpt(c, npart)
    print statistics(parts)
    for i,p in enumerate(parts):
      s.update({i : p[1]})
  else:
    del r
  #send each partition to the proper rank
  local = all2all(s)
  del s
  if rank == 0:
    print "load_bal time %g" % (h.startsw()-elapse)
  return local[0]

if __name__ == '__main__':
  from util import serialize, finish
  if True:
    cx = [(10*rank+i, 10*rank+i) for i in range(1,5)]
    print cx
    cx = load_bal(cx, nhost)
    for r in serialize():
      print 'rank %d '%rank, cx
  if nhost > 0:
    finish()

