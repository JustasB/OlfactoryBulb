'''
Load balance requires an at present unknown gid distribution that cannot
be calculated til the connections are known. In particular, the complexity
of granules are dominated by their number of MGRS.
Calculation of all the connections and complexities can be accomplished
in parallel if we temporarily use a whole cell gid distribution in which
rank is easily derivable from gid, e.g. rank = gid%nhost. Then it is easy
to communicate the information needed to each rank.
'''

import params
from common import *
from all2all import all2all
import util

t_begin = h.startsw()

def gid2rank(gid):
  return gid%nhost

#set of gids on this rank (default round-robin)
def round_robin_distrib(model):
  model.gids = set(range(rank, ncell, nhost))
  model.mitral_gids = set(range(rank, nmitral, nhost))
  model.granule_gids = model.gids - model.mitral_gids

round_robin_distrib(getmodel())

'''
In this section, presume connections determined by m2g_connections.py.
I.e. mitral statistics controlled and cause unknown granule statistics.
'''
import mkmitral

def mk_mitrals(model):
  ''' Create all the mitrals specified by mitral_gids set.'''
  model.mitrals = {}
  for gid in model.mitral_gids:
    m = mkmitral.mkmitral(gid)
    model.mitrals.update({gid : m})
  util.elapsed('%d mitrals created and connections to mitrals determined'%int(pc.allreduce(len(model.mitrals),1)))

def mk_gconnection_info_part1(model):
  ''' after mk_gconnection_info_part2()
      rank_gconnections is the connection info for granules on rank ggid%nhost
      also granule_gids are the granules on this rank (granules with no
      connection will not exist)
  '''
  model.rank_gconnections = {}
  for cilist in model.mconnections.values():
    for ci in cilist:
      ggid = ci[3]
      r = gid2rank(ggid)
      if not model.rank_gconnections.has_key(r):
        model.rank_gconnections.update({r : []})
      model.rank_gconnections[r].append(ci)

def mk_gconnection_info_part2(model):
  #transfer the gconnection info to the proper rank and make granule_gids set
  model.rank_gconnections = all2all(model.rank_gconnections)
  util.elapsed('rank_gconnections known')
  model.granule_gids = set([i[3] for r in model.rank_gconnections for i in model.rank_gconnections[r]])
  util.elapsed('granule gids known on each rank')

def mk_gconnection_info(model):
  mk_gconnection_info_part1(model)
  mk_gconnection_info_part2(model)
  util.elapsed('mk_gconnection_info (#granules = %d)'%int(pc.allreduce(len(model.granule_gids),1)))


if __name__ == '__main__':
  model = getmodel()
  mk_mitrals(model)
  mk_mconnection_info(model)
  mk_gconnection_info_part1(model)

  sizes = all2all(model.rank_gconnections, -1)
  for r in util.serialize():
    print rank, " all2all sizes ", sizes

if rank == 0: print "determine_connections ", h.startsw()-t_begin

