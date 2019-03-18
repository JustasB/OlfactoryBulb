'''
Load balanced (either whole cell, or multisplit) distribution of
mitrals, and granules along with all necessary connection information
on each rank. To do this we start with a
1)round-robin distribution of mitrals and granule
2)create the mitrals
3)determine the connection info for each mitral
4)transfer the connection info appropriate for each round-robin located granule
5)create the granules
6)create the MGRS
7)determine complexity of the mitrals and granules
8)destroy the model, retaining the connection info for mitrals and granules
9)load balance the mitrals (or pieces) and granules
10)create the mitrals (multisplit pieces if needed) and granules
11)register the mitral and granule gids with their spike generation locations
11)transfer the connection info to the proper ranks
12)make the connections

Note that 1-6 is accomplished by importing the net_mitral_centric model
and executing nmc.build_net_round_robin(getmodel())
'''

from common import *
import net_mitral_centric as nmc
from complexity import determine_complexity
from balance import load_bal
from destroy_model import destroy_model
from all2all import all2all
import util
import mgrs

def whole_cell_distrib(model):
  enter = h.startsw()
  cx = determine_complexity(model)
  #start over
  destroy_model(model)
  # but we still have model.mconnections and model.rank_gconnections
  # from the round-robin distribution perspective.
  # we will need to tell the ranks  where to distribute that information

  cx = load_bal(cx, nhost) #cx is the list of (cx,gid) we want on this process
  # the new distribution of mitrals and granules
  model.gids = set([item[1] for item in cx])
  model.mitral_gids = set([gid for gid in model.gids if gid < params.gid_granule_begin])
  model.granule_gids = model.gids - model.mitral_gids

  # get the correct mconnections and gconnections
  # Round-robin ranks that presently have the connection info for the gids
  rr = {}
  for gid in model.gids:
    r = gid%nhost
    if not rr.has_key(r):
      rr.update({r:[]})
    rr[r].append(gid);
  rr = all2all(rr)
  # rr is now the ranks for where to send the synapse information
  # all the gids in rr were 'owned' by the round-robin distribution
  
  mc = model.mconnections
  gc = model.rank_gconnections
  # construct a ggid2connection dict
  ggid2connection = {}
  for r in gc:
    for ci in gc[r]:
      ggid = ci[3]
      if not ggid2connection.has_key(ggid):
        ggid2connection.update({ggid:[]})
      ggid2connection[ggid].append(ci)
  for r in rr:
    gids = rr[r]
    mgci = []
    rr[r] = mgci
    for gid in gids:
      if mc.has_key(gid):
        mgci.append(mc[gid])
      else:
        mgci.append(ggid2connection[gid])
  mgci = all2all(rr)
  # mgci contains all the connection info needed by the balanced distribution
        
  # create mitrals and granules and register and create synapses
  nmc.dc.mk_mitrals(model)
  nmc.register_mitrals(model)
  nmc.build_granules(model)
  nmc.register_granules(model)
  # build_synapses() ... use mgci to build explicitly
  model.mgrss = {}
  for r in mgci:
    for cil in mgci[r]:
      for ci in cil:
        if not model.mgrss.has_key(mgrs.mgrs_gid(ci[0], ci[3], ci[6])):
          rsyn = mgrs.mk_mgrs(*ci[0:7])
          if rsyn:
            model.mgrss.update({rsyn.md_gid : rsyn})
  nmultiple = int(pc.allreduce(mgrs.multiple_cnt(), 1))
  if rank == 0:
    print 'nmultiple = ', nmultiple
  detectors = h.List("ThreshDetect")
  util.elapsed('%d ThreshDetect for reciprocalsynapses constructed'%int(pc.allreduce(detectors.count(),1)))
  if rank == 0: print 'whole_cell_distrib time ', h.startsw() - enter
