from common import *
from loadbalutil import lb
import net_mitral_centric as nmc
from balance import load_bal
from destroy_model import destroy_model
from all2all import all2all
import split
import util
import mgrs

def determine_multisplit_complexity(model):
  ''' single phase '''
  cxlist = []
  for gid in model.mitrals:
    cxm = split.mitral_complexity(model.mitrals[gid])
    cxlist.append((cxm[1],(gid,-1)))
    for i,cx in enumerate(cxm[2]):
      cxlist.append((cx,(gid,i)))
  for gid in model.granules:
    cxlist.append((lb.cell_complexity(sec=model.granules[gid].soma),(gid,-1)))
  return cxlist

def multisplit_distrib(model):
  enter = h.startsw()
  cxlist = determine_multisplit_complexity(model)
  #start over
  destroy_model(model)
  # but we still have model.mconnections and model.rank_gconnections
  # from the round-robin distribution perspective.
  # we will need to tell the ranks  where to distribute that information

  cxlist = load_bal(cxlist, nhost) #cxlist is the list of (cx,(gid,piece)) we want on this process
  # the new distribution of mitrals and granules
  model.gids = set([item[1][0] for item in cxlist])
  model.mitral_gids = set([gid for gid in model.gids if gid < params.gid_granule_begin])
  model.granule_gids = model.gids - model.mitral_gids
  # for splitting need gid:[pieceindices]
  gid2pieces = {}
  for item in cxlist:
    gid = item[1][0]
    piece = item[1][1]
    if not gid2pieces.has_key(gid):
      gid2pieces.update({gid:[]})
    gid2pieces[gid].append(piece)

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
  # may wish to revisit so that only synapse info for relevant pieces is
  # scattered.
  
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
        
  # create mitrals and granules, split and register and create synapses
  nmc.dc.mk_mitrals(model) # whole cells 
  nmc.build_granules(model)
  for gid in gid2pieces:
    if gid < params.Nmitral:
      split.splitmitral(gid, model.mitrals[gid], gid2pieces[gid])
  pc.multisplit()
  nmc.register_mitrals(model)
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
  detectors = h.List("AmpaNmda")
  util.elapsed('%d ampanmda for reciprocalsynapses constructed'%int(pc.allreduce(detectors.count(),1)))
  detectors = h.List("FastInhib")
  util.elapsed('%d fi for reciprocalsynapses constructed'%int(pc.allreduce(detectors.count(),1)))
  detectors = h.List("ThreshDetect")
  util.elapsed('%d ThreshDetect for reciprocalsynapses constructed'%int(pc.allreduce(detectors.count(),1)))
  if rank == 0: print 'multisplit_distrib time ', h.startsw() - enter

