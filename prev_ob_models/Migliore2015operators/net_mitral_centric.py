'''
Construct bulb network with round-robin gid distribution and
mitral-centric synapse distribution.
The mitrals have already been constructed.
'''

from common import *
from util import elapsed
from split import wholemitral, mpiece_exists

t_begin = h.startsw()
import determine_connections as dc
h.load_file("granule.hoc")
import mgrs
elapsed('net_mitral_centric after import mgrs')

import params
import granules

def register_mitrals(model):
  '''register mitrals'''
  for gid in model.mitrals:
    if h.section_exists("initialseg", model.mitrals[gid]):
      s = model.mitrals[gid].initialseg
      pc.set_gid2node(gid, rank)
      pc.cell(gid, h.NetCon(s(1)._ref_v, None, sec=s))
      if not mpiece_exists(gid): # must not be doing multisplit
        wholemitral(gid, model.mitrals[gid])
  elapsed('mitrals registered')


def mkgranule(gid):
  # calculate the dendrite length
#  def granule_priden_length(ggid):
#    from misc import distance, Ellipsoid
#    pos = granules.ggid2pos[ggid]
#    prj = Ellipsoid(params.bulbCenter, params.somaAxis[1]).project(pos)
#    return distance(pos, prj)

  g = h.Granule()
  # set length
#  g.priden2[0].L = params.granule_priden2_len
#  g.priden2[0].nseg = int(g.priden2[0].L / 10.) + 1
#  g.priden.L = granule_priden_length(gid)
#  g.priden.nseg = int(g.priden.L / 15.) + 1
#  g.memb()
  
  return g

def build_granules(model):
  '''build granules'''
  model.granules = {}
  for gid in model.granule_gids:
    g = mkgranule(gid)    
    model.granules.update({gid : g})
  elapsed('%d granules built'%int(pc.allreduce(len(model.granules),1)))

def register_granules(model):
  for gid in model.granules:
    g = model.granules[gid]
    pc.set_gid2node(gid, rank)
    pc.cell(gid, h.NetCon(g.soma(.5)._ref_v, None, sec=g.soma))
  elapsed('granules registered')

def build_synapses(model):
  '''construct reciprocal synapses'''
  model.mgrss = {}
  for r in model.rank_gconnections:
    for ci in model.rank_gconnections[r]:
      rsyn = mgrs.mk_mgrs(*ci[0:7])
      if rsyn:
        model.mgrss.update({rsyn.md_gid : rsyn})
  for mgid in model.mconnections:
    for ci in model.mconnections[mgid]:
      #do not duplicate if already built because granule exists on this process
      if not model.mgrss.has_key(mgrs.mgrs_gid(ci[0], ci[3], ci[6])):
        rsyn = mgrs.mk_mgrs(*ci[0:7])
        if rsyn:
          model.mgrss.update({rsyn.md_gid : rsyn})
  nmultiple = int(pc.allreduce(mgrs.multiple_cnt(), 1))
  if rank == 0:
    print 'nmultiple = ', nmultiple
  detectors = h.List("ThreshDetect")
  elapsed('%d ThreshDetect for reciprocalsynapses constructed'%int(pc.allreduce(detectors.count(),1)))


def read_mconnection_info(model, connection_file):
  #model.mconnections = { mgid:[] for mgid in model.mitral_gids }
  from struct import unpack
  fi = open(connection_file, 'rb')
  rec = fi.read(22)
  while rec:
    md_gid, mgid, isec, xm, ggid, xg = unpack('>LLHfLf', rec)
    if mgid in model.mitral_gids:
      slot = mgrs.gid2mg(md_gid)[3]
      cinfo = (mgid, isec, xm, ggid, 0, xg, slot, (0.,0.,0.))
      if not model.mconnections.has_key(mgid):
        model.mconnections.update({ mgid:[] })
      model.mconnections[mgid].append(cinfo)
    rec = fi.read(22)
  fi.close()
  
def build_net_round_robin(model, mitral_gids=set(range(params.Nmitral)), connection_file=''):
  model.clear()
  
  enter = h.startsw()
  for mgid in mitral_gids:
    
    if rank != (mgid % nhost):
      continue

    model.mitral_gids.add(mgid) # add mitral
  # gen. mitrals
  dc.mk_mitrals(model)
  if len(connection_file) > 0:
    # read connection file
    read_mconnection_info(model, connection_file)
  else:
    dc.mk_mconnection_info(model)  
  build_round_robin(model)  
  if rank == 0: print 'build_subnet_round_robin ', h.startsw()-enter
  

def build_round_robin(model):
  dc.mk_gconnection_info(model)
  model.gids = model.mitral_gids.copy()
  model.gids.update(model.granule_gids)
  register_mitrals(model)
  build_granules(model)
  register_granules(model)
  build_synapses(model)
  elapsed('build_net_round_robin')
  if rank == 0: print "round robin setuptime ", h.startsw() - t_begin

#build_net_round_robin(getmodel())

if __name__ == '__main__':
  from util import serialize
  model = getmodel()
  for r in serialize():
    print "rank %d  %d mitrals  %d granules  %d MGRS nmultiple=%d max_multiple=%d" % (r,len(model.mitrals),len(model.granules), len(mgrss), mgrs.nmultiple, mgrs.max_multiple)
