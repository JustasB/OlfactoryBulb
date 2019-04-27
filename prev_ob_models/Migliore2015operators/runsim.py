from common import *
import params
import util
import parrun
import weightsave
import vrecord as vr
import net_mitral_centric as nmc

def build_part_model(gloms, mitrals, dicfile=''):

  model = getmodel()
  model.clear()

  # gids
  gids = set()
  for glomid in gloms:
    gids.update(range(glomid * params.Nmitral_per_glom, (glomid+1) * params.Nmitral_per_glom))
  gids.update(mitrals)
  
  # distribute
  nmc.build_net_round_robin(model, gids, dicfile)

  build_model()
  
def build_complete_model(dicfile=''):
  build_part_model(range(params.Ngloms), [], dicfile)

def build_model():
  import distribute
  import multisplit_distrib
  multisplit_distrib.multisplit_distrib(distribute.getmodel())

  # set initial weights
  if len(params.initial_weights) > 0:
    weightsave.weight_load(params.initial_weights)


  # print sections
  nc = h.List("NetCon")
  nc = int(pc.allreduce(nc.count(),1))
  if rank == 0: print "NetCon count = ", nc
  nseg = 0
  for sec in h.allsec():
    nseg += sec.nseg
  nseg = int(pc.allreduce(nseg, 1))
  if rank == 0: print "Total # compartments = ", nseg

  pc.spike_record(-1, parrun.spikevec, parrun.idvec)
  util.show_progress(200)

  from odorstim import OdorSequence
  odseq = OdorSequence(params.odor_sequence)

  # record
  for rec in params.sec2rec:
    vr.record(*rec)

  if rank == 0: print 'total setup time ', h.startsw()-startsw

h("proc setdt(){}")
h.dt = 1./64. + 1./128.

def run():
  
  parrun.prun(params.tstop)
  weightsave.weight_file(params.filename + '.weight.dat')
  vr.out()
  util.finish()
