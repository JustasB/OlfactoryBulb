from common import *
import params
import util
import parrun
import weightsave
import net_mitral_centric as nmc

  
def build_complete_model(connection_file):
  # distribute
  nmc.build_net_round_robin(getmodel(), connection_file)
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

  util.show_progress(200)

  from odorstim import OdorSequence
  odseq = OdorSequence(params.odor_sequence)

  if rank == 0: print 'total setup time ', h.startsw()-startsw

h("proc setdt(){}")
h.dt = 1./64. + 1./128.

def run():
  
  parrun.prun(params.tstop)
  weightsave.weight_file(params.filename + '.weight.dat')
  util.finish()
