'''
OdorStim supplies an odors[name] stimulus to each mitral tuft dendrite
defined by odors[name].glom_weights. Thus there is a separate NetCon for
each tuft dendrite on this process.
'''
from common import *
from gidfunc import *
import fileinput
import params
from odors import odors

class OdorStim():
  def __init__(self, od, start, dur, rel_conc=1.):
    ''' Specifies the odor for an OdorStim. Note that the OdorStim is
        activated with setup which can only be called after the mitrals
        dict exists (usually from determine_connections.py).
    '''
    # set odor weights
    if type(od) == str:
      self.odor = odors[od]
    else:
      self.odor = od

    self.rel_conc = rel_conc
    self.verbose = True
    self.tstop = start + dur

    
    mitrals = getmodel().mitrals
    self.netcons = {}
    self.rng_act = params.ranstream(0, params.stream_ods_act)
    self.rng_act.uniform(params.ods_freql, params.ods_freqh)
    
    for gid in mitrals:
      m = mitrals[gid]
      
      # in case of multisplit
      if not h.section_exists("tuftden", 0, m):
        continue
      
      iglom = mgid2glom(gid)
      w = self.odor.glom_weights[iglom]
      for i in range(int(m.synls.count())):
        nc = h.NetCon(None, m.synls.o(i))
        # rng for weights
        rw = params.ranstream(gid, params.stream_ods_w + i)
        nc.weight[0] = w * self.rel_conc * rw.uniform(params.ods_wl, params.ods_wh)

        self.netcons.update({(gid, i):(nc, rw)})

    self.fih = h.FInitializeHandler(0, (self.init_ev, (start,)))

  def init_ev(self, start):
    ''' first event at start.
        In principle this can be called by the user during a simulation
        but if h.t < the previous stop, the randomness will be mixed
        between the multiple (start,stop) intervals.
    '''
    if params.sniff_invl == None:
      self.nxt_invl = self.rng_act.repick()
      start += self.nxt_invl
    else:
      self.nxt_invl = params.sniff_invl
    h.cvode.event(start, (self.ev, (start,)))

  #def ev(self, time, interval, stop):

  def ev(self, time):
    ''' time is the standard time with no randomness.
        h.t may be before or after and each syapse will receive its event
        at h.t + individual netcon delay
    '''

    # update weights
    for key in self.netcons:
      iglom = mgid2glom(key[0])
      w = self.odor.glom_weights[iglom]
      nc, rw = self.netcons[key]
      nc.weight[0] = w * self.rel_conc * rw.repick()
      nc.delay = 0.

      # call event queue
      nc.event(h.t)

    # fix next activations
    if params.sniff_invl == None:
        self.nxt_invl = 1000. / self.rng_act.repick()

    if time + self.nxt_invl < self.tstop:
      if rank == 0 and self.verbose: print 'activation of %s at %.3g (ms)\tinterval %.3g' % (self.odor.name, h.t, self.nxt_invl)
      h.cvode.event(time + self.nxt_invl, (self.ev, (time + self.nxt_invl,)))

# create odorseq
def OdorSequence(seq):
  odseq = []
  
  if type(seq) == str:
    for i, line in enumerate(fileinput.input(seq)):
      tk = line.split()
      
      if len(tk) < 4 and rank == 0:
        print 'line %i of %s was ignored' % (i, seq)
        continue
      
      name = tk[0]
      init = float(tk[1])
      dur = float(tk[2])
      conc = float(tk[3])
      
      odseq.append(OdorStim(name, init, dur, conc))
      
  elif type(seq) == list:
    for odinfo in seq:
      odseq.append(OdorStim(*odinfo))
      
  return odseq

if __name__ == '__main__':
  h.load_file("nrngui.hoc")
  import common
  common.nmitral = 2
  common.ncell = 10
  import determine_connections
  import odors
  ods = OdorStim(odors.odors['Apple'])
  ods.setup(determine_connections.mitrals, 10., 20., 100.)
  h.tstop = 150
  h.run()
  print 't=', h.t
