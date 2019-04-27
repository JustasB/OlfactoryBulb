from common import *
h.load_file("spike2file.hoc")


idvec = h.Vector()
idvec.buffer_size(5000000)
spikevec = h.Vector()
spikevec.buffer_size(5000000)

n_spkout_files = max(nhost/64, 1) # each file contains spikes from 64 ranks
n_spkout_sort = min(n_spkout_files*8, nhost) #each file serializes from 8 ranks
# so each sorting rank gathers spikes from nhost/n_spkout_sort ranks
checkpoint_interval = 50000.

def prun(tstop):
  cvode = h.CVode()
  cvode.cache_efficient(1)
  #pc.spike_compress(0,0,1)
  pc.setup_transfer()
  mindelay = pc.set_maxstep(10)
  if rank == 0: print 'mindelay = %g'%mindelay
  runtime = h.startsw()
  exchtime = pc.wait_time()

  inittime = h.startsw()
  h.stdinit()
  inittime = h.startsw() - inittime
  if rank == 0: print 'init time = %g'%inittime
  
  while h.t < tstop:
    told = h.t
    tnext = h.t + checkpoint_interval
    if tnext > tstop:
      tnext = tstop
    pc.psolve(tnext)
    if h.t == told:
      if rank == 0:
        print "psolve did not advance time from t=%.20g to tnext=%.20g\n"%(h.t, tnext)
      break
  
    # save spikes and dictionary in a binary format to
    # make them more comprimibles
    import binsave
    binsave.save(params.filename, spikevec, idvec)
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    
    h.spike2file(params.filename, spikevec, idvec, n_spkout_sort, n_spkout_files)
  
  runtime = h.startsw() - runtime
  comptime = pc.step_time()
  splittime = pc.vtransfer_time(1)
  gaptime = pc.vtransfer_time()
  exchtime = pc.wait_time() - exchtime
  if rank == 0: print 'runtime = %g'% runtime
  printperf([comptime, exchtime, splittime, gaptime])

def printperf(p):
  avgp = []
  maxp = []
  header = ['comp','spk','split','gap']
  for i in p:
    avgp.append(pc.allreduce(i, 1)/nhost)
    maxp.append(pc.allreduce(i, 2))
  if rank > 0:
    return
  b = avgp[0]/maxp[0]
  print 'Load Balance = %g'% b
  print '\n     ',
  for i in header: print '%12s'%i,
  print '\n avg ',
  for i in avgp: print '%12.2f'%i,
  print '\n max ',
  for i in maxp: print '%12.2f'%i,
  print ''
 
if __name__ == '__main__':
  import common
  import util
  common.nmitral = 1
  common.ncell = 2
  import net_mitral_centric as nmc
  nmc.build_net_roundrobin(getmodel())
  pc.spike_record(-1, spikevec, idvec)
  from odorstim import OdorStim
  from odors import odors
  ods = OdorStim(odors['Apple'])
  ods.setup(nmc.mitrals, 10., 20., 100.)
  prun(200.)
  util.finish()
