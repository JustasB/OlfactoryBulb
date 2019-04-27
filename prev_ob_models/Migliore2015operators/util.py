from common import *
import sys

def serialize():
  ''' Execute body of for loop for each rank in order. Hopefully
      printing will also be segregated and in order
  '''
  for r in range(nhost):
    pc.barrier()
    if r == rank:
      yield r
      sys.stdout.flush()
  pc.barrier()

def group_serialize(ngroup=nhost):
  ''' Execute body of for loop for ngroups of contiguous ranks. 
      The ranks in each group is sequentially executed. Within each barrier,
      one rank in every group is executed in parallel.
  '''
  for r in range(ngroup):
    pc.barrier()
    if r == rank%ngroup:
      yield (rank/ngroup, r)
      sys.stdout.flush()
  pc.barrier()

def finish():
  ''' proper way to quit '''
  if nhost > 0:
    pc.runworker()
    pc.done()
    print 'total elapsed time ', h.startsw()-startsw
    h.quit()
    
elapsedtime=h.startsw()
def elapsed(message):
  ''' Rank 0 prints message and walltime elapsed since
      previous call to this function.
  '''
  global elapsedtime
  if rank == 0:
    print "%s elapsedtime %g"% (message, h.startsw() - elapsedtime)
  elapsedtime = h.startsw()


def progress(pinvl, swlast):
  sw = h.startsw()
  print "t=%g wall interval %g"% (h.t, sw-swlast)
  h.cvode.event(h.t+pinvl, (progress, (pinvl , sw)))

def show_progress(invl):
  global fih
  if rank == 0:
    fih = h.FInitializeHandler(2, (progress, (invl, h.startsw())))

if __name__ == '__main__':
  h.tstop = 1000
  show_progress(200)
  h.run()
