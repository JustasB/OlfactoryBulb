import custom_params
custom_params.filename='g37e1i002'
from common import *
import granules
import misc
import params

from all2all import all2all

gdist_std = .1
doesn=0


''' convert the arc s to xyz pos '''
def xyz(sec, s):
  sec.push()
  n0 = int(h.n3d())
  x = h.Vector(n0) ; y = h.Vector(n0) ; z = h.Vector(n0)
  s0 = h.Vector(n0)
  for i in range(n0):
    x.x[i], y.x[i], z.x[i], s0.x[i] = (h.x3d(i), h.y3d(i), h.z3d(i), h.arc3d(i))
  x.interpolate(s, s0)
  y.interpolate(s, s0)
  z.interpolate(s, s0)
  h.pop_section()
  return (x, y, z)

''' arc of connections position '''
def secden_connection_positions(sec, rng):
  # Assume poisson distribution of positions of synapses on secondary dendrite.
  rng.negexp(1.)
  s = h.Vector()
  x = 0
  L = sec.L
  while True:
    x += params.mean_synapse_interval #* rng.repick()
    if (x > L):
      break
    s.append(x)
  return s      
  
''' connect a dendrite '''
def determine_secden_target(mgid, isec, sec, rng):
  s = secden_connection_positions(sec, rng) # arc positions
  p = xyz(sec, s) # from arc to 3d pos
  s.mul(1 / sec.L)
  
  connection_info = []
  for i, x in enumerate(s):
    connection_info.append((mgid, isec, x, None, 0, None, 0, (p[0][i], p[1][i], p[2][i])))
  return connection_info
  

''' determine all connections info '''
def determine_mitral_target(mgid, mitral):
  connection_info = []
  for i in range(int(mitral.nsecden)):
    if h.section_exists('secden', i, mitral):
      sec = mitral.secden[i]
      rng = params.ranstream(mgid, params.stream_mitraldendconnect + i)
      connection_info += determine_secden_target(mgid, i, sec, rng)
  return connection_info

''' select the granules '''
def connect_to_granule(x, y, z, rng, gconnected=set()):
  gpos = list(granules.get_below([x, y, z]).difference(gconnected))

  eup = misc.Ellipsoid(params.bulbCenter, params.somaAxis[0])
  edw = misc.Ellipsoid(params.bulbCenter, params.granAxisInf)

  gi = int(rng.discunif(0, len(gpos) - 1)) # pick the random granules  
  while len(gpos) > 0:
    if eup.normalRadius(gpos[gi]) < 1 and edw.normalRadius(gpos[gi]) > 1: # inside layer
      break

    del gpos[gi]
    gi = int(rng.discunif(0, len(gpos) - 1)) # pick the random granules 

  # can't connect
  if len(gpos) == 0:
    return -1, -1, ()
  
  gpos = gpos[gi]
  ggid = granules.pos2ggid[gpos]


  gprj = eup.project(gpos) # granule projected pos

  
  u = misc.versor(gprj, params.bulbCenter)
              
  gx = u[0]*(x-gprj[0])+u[1]*(y-gprj[1])+u[2]*(z-gprj[2]) # near pos. on line

  # randomize
  gx += rng.normal(0, gdist_std**2)
  
  # convert to arc
  gx /= params.granule_priden2_len
  
  # force if out of segment
  if gx > 1 or gx < 0: gx = rng.uniform(0, 1)

  return ggid, gx, gpos



''' wait '''
def connections_wait():
  iters = int(pc.allreduce(0, 2)) * params.Nmitral_per_glom 
  for i in range(iters): all2all({})
  
''' connect the target points to the granules '''
def determine_mitral_connections(mgid, mitral):

  ctarget = determine_mitral_target(mgid, mitral) # mitral dendrite's segment
                                             # to connect

  rng = params.ranstream(mgid, params.stream_m2g)
  
  gconnected = set() # granules within the same glom
  ci = []
  if params.m2g_multiple: # multiple conn. within glom
    
    for i in range(len(ctarget)):
      x, y, z = ctarget[i][-1]
      ggid, gx, gpos = connect_to_granule(x, y, z, rng, gconnected)
      if ggid > 0:
        gconnected.add(gpos)
        ci.append(ctarget[i][:3] + (ggid, 0, gx) + ctarget[i][6:]) # update conn. info
        
  else:
    iters = int(pc.allreduce(len(ctarget), 2)) # max iterations is the maximal number of connections
    nmg = params.Nmitral_per_glom
    rbase = (rank-mgid%nmg)
    
    for i in range(iters):
      for mi in range(nmg):
        
        # no more connections needs
        if i >= len(ctarget):
          all2all({})
        else:
          newg = {}
          # current mitral to connect
          if mi == rank%nmg:
            x, y, z = ctarget[i][-1]
            ggid, gx, gpos = connect_to_granule(x, y, z, rng, gconnected)
            

            if ggid > 0:
              gconnected.add(gpos)
              ci.append(ctarget[i][:3] + (ggid, 0, gx) + ctarget[i][6:])

              # fill the message to sent
              for r in range(nmg):
                r += rbase
                if r == rank:
                  continue
                newg.update({ r:gpos })
                
            # send it!  
            all2all(newg)
            
          else:
            # read new connection
            newg = all2all({})
            if newg.has_key(rank):
              gconnected.add(newg[rank])
  return ci


''' serial version of determine_mitral_connections '''
def determine_glom_connections(glomid):
  global doesn
  import mkmitral
  mtarget = [ None ] * params.Nmitral_per_glom
  rng = [ None ] * params.Nmitral_per_glom
  gconnected = [ None ] * params.Nmitral_per_glom
  
  for i in range(params.Nmitral_per_glom):
    mgid = i + glomid * params.Nmitral_per_glom
    mtarget[i] = determine_mitral_target(mgid, mkmitral.mkmitral(mgid))
    rng[i] = params.ranstream(mgid, params.stream_m2g)

    # connection
    if params.m2g_multiple:
      gconnected[i] = set()
    else:
      if i == 0:
        gconnected[i] = set()
      else:
        gconnected[i] = gconnected[0]
      
    
  ciout = []

  hasmore = True
  j = 0
  while hasmore:
    hasmore = False
    for i in range(len(mtarget)):
      if j < len(mtarget[i]):
        hasmore = True
        ci = mtarget[i][j]
        pos = ci[-1]
        ggid, gx, gpos = connect_to_granule(pos[0], pos[1], pos[2], rng[i], gconnected[i])
        if ggid > 0:
          gconnected[i].add(gpos)
          ciout.append(ci[:3] + (ggid, 0, gx) + ci[6:])
        else:
          doesn+=1
          print 'doesn',doesn
    j += 1
  return ciout

if __name__ == '__main__':
  import mkmitral

  if rank == 0:
    m = mkmitral.mkmitral(0)
    print 'rank:',rank
    ci = determine_mitral_connections(0, m)
  elif rank == 1:
    m = mkmitral.mkmitral(1)
    print 'rank:',rank
    ci = determine_mitral_connections(1, m)
  elif rank == 2:
    m = mkmitral.mkmitral(2)
    print 'rank:',rank
    ci = determine_mitral_connections(2, m)
  elif rank == 3:
    m = mkmitral.mkmitral(3)
    print 'rank:',rank
    ci = determine_mitral_connections(3, m)
  elif rank == 4:
    m = mkmitral.mkmitral(4)
    print 'rank:',rank
    ci = determine_mitral_connections(4, m)
  else:
    connections_wait()
  print 'finished', rank
  #print ci


  
