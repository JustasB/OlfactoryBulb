# granules
import params
import misc

gid2pos = {}
pos2gid = {}
ngranule = 0

#gsyn_cnt = {}
#gsyn_samples = params.Nmitral * params.mg_max_count
#gsyn_mu = 50.
#assert gsyn_mu < gsyn_samples, 'You have to set the synmean value less than synsamples' 
#gsyn_p = gsyn_mu / gsyn_samples


def initgranules():
    global ngranule

    gid2pos.clear()
    pos2gid.clear()
    #gsyn_cnt.clear()
    
    eup = misc.Ellipse(params.bulbCenter, params.somaAxis[0])
    edw = misc.Ellipse(params.bulbCenter, params.granAxisInf)

    for gindex in range(params.Nx_granule * params.Ny_granule * params.Nz_granule):
        pos = [ 0. ] * 3
        pos[0] = ((gindex % (params.Nx_granule * params.Ny_granule)) % params.Nx_granule) * params.grid_dim + params.granule_origin[0]
        pos[1] = int((gindex % (params.Nx_granule * params.Ny_granule)) / params.Nx_granule) * params.grid_dim + params.granule_origin[1]
        pos[2] = int(gindex / (params.Nx_granule * params.Ny_granule)) * params.grid_dim + params.granule_origin[2]
        if eup.normalRadius(pos) < 1. and edw.normalRadius(pos) > 1.:
            pos = tuple(pos)
            gid = len(gid2pos) + params.gid_granule_begin
            gid2pos.update({ gid:pos })
            pos2gid.update({ pos:gid })
            
            ##### don't forget the stream
            #r = params.ranstream(gid, params.stream_dummy_nsyn)
            #n = int(r.binomial(gsyn_samples, gsyn_p))
            #while n == 0:
            #  n = int(r.repick())
            #gsyn_cnt.update({ gid:n })

    ngranule = len(gid2pos)

initgranules()

# this code is used to searching the granule's soma near to the segment

moves = set() # moves table stored to speed up

def initmoves():
  depth = int(round(params.granule_field_radius / params.grid_dim))
  d = params.grid_dim
  
  lastmoves = set([ (0., 0., 0.) ]) 
  for i in range(depth):
    newmoves = set()
    
    for lm in lastmoves:
      for dx in range(-d, d + 1, d):
        for dy in range(-d, d + 1, d):
          for dz in range(-d, d + 1, d):
            newmoves.add((lm[0] + dx, lm[1] + dy, lm[2] + dz))
    moves.update(newmoves)
    lastmoves = newmoves
    
initmoves() 

#out = False
def granule_voxels(p1, p2):
  u = misc.versor(p2, p1)
  def nears(q):
    nn = []
    for dx, dy, dz in moves:
      nn.append((q[0] + dx, q[1] + dy, q[2] + dz))
    return nn


  def pt(x):
    
    p = ( p1[0] + x * u[0], p1[1] + x * u[1], p1[2] + x * u[2])
    p = ( int(round((p[0] - params.granule_origin[0]) / params.grid_dim)) * params.grid_dim + params.granule_origin[0],  \
          int(round((p[1] - params.granule_origin[1]) / params.grid_dim)) * params.grid_dim + params.granule_origin[1],  \
          int(round((p[2] - params.granule_origin[2]) / params.grid_dim)) * params.grid_dim + params.granule_origin[2])
    return p


  L = misc.distance(p1, p2)
  dx = L / params.grid_dim
  
  visited = set()
  x = 0.
  while x <= L:
    visited.add(pt(x))
    x += dx


  nnpts = set() # near points
  for q in visited:
    nnpts.update(nears(q))

  # return ggids
  ggids = set()
  for q in nnpts:
    if pos2gid.has_key(q):
      ggids.add(pos2gid[q]) 

  return list(ggids)

def granule_position_orientation(gid):
    from misc import versor, ellipseLineIntersec as eli
    pos = list(gid2pos[gid])
    u = versor(pos, params.bulbCenter)
    proj = eli(u, pos, params.bulbCenter, params.somaAxis[1])
    return pos, u, proj

