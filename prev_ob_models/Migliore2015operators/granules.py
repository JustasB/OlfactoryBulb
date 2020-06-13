import params
import misc

ggid2pos = {}
pos2ggid = {}
Ngranule = 0
moves = set([ (0., 0., 0.) ])

# init granules
def initgranules():
  global Ngranule

  ggid2pos.clear()
  pos2ggid.clear()

  eup = misc.Ellipsoid(params.bulbCenter, params.somaAxis[0])
  edw = misc.Ellipsoid(params.bulbCenter, params.granAxisInf)

  for index in range(params.Nx_granule * params.Ny_granule * params.Nz_granule):
    pos = [ 0. ] * 3
    pos[0] = ((index % (params.Nx_granule * params.Ny_granule)) % params.Nx_granule) * params.grid_dim + params.granule_origin[0]
    pos[1] = int((index % (params.Nx_granule * params.Ny_granule)) / params.Nx_granule) * params.grid_dim + params.granule_origin[1]
    pos[2] = int(index / (params.Nx_granule * params.Ny_granule)) * params.grid_dim + params.granule_origin[2]

    
    if eup.normalRadius(pos) < 1. and edw.normalRadius(pos) > 1.:
      pos = tuple(pos)
      ggid = len(ggid2pos) + params.gid_granule_begin
      
      ggid2pos.update({ ggid:pos })
      pos2ggid.update({ pos:ggid })

  Ngranule = len(ggid2pos)



# init moves
def initmoves():
  depth = 1 #int(round(params.granule_field_radius / params.grid_dim))
  #print depth
  d = params.grid_dim

  oldmoves = set([ (0.,0.,0.) ])
  for i in range(depth):
    newmoves = set()
    for m in oldmoves:
      for dx in range(-d, d+1, d):
        for dy in range(-d, d+1, d):
          for dz in range(-d, d+1, d):
            p = (m[0]+dx, m[1]+dy, m[2]+dz)
            newmoves.add(p)
            moves.add(p)
    oldmoves = newmoves

# get granules under position
def get_below(p):

  def get_nears(q):
    npt = []
    for dx, dy, dz in moves:
      npt.append((q[0]+dx, q[1]+dy, q[2]+dz))
    return npt

  def get_p(x, u):
    dim = params.grid_dim
    o = params.granule_origin
    q = ()
    for i in range(3):
      q += (int(round((p[i] + x * u[i] - o[i])/dim)*dim+o[i]),)
    return q

  # scroll the line from points to the center
  u = misc.versor(params.bulbCenter, p)
  L = misc.distance(p, params.bulbCenter)
  dx = L / params.grid_dim

  nnpts = set()
  x = 0.
  while x <= L:
    q = get_p(x, u)
    if q not in nnpts:
      nnpts.update(get_nears(q))
    x += dx
    
  return nnpts

initmoves()
initgranules()

