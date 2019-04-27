filedict = 'c10.dic'
fileweights = 'fig7.weight.dat'

try:
  import custom_params
  custom_params.filename='fig7'
except:
  print 'You can\'t import the params file'



from colors import palette
from misc import Spherical
from math import sqrt, sin, cos, pi



def ConvertDirection(phi, theta, phir, thetar, transf):
  def mult(m, v):
    u = [ 0 ] * len(m)
    for i in range(len(m)):
      for j in range(len(m[i])):
        u[i] += v[j] * m[i][j]
    return u

  v = Spherical.xyz(1., phi, theta, [ 0 ] * 3)
  if transf:
    Rz = [ [ cos(phir), -sin(phir), 0. ],
    [ sin(phir), cos(phir), 0. ],
    [ 0., 0., 1. ] ]
    Ry = [ [ cos(thetar), 0., sin(thetar) ],
    [ 0., 1., 0. ],
    [ -sin(thetar), 0., cos(thetar) ] ]
    u = mult(Ry, v)
    u = mult(Rz, u)
  else:
    Rz = [ [ cos(-phir), -sin(-phir), 0. ],
    [ sin(-phir), cos(-phir), 0. ],
    [ 0., 0., 1. ] ]
    Ry = [ [ cos(-thetar), 0., sin(-thetar) ],
    [ 0., 1., 0. ],
    [ -sin(-thetar), 0., cos(-thetar) ] ]
    u = mult(Rz, v)
    u = mult(Ry, u)

  _rho, _phi, _theta = Spherical.to(u, [ 0 ] * 3)
  return _phi, _theta

weights = {}

def loadweights(fweight):
  global weights
  weights = {}
  with open(fweight) as f:
    line = f.readline()
    while line:
      tokens = line.split()
      gid = int(tokens[0]); w = float(tokens[1])
      weights.update({ gid:w })
      line = f.readline()
  
def gcolors():
  gcol = {}
  for gid in weights.keys():
    if gid % 2 == 1 and gid >= params.gid_granule_begin + granules.ngranule:
      ggid = gd.gid_dict[gid + 1][3]
      w_new = weights[gid]
      if gcol.has_key(ggid):
        w_old = gcol[ggid]
        if w_old < w_new:
          gcol[ggid] = w_new
      else:
        gcol.update({ ggid:w_new })
  return gcol



import bindict as gd
gd.load(filedict)
from mayavi.mlab import figure, text
fig = figure(bgcolor=(0, 0, 0))

from mayavi import mlab
points3d=mlab.points3d
import params
import granules

class GranulesManager:

  def __init__(self):
    self.gran = set()
    self.actor1 = None
    self.gran_color = (0., 157 / 255., 157 / 255.)
    self.mitrals = None
    self.colors = None
    self.projected = False
    self._gran = None

  def setup_colors(self, colors):
    self.colors = colors.copy()
  
  def show(self):          
    if self.actor1:
      self.actor1.remove()
    from granules import granule_position_orientation as gpo
    from params import granule_diam as diam
    x = []
    y = []
    z = []
    s = []
    for gid in self.gran:
      p = gpo(gid)[0]
      s.append(self.colors[gid])
      r = params.ranstream(gid, 0)
      r.uniform(-params.grid_dim * .5, params.grid_dim * .5)          
      x.append(p[0] + r.repick())
      y.append(p[1] + r.repick())
      z.append(p[2] + r.repick())
    self.actor1 = points3d(x, y, z, s, scale_factor=diam, scale_mode='none', vmin=0, vmax=100)
    
  def add(self, ci):
    new = set()
    for _ci in ci:
      new.add(_ci[3])
    self.gran.symmetric_difference_update(new)
    self.gran = set(sorted(self.gran))


import misc

from tvtk.api import tvtk
cone_factor = 2.

def vtkCone(p, q):
    from math import pi

    phi_base, theta_base = misc.Spherical.to(q, p)[1:]

    quads = tvtk.CellArray() #vtk.vtkCellArray()
    points = tvtk.Points()   #vtk.vtkPoints()
    
    for i in range(11):
        # rotate
        phi, theta = ConvertDirection((i % 10) * 2 * pi / 10, pi * .5, phi_base, theta_base, True)
        
        # generate  new points
        _p = tuple(misc.Spherical.xyz(p[3] * .5 * cone_factor, phi, theta, p[0:3]))
        _q = tuple(misc.Spherical.xyz(q[3] * .5 * cone_factor, phi, theta, q[0:3]))
        
        # insert points
        points.append(_p)
        points.append(_q)
        
        if i >= 1:
            # create a face            
            quad = tvtk.Quad()
            n = points.number_of_points - 1
            
            quad.point_ids.set_id(0, n - 3) # p
            quad.point_ids.set_id(1, n - 2) # q
            quad.point_ids.set_id(2, n) # q
            quad.point_ids.set_id(3, n - 1) # p
            
            # insert the new face
            quads.insert_next_cell(quad)

    # create the actor
    polydata = tvtk.PolyData(points=points, polys=quads)
    mapper = tvtk.PolyDataMapper(input=polydata)
    actor = tvtk.Actor(mapper=mapper)
    return actor

class vtkMitral:

  
  def __vtkconvert(self):
    def drawsoma():
      pts = self.mitral.soma.points
      center = misc.centroid(pts)

      # calc. soma radius
      radius = 0.
      for p in pts:
        radius += misc.distance(p, center)
      radius /= len(pts)

      radius *= cone_factor

      # versor
      u = tuple(misc.versor(self.mitral.apic.points[0], self.mitral.apic.points[1]))

      src = tvtk.ConeSource(center=tuple(center[0:3]), radius=radius, height=radius, direction=u, resolution=20)
      mapper = tvtk.PolyDataMapper(input=src.output)
      actor = tvtk.Actor(mapper=mapper)
      fig.scene.add_actor(actor)

      actor.property.color = self.soma_color 
      return actor

    # create a colored segment on the scene
    def mksegment(p1, p2, color):
      actor = vtkCone(p1, p2)      
      actor.property.color = color
      fig.scene.add_actor(actor)
      return actor
          
    def drawsection(pts):
      section = []
      for i in range(1, len(pts)):
        section.append(mksegment(pts[i - 1], pts[i], self.section_color))
      return section
    
    #fig.scene.disable_render = True
    self.soma = drawsoma()
    self.apic = drawsection(self.mitral.apic.points)
    for i in range(len(self.mitral.tuft)):
      self.tuft.append(drawsection(self.mitral.tuft[i].points))


    # gid associated to a segment
    self.dend_gids = []
    for i in range(len(self.mitral.dend)):
      self.dend.append(drawsection(self.mitral.dend[i].points))
      
      aux = []
      for i in range(len(self.dend[-1])):
        aux.append(set())
      self.dend_gids.append(aux)


  def __interp_rules(self, tointerp):

    def move(key, forward):
      newkey = set()
      
      for _key in key:
        sectype, isec, iseg = _key
        if forward:
          iseg += 1
          if iseg < len(self.dend[isec]):
            newkey.add((sectype, isec, iseg))
            
          elif len(self.mitral.dend[isec].sons) > 0:
            for _son in self.mitral.dend[isec].sons:
              newkey.add((sectype, self.mitral.dend.index(_son), 0))
              
        else:
          iseg -= 1
          if iseg >= 0:
            newkey.add((sectype, isec, iseg))
            
          elif self.mitral.dend[isec].parent != self.mitral.soma:
            isec = self.mitral.dend.index(self.mitral.dend[isec].parent)
            iseg = len(self.dend[isec])-1
            newkey.add((sectype, isec, iseg))
      return newkey
    
    def look(keys, forward=True):
      retkeys = set()
      while len(keys) > 0:
        _keys = set()
        
        for k in move(keys, forward):
          if k in tointerp:
            _keys.add(k)
          else:
            retkeys.add(k)
            
        keys = _keys        
      return retkeys
        
    def look_previous(key): return look(key, False)
    
    def look_forward(key): return look(key)

    interprules = {}
    for key in tointerp:
      interprules.update({key:look_forward(set([key])).union(look_previous(set([key])))})
    return interprules
      
  
  def __init__(self, mgid):
    self.mgid = mgid

    from getmitral import getmitral
    self.mitral = getmitral(mgid)


    self.soma = None
    self.apic = None
    self.dend = []
    self.tuft = []

    self.soma_color = (250. / 255, 210. / 255, 51. / 255)
    self.section_color = (1., 1., 1.)

    self.__vtkconvert()
    
    self.conn_info = []

    tointerp = set()
    for isec in range(len(self.dend)):
      for iseg in range(len(self.dend[isec])):
        tointerp.add((2, isec, iseg))
        
    for gid in gd.mgid_dict[mgid]:
      if gid >= params.gid_granule_begin + granules.ngranule:
        if gid % 2 == 0:
          self.conn_info.append(gd.gid_dict[gid])
          isec, x = gd.gid_dict[gid][1:3]
          iseg = int(x*len(self.dend[isec]))
          if iseg >= len(self.dend[isec]):
            iseg = len(self.dend[isec])-1
          tointerp.discard((2, isec, iseg))

    self.__tointerp = self.__interp_rules(tointerp)
    self.__show_weights = False


      
  def __set_segment(self, info):
    for secinfo, ind_color in info.items():

      # interpolate
      sec_type, isec, iseg = secinfo

      if sec_type == -1: # soma
        o = self.soma
      elif sec_type == 0: # tuft
        o = self.tuft[isec][iseg]
      elif sec_type == 1: # apical
        o = self.apic[iseg]
      elif sec_type == 2: # dendrites
        o = self.dend[isec][iseg]

      # set
      o.property.color = palette[ind_color]

    # interpolate
    for key1, linked in self.__tointerp.items():
      ind_color = 0
      for key2 in linked:
        ind_color += info[key2]
      ind_color /= len(linked)
      self.dend[key1[1]][key1[2]].property.color = palette[ind_color]

  def show_weights(self, excit):
    
    if len(weights) == 0:
      return

    self.__show_weights = True
    self.__excit = excit

    w = {}
    for gid in gd.mgid_dict[self.mgid]:
      if gid >= params.gid_granule_begin + granules.ngranule:

        if gid % 2:
          continue

        isec, x = gd.gid_dict[gid][1:3]

        if x >= 1:
          iseg = len(self.dend[isec]) - 1
        else:
          iseg = int(x * len(self.dend[isec]))

          
        if not excit:
          gid -= 1

        try:
          wsym = weights[gid]
          try:
            w[(2, isec, iseg)].append(wsym)
          except KeyError:
            w.update({ (2, isec, iseg):[ wsym ] })
        except KeyError:
          print 'gid %d not found' % gid



    max_steps = 100.

    for k, x in w.items():
      w[k] = int(max(x) / max_steps * (len(palette) - 1))

    self.__set_segment(w)
          

  def clean_weight(self):
    self.__show_weights = False
    
    # color all black
    for dnd in self.dend:
      for seg in dnd:
        seg.property.opacity = 1.

    self.clean()
    
  def __color_section(self, sec, color):
    for s in sec:
      s.property.color = color
      
  def clean(self):
    if self.__show_weights:
      self.show_weights(self.__excit)
    else:
      for sec in self.dend:
        self.__color_section(sec, self.section_color)
    for sec in self.tuft:
      self.__color_section(sec, self.section_color)
    self.__color_section(self.apic, self.section_color)
    self.soma.property.color = self.soma_color


try:
  from enthought.traits.api import HasTraits, Range, String, Button, Int, Bool, Str, Float
  from enthought.traits.ui.api import View, Item
except:
  from traits.api import HasTraits, Range, String, Button, Int, Bool, Str, Float
  from traitsui.api import View, Item
  
class BulbGUI(HasTraits):
  
  w_excit = Button('Weights Excit.')
  w_inhib = Button('Weights Inhib.')
  w_clean = Button('Weights Clean')

  t_stop = Float
  t_win = Float
  max_freqs = Float
  
  view = View(Item(name='w_excit'), Item(name='w_inhib'), Item(name='w_clean'))
  
  def __init__(self, mbp):
    self.edit_traits()
    self.mbp = mbp


  def _w_excit_fired(self):
    self.mbp.show_weights(True)
    fig.scene.render()

  def _w_inhib_fired(self):
    self.mbp.show_weights(False)
    fig.scene.render()

  def _w_clean_fired(self):
    self.mbp.clean_weights()
    fig.scene.render()


class mayaBulbPlot:
  def __init__(self):
    self.sel_descriptor = None
    self.mitrals = {}
    self.granules = GranulesManager()

  def draw_mitral(self, mgid):
    if self.mitrals.has_key(mgid):
      return

    m = vtkMitral(mgid)
    m.sel_color = (0,1,0)
    m.granules = self.granules
    self.mitrals.update({ mgid:m })
    self.granules.add(m.conn_info) # draw granules
    self.granules.show()

  def show_weights(self, excit):
    for m in self.mitrals.values():
      m.show_weights(excit)

  def clean_weights(self):
    for m in self.mitrals.values():
      m.clean_weight()

      
mbp = mayaBulbPlot()
mbp.granules.mitrals = mbp.mitrals

loadweights(fileweights)
mbp.granules.setup_colors(gcolors())

gui = BulbGUI(mbp)
fig.scene.disable_render=True
mbp.draw_mitral(185)
mbp.draw_mitral(202)
mbp.draw_mitral(137)

fig.scene.disable_render=False


import BulbSurf
import odordisp as od
l,val=od.OdorsInput('input-odors.txt')
ii=l.index('Mint')
b=BulbSurf.Bulb3d(fig)
for i in [187/5,137/5,202/5]:
  b.real_h[i].property.opacity=0.5
  b.real_h[i].property.color = palette[val[ii][i]]
  
fig.scene.camera.view_up=[ 0.01551967,  0.99458581,  0.10275311]
fig.scene.camera.position=[ 3521.84237552,   621.77602255,  4977.47240287]
