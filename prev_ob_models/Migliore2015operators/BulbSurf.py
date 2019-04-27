
from mayavi.mlab import points3d, mesh
import params
from Glom import *
from grow import *

from tvtk.api import tvtk

try:
  from enthought.traits.api import HasTraits, Range
  from enthought.traits.ui.api import View, Item
except:
  from traits.api import HasTraits, Range
  from traitsui.api import View, Item
  
import numpy
pi = numpy.pi
cos = numpy.cos
sin = numpy.sin

class Bulb3d(HasTraits):

    falseGl = Range(0., 1., 1., desc='False Glomerulus Opacity')
    realGl = Range(0., 1., 1., desc='Real Glomerulus Opacity')
    surf = Range(0., 1., 1., desc='Surface Opacity')
    view = View(Item(name='falseGl'), Item(name='realGl'), Item(name='surf'))
    
    def _surf_changed(self):
      self.bulb_ellipse.actor.property.opacity = self.surf
      
    def _falseGl_changed(self):
      self.false_h.actor.property.opacity = self.falseGl
      
    def _realGl_changed(self):
        flag = self.fig.scene.disable_render
        
        self.fig.scene.disable_render = True
        
        for x in self.real_h:          
          x.property.opacity = self.realGl
            
        if not flag:
          self.fig.scene.disable_render = False
    
    def _cutSurface(self, x, y, z):
        for i in range(len(y)):
            for j in range(len(y[i])):
                if y[i][j] > cuttingY:
                    y[i][j] = cuttingY


    def __mesh_pts(self, depth):
        center = params.bulbCenter
        dphi = 2 * pi / 250
        dtheta = pi / 250
        [phi, theta] = numpy.mgrid[0:2 * pi + dphi:dphi, 0:pi + dtheta:dtheta]
        x = (params.bulbAxis[0] / 2 - depth) * cos(phi) * sin(theta) + center[0]
        y = (params.bulbAxis[1] / 2 - depth) * sin(phi) * sin(theta) + center[1]
        z = []
        for i in range(len(theta)):
            zaux = []
            for j in range(len(theta[i])): zaux.append((bulbHalfAxisZ(theta[i][j]) - depth) * cos(theta[i][j]) + center[2])
            z += [ zaux ]
        self._cutSurface(x, y, z)
        return x, y, z

                
    def __drawEllipse(self, depth, col, op):
        x, y, z = self.__mesh_pts(depth)
        return mesh(x, y, z, color=col, opacity=op)


    def _drawBulb(self):
        self.bulb_ellipse = self.__drawEllipse(0., (0., 1., 0.), self.surf)
        
    def _drawSoma(self):
        self.soma_ellipse = self.__drawEllipse(900, (0., 0., 0.), 1.)
   
        
    def _drawRealGloms(self):
        ls = []
        for pos in glomRealCoords:
            src = tvtk.SphereSource(center=tuple(pos), radius=params.GLOM_RADIUS)
            mapper = tvtk.PolyDataMapper(input=src.output)
            actor = tvtk.Actor(mapper=mapper)
            actor.property.color = (1., 0., 0.)
            actor.property.opacity = self.realGl
            self.fig.scene.add_actor(actor)
            ls.append(actor)
        return ls

    def _drawFalseGloms(self):
        x = []; y = []; z = []
        for pos in glomFalseCoords:
            x.append(pos[0]); y.append(pos[1]); z.append(pos[2])
        return points3d(x, y, z, color=(0., 1., 0.), scale_factor=params.GLOM_RADIUS * 2., opacity=self.falseGl)
    
    def __init__(self, fig):
        self.edit_traits()
        self.fig = fig
        
        # draw real gloms
        self.real_h = self._drawRealGloms()
        
        # draw false gloms
        self.false_h = self._drawFalseGloms()
        
        # draw surface
        self._drawSoma()
        self._drawBulb()

        self.falseGl = 0.
        self.surf = .1
        
        self.sel_descriptor = None
        self.__oldsel = None

      
    def pick_callback(self, picker):
      if self.sel_descriptor:
        if self.__oldsel:
          iglom, col = self.__oldsel
          self.real_h[iglom].property.color = col
        try:
          iglom = self.real_h.index(picker.actor)
          self.sel_descriptor.sel_description = 'glomerulus %d' % iglom
          self.__oldsel = (iglom, self.real_h[iglom].property.color) # store
          self.real_h[iglom].property.color = (1., .5, 1.)
        except:
          self.__oldsel = None
          pass
        
if __name__ == '__main__':
  from mayavi.mlab import figure
  fig = figure(bgcolor=(0.,0.,0.))
  loadGloms()
  b = Bulb3d(fig)
