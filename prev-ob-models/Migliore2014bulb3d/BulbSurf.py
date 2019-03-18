
import custom_params
if len(custom_params.filename) == 0:
  custom_params.filename='fig7'

import params
from glom import *

from tvtk.api import tvtk

try:
  from enthought.traits.api import HasTraits, Range
  from enthought.traits.ui.api import View, Item
except:
  from traits.api import HasTraits, Range
  from traitsui.api import View, Item

  
class Bulb3d:

   
        
    def _drawRealGloms(self):
        ls = []
        for pos in glomRealCoords:
            src = tvtk.SphereSource(center=tuple(pos), radius=params.GLOM_RADIUS)
            mapper = tvtk.PolyDataMapper(input=src.output)
            actor = tvtk.Actor(mapper=mapper)
            actor.property.color = (1., 0., 0.)
            actor.property.opacity = 0
            self.fig.scene.add_actor(actor)
            ls.append(actor)
        return ls
    
    def __init__(self, fig):
        self.fig = fig
        
        # draw real gloms
        self.real_h = self._drawRealGloms()
                
   
