'''
  Everyone should be using the unique instance, ie model()
'''

class ModelData():
  def __init__(self):
    self.gids = set([])
    self.mitral_gids = set([])
    self.granule_gids = self.gids - self.mitral_gids
    self.mitrals = {}
    self.mconnections = {}
    self.granules = {}
    self.rank_gconnections = {}
    self.mgrss = {}
    self.mgid2piece = {} # see split.py
    self.dummy_syns = []

  def clear(self):
    self.gids.clear()
    self.mitral_gids.clear()
    self.granule_gids.clear()
    self.mitrals.clear()
    self.mconnections.clear()
    self.granules.clear()
    self.rank_gconnections.clear()
    self.mgrss.clear()
    self.mgid2piece.clear()
    self.dummy_syns = []    

modeldata = ModelData()
def getmodel():
  return modeldata
