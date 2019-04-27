from common import *

''' destroy mitrals, granules, MGRS, and clear internal gid maps. '''
def destroy_model(model):
  pc.gid_clear()
  model.mgrss = {}
  model.granules = {}
  model.mitrals = {}
  model.mgid2piece = {}
