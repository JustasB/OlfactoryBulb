'''
The concept of gid has two uses. One is administrative so that, given a
gid, one can determine the (pieces that exist) on this machine. The other
is for purposes of spike exchange and is associated with the spike output.
The overall design is to have the necessary ParallelContext.gid_exist(mgid)
associated the mitral cell spike detector and, anytime some piece of the
mitral cell exists on a process, there is a {mgid:cell} entry in mgid2piece.

Because of reciprocal synapses, a pair of spike gids must be reserved
for each synapse and it must be easy to determine from a synapse gid, the
appropriate mitral/granule, and vice versa. Furthermore, splitting consists
of splitting all of a variable number of secondary dendrites from the
soma,primary,tuft,axon region. So we need a single sid for each mitral
which can be the same as the mitral gid.  Python nicely allows us to
have a dict, independent of that supplied by ParallelContext, that says
that one or more pieces exist on this machine.
'''
'''
pc.gid_exists(mgid) returns a positive number only if the soma-priden-tuft-axon
piece exists on this cpu. For other existence questions, use
mpiece_exists, mgid2pieces, and msecden.
Although granules are not currently split, we have their equivalents:
  gpiece_exists, ggid2pieces, and gpriden.
'''
from neuron import h
pc = h.ParallelContext()

try:
  from loadbalutil import lb
except ImportError:
  pass
#dictionary for mgid, mcell (whats left of it)
#model.mgid2piece
from modeldata import getmodel
model = getmodel()

def mpiece_exists(mgid):
  return model.mgid2piece.has_key(mgid)

def mgid2pieces(mgid):
  ''' return cell with existing pieces for mgid; None if does not exist.'''
  if mpiece_exists(mgid):
    return model.mgid2piece[mgid]
  return None

def msecden(mgid, i):
  ''' return secondary dendrite if it exists, otherwise None.'''
  c = mgid2pieces(mgid)
  if c and h.section_exists("secden", i, c):
    return c.secden[i]
  return None

#granules are presently not split but we provide equivalents to the above
def gpiece_exists(ggid):
  return pc.gid_exists(ggid) > 0.0
def ggid2pieces(ggid):
  if gpiece_exists(ggid):
    return pc.gid2cell(ggid)
  return None  
def gpriden(ggid, i):
  c = ggid2pieces(ggid)
  if c and h.section_exists("priden2", i, c):
    return c.priden2[i]
  return None

def secparent(sec):
  sr = h.SectionRef(sec = sec)
  if sr.has_parent():
    return sr.parent
  else:
    return None

def secden_indices_connected_to_soma(cell):
  #list of secden indices connected to soma
  isecden = []
  for i, s in enumerate(cell.secdens):
    if secparent(s) == cell.soma:
      isecden.append(i)
  return isecden

def wholemitral(mgid, cell):
  ''' unsplit mitral cell '''
  model.mgid2piece.update({mgid:cell})

def splitmitral(mgid, cell, piecelist):
  ''' split a mitral cell into secondary dendrites and the soma/priden/axon
      and destroy pieces not on this cpu and connect pieces with multisplit.
      Note that -1 is the piece that includes the soma.
      Also note that secondary dendrites have branches and the piecelist is
      for the indices of secondary dentrites that connect to the soma.
      The {mgid:cell} is added to mgid2piece so the assumption is that
      piecelist is not empty.
  '''
  isecden = secden_indices_connected_to_soma(cell)
  #disconnect all secden and destroy what is not supposed to exist
  for i in isecden:
    s = cell.secden[i]
    h.disconnect(sec = s)
    if i not in piecelist:
      subtree = h.SectionList()
      subtree.wholetree(sec = s)
      for ss in subtree:
        h.delete_section(sec = ss)
  rest = h.SectionList()
  rest.wholetree(sec = cell.soma)
  if -1 not in piecelist:
    for s in rest:
      h.delete_section(sec = s)
  #multisplit connect using mgid
  for i in piecelist:
    if i == -1:
      pc.multisplit(0.5, mgid, sec = cell.soma)
    else:
      pc.multisplit(0.0, mgid, sec=cell.secden[i])
  # add to piece dictionary
  model.mgid2piece.update({mgid:cell})

def subset_complexity(subset):
  cx = 0.
  for sec in subset:
    cx += lb.sec_complexity(sec = sec)
  return cx

def mitral_complexity(cell):
  '''For a full mitral cell return a three-tuple. Total complexity,
     complexity of soma, priden, tuft, axon, and finally a list of
     secondary dendrite subtree complexities.
  '''

  total_cx = 0.

  soma_etc = h.SectionList()
  soma_etc.wholetree(sec = cell.soma)
  soma_etc.remove(cell.secdens)
  soma_etc_cx = subset_complexity(soma_etc)
  total_cx += soma_etc_cx

  isecden = secden_indices_connected_to_soma(cell)
  secden_cx = []
  for i in isecden:
    subtree = h.SectionList()
    subtree.subtree(sec = cell.secden[i])
    secden_cx.append(subset_complexity(subtree))
    total_cx += secden_cx[-1]

  return (total_cx, soma_etc_cx, secden_cx)

def msoma(mgid):
  c = mgid2pieces(mgid)
  if c and h.section_exists('soma', c):
    return c.soma
  return None

def gsoma(ggid):
  c = ggid2pieces(ggid)
  if c and h.section_exists('soma', c):
    return c.soma
  return None

if __name__ == "__main__":
  from mkmitral import mkmitral
  gid = 259
  mcell = mkmitral(gid) # according to mkmitral.py this has tertiary branches
  print "mitral_complexity ", mitral_complexity(mcell)
  print "cell_complexity = ", lb.cell_complexity(mcell)
  pieces = secden_indices_connected_to_soma(mcell)
  pieces.append(-1)
  splitmitral(gid, mcell, pieces)
  h.topology()
  print "mgid2piece ", model.mgid2piece
