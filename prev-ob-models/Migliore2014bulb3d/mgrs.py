'''
mitral-granule reciprocal synapse

patterned after mgrs.hoc of the bulb3test model but allow any number of
secondary dendrite processes (indexed by mitral.secden[i]).  Ie a
connection is defined (in python) by the 6 tuple (mitral_gid,
secden_index, x, granule_gid, priden_index, x). Connection algorithms
allow more than one mgrs with the same mitral and granule.
Therefore, when the function map from (mgid, ggid) to synapse_gid is used,
it may be necessary to do futher disambiguation.
'''

from common import *
import split
nmitral = params.Nmitral
ngranule = granules.ngranule
gid_mgrs_begin = params.gid_granule_begin + ngranule

# it's used to generate excitatory sinapses pair
# and  inhibitory odd
if gid_mgrs_begin % 2 == 0:
  gid_mgrs_begin += 1

''' 20 slot pairs allowed for multiple MGRS with same mgid,ggid. '''
slot2 = 2 * 5

def mgrs_gid(gid_source, gid_target, slot=0):
  ''' Global index for the ThreshDetect object of the reciprocal synapse. '''
  # note: MGRS below uses the explicit assumption that gd_gid = md_gid - 1
  if (gid_source < nmitral): #detector on mitral
    i = (gid_target*nmitral + gid_source + 1)*slot2 + 2*slot + 1 + gid_mgrs_begin
  else: #detector on granule
    i = (gid_source*nmitral + gid_target + 1)*slot2 + 2*slot + gid_mgrs_begin
  return i

def gid2mg(syngid):
  ''' return (mgid, ggid, source_is_mitral, slot) '''
  sgid = syngid
  sgid -= gid_mgrs_begin
  m2g = sgid%2  # 1 if source on mitral
  i = sgid - m2g
  slot = (i%slot2)/2
  i /= slot2
  i -= 1
  mgid = i%nmitral
  ggid = i/nmitral # - nmitral
  if m2g == 1:
    if syngid != mgrs_gid(mgid, ggid, slot):
      print syngid, mgrs_gid(mgid, ggid, slot), mgid, ggid, slot, m2g
    assert(syngid == mgrs_gid(mgid, ggid, slot))
  else:
    assert(syngid == mgrs_gid(ggid, mgid, slot))
  return (mgid, ggid, m2g==1, slot)

class MGRS:
  '''From a mitral and granule synapse location, and consistent with what
     exists on this process, construct the 5 parts of the reciprocal synapse.
     If the granule location exists, then a spine, ThreshDetect, and
     AmpaNmda synapse will be created. If the mitral location exists, then
     a ThreshDetect and FastInhib synapse will be created. The appropriate gid
     for the ThreshDetect instances will be registered. And the appropriate
     NetCons will connect to the synapses.
  '''
  '''
     To allow use of the FastInhibSTDP synapse on the Mitral side of the
     MGRS, there is an additional part which is a negative weight netcon
     connecting from the mitral side ThreshDetect to the FastInhibSTDP
     instance which provides the post synaptic spike timing information.
     There are major administrative differences due to the weight vector
     differences between FastInhib and FastInhibSTDP.
  '''

  def __init__(self, mgid, isec, xm, ggid, ipri, xg, slot):
    self.mgid = mgid
    self.ggid = ggid
    self.slot = slot
    self.xm = xm
    self.xg = xg
    self.isec = isec
    self.ipri = ipri
    
    self.msecden = split.msecden(mgid, isec)
    self.gpriden = split.gpriden(ggid, ipri)
    self.md_gid = mgrs_gid(mgid, ggid, slot)
    self.gd_gid = mgrs_gid(ggid, mgid, slot)
    self.md = None #ThreshDetect on mitral
    self.gd = None #ThreshDetect on granule


    self.fi = None #FastInhib on mitral
    self.ampanmda = None #AmpaNmda on granule
    self.gd2fi = None #NetCon to fi
    self.md2ampanmda = None #NetCon to ampanmda

    if pc.gid_exists(self.md_gid) > 0. or pc.gid_exists(self.gd_gid) > 0.:
      print "md_gid=%d and/or gd_gid already registered" % (self.md_gid, self.gd_gid)
      raise RuntimeError

    if self.msecden:
      self.md = h.ThreshDetect(self.msecden(xm))
      self.fi = h.FastInhib(self.msecden(xm))
      self.fi.gmax = params.inh_gmax
      self.fi.tau1 = params.fi_tau1
      self.fi.tau2 = params.fi_tau2
      pc.set_gid2node(self.md_gid, pc.id())
      pc.cell(self.md_gid, h.NetCon(self.md, None), 1)

    if self.gpriden:
      self.spine = h.GranuleSpine()
      self.spine.neck.connect(self.gpriden(xg))
      self.gd = h.ThreshDetect(self.spine.head(0.5))
      self.ampanmda = h.AmpaNmda(self.spine.head(0.5))
      self.ampanmda.gmax = params.exc_gmax
      pc.set_gid2node(self.gd_gid, pc.id())
      pc.cell(self.gd_gid, h.NetCon(self.gd, None), 1)

    # Cannot be done above because output ports must exist prior to using 
    # an output gid as an input port on the same process.
    if self.fi:
      self.gd2fi = pc.gid_connect(self.gd_gid, self.fi)
      self.gd2fi.weight[0] = 1 # normalized
      self.gd2fi.weight[1] = 0
      self.gd2fi.delay = 1
    if self.ampanmda:
      self.md2ampanmda = pc.gid_connect(self.md_gid, self.ampanmda)
      self.md2ampanmda.weight[0] = 1 #normalized
      self.md2ampanmda.weight[1] = 0
      self.md2ampanmda.delay = 1

  def pr(self):
    print "%d %d <-> %d %d"%(self.mgid, self.md_gid, self.gd_gid, self.ggid)
    if self.msecden:
      print self.msecden.name(), self.md.hname(), self.fi.hname(), self.gd2fi.hname(), " ", int(self.gd2fi.srcgid())
    if self.gpriden:
      print self.gpriden.name(), self.gd.hname(), self.ampanmda.hname(), self.md2ampanmda.hname(), " ", int(self.md2ampanmda.srcgid())


  def mg_dic_str(self):
    s = ''
    if self.gd2fi:
      s += '%d %d %d %g %d\n' % (self.gd_gid, self.ggid, self.ipri, self.xg, self.slot)
    if self.md2ampanmda:
      s += '%d %d %d %g %d\n' % (self.md_gid, self.mgid, self.isec, self.xm, self.slot)
    return s

  def wstr(self):
    ''' return string in proper wsfile format '''
    s = ''
    if self.gd2fi:
      s += '%d %g %g\n'%(self.gd_gid, self.sm(), self.wm())
    if self.md2ampanmda:
      s += '%d %g %g\n'%(self.md_gid, self.sg(), self.wg())
    return s

  def sm(self):
    return self.gd2fi.weight[1]
  def sg(self):
    return self.md2ampanmda.weight[1]
  def wm(self):
    return self.gd2fi.weight[2] * self.fi.gmax
  def wg(self):
    return self.md2ampanmda.weight[2] * self.ampanmda.gmax

def mk_mgrs(mgid, isec, xm, ggid, ipri, xg, slot):
  ''' Return MGRS instance if at least on half exists, otherwise None.'''
  if split.msecden(mgid, isec) or split.gpriden(ggid, ipri):
    return MGRS(mgid, isec, xm, ggid, ipri, xg, slot)
  return None
    
def multiple_cnt():
  cnt = 0;
  for mgrs in getmodel().mgrss.values():
    if mgrs.slot > 0:
      if mgrs.gd: cnt += 1
      if mgrs.md: cnt += 1
  return cnt

if __name__ == "__main__":
  import mkmitral, split
  h.load_file("granule.hoc")

  m = mkmitral.mkmitral(1)
  pieces = split.secden_indices_connected_to_soma(m)
  pieces.append(-1)
  split.splitmitral(1, m, pieces)
  pc.set_gid2node(1, pc.id())
  pc.cell(1, h.NetCon(m.soma(.5)._ref_v, None, sec=m.soma))

  g = h.Granule()
  pc.set_gid2node(10000, pc.id())
  pc.cell(10000, h.NetCon(g.soma(.5)._ref_v, None, sec=g.soma))

  mgrs = MGRS(1, 0, .8, 10000, 0, .1)
  mgrs.pr()
  mgrs2 = MGRS(1, 0, .8, 10000, 0, .1)


