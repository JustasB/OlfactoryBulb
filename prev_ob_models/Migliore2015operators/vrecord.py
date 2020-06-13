
from neuron import h
import params
import modeldata as md
from common import *
import split

if rank == 0:
  print 'vrecord active'
  
tvec = None
trajec = []
trajec_descr = []
    
def record(gid, secid, arc, *_prefix):    
  sec = None
  
  if gid < params.Nmitral: # it's a mitral
    if secid >= 0:
      sec = split.msecden(gid, secid)
    else: #if md.getmodel().mitrals.has_key(gid) and md.getmodel().mitrals[gid].soma:
      sec = split.msoma(gid) #md.getmodel().mitrals[gid].soma
      arc = 0.5
  elif gid < granules.Ngranule + params.Nmitral: # it's a granule
    if secid >= 0:
      sec = split.gpriden(gid, secid)
    else: #if md.getmodel().granules.has_key(gid) and md.getmodel().granules[gid].soma:
      sec = split.gsoma(gid) #md.getmodel().granules[gid].soma
      arc = 0.5

  if sec:
    vec = h.Vector()
    vec.record(sec(arc)._ref_v)
    trajec.append(vec)
    
    if len(_prefix) > 0:
      descr = _prefix[0] + '-'
    else:
      descr = ''
      
    descr += params.filename + params.sym_params_descr()
    descr += ('-%d-%d-%d' % (gid, secid, arc * 1000))
    
    trajec_descr.append(descr)
    global tvec
    if not tvec:
      tvec = h.Vector()
      tvec.record(h._ref_t)


def out():
  if rank == 0:
    print '\tvrecord: output'
  for i, descr in enumerate(trajec_descr):
    f = open('vrec_' + descr + '.dat', 'w')
    for j in range(int(trajec[i].size())):
      f.write(str(tvec.x[j]) + ' ' + str(trajec[i].x[j]) + '\n')
    f.close()
      
  
  
