# -*- coding: cp1252 -*-
import custom_params
if len(custom_params.filename) == 0: custom_params.filename = 'fig7'
from params import Nmitral



class Section:
  def __init__(self):
    self.index = -1
    self.points = []
    self.sons = []
    self.parent = None

class Neuron:
  def __init__(self): 
    self.dend = []
    self.apic = None
    self.tuft = []
    self.soma = None


def getmitral(mgid):
  from struct import unpack

  def sec_read():
    sec = Section()
    parent_index, section_index, n = unpack('>hhH', fi.read(6))
    for i in range(n):
      sec.points.append(list(unpack('>ffff', fi.read(16))))
    return parent_index, section_index, sec



  fi = open('mitral.dump', 'rb')
  offset = [ None ] * Nmitral
  s = 2*Nmitral
  for i in range(Nmitral):
    l = unpack('>H', fi.read(2))[0]
    offset[i] = s
    s += l

  fi.seek(offset[mgid])

  nrn = Neuron()
  nrn.soma = sec_read()[2]
  nrn.apic = sec_read()[2]
  nrn.apic.parent = nrn.soma
  nrn.soma.sons.append(nrn.apic)

  parent_index, section_index, sec = sec_read()
  while section_index != 0:

    if section_index < 0:
      nrn.tuft.append(sec)
      sec.parent = nrn.apic
      nrn.apic.sons.append(sec)
    else:
      nrn.dend.append(sec)
      sec.index = section_index - 1
      if parent_index == 0:
        sec.parent = nrn.soma
        nrn.soma.sons.append(sec)
      else:
        secpar = nrn.dend[parent_index - 1]
        secpar.sons.append(sec)
        sec.parent = secpar

    parent_index, section_index, sec = sec_read()       

  fi.close()
  return nrn




