from util import *

import fileinput
from common import getmodel

def weight_load(filename):
  model = getmodel()
  for l in fileinput.input(filename):
    tk = l.split()
    gid = int(tk[0])
    s = int(tk[1])
    
    # inhib check    
    if gid % 2 != 0:
      gid += 1
      inhib = True
    else:
      inhib = False

    # has key
    if model.mgrss.has_key(gid):
      rsyn = model.mgrss[gid]
      
      if inhib and rsyn.gd2fi:
        rsyn.gd2fi.weight[1] = s
      elif not inhib and rsyn.md2ampanmda:
        rsyn.md2ampanmda.weight[1] = s

  fileinput.close()
        
    
#mg_dict_filename = 'mg_dict.txt'
def weight_file(prefix):
  wtime = h.startsw()
  mingroupsize = max(nhost/64, 1)
  ng = nhost/mingroupsize
  for r in group_serialize(ng):
    name = prefix + '.' + str(r[0])
    
    if r[1]:
      f = open(name, 'a')
      #fdic = open(mg_dict_filename + '.' + str(r[0]), 'a')
    else:
      f = open(name, 'w')
      #fdic = open(mg_dict_filename + '.' + str(r[0]), 'w')
      
    vs = getmodel().mgrss.values
    for rs in vs():
      s = rs.wstr()
      f.write(s)
      #sdic = rs.mg_dic_str()
      #fdic.write(sdic)
    f.close()
    #fdic.close()
    
  if rank == 0 : print "weight_files %s.[0:%d] write time %g s" % (prefix, ng, h.startsw()-wtime)
