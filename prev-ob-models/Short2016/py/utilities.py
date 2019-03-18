# utilities.py
"""Contains file reading and writing functions

vector_list=read_vec("filename.dat")
returns a vector_list that then can be set to a neuron vector
nrn_vec = h.Vector(read_vec("filename.dat"))

write_vec("filename.dat", nrn_vec)
will write the nrn_vec to a file in a format readable by this program
matlab, excel, etc.
"""
from neuron import h, gui

def read_vec(filename):
  vector_list=[]
  for line in open(filename,'r'):
    vector_list.append(eval(line))
  return vector_list

def read_nrn_vec(filename):
  return h.Vector(read_vec(filename))

def write_vec(filename, nrn_vec):
  length=nrn_vec.size()
  fid=open(filename,'w')
  for x in nrn_vec:
    fid.write('%G\n' % x)
  fid.close()
