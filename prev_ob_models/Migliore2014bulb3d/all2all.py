'''
input: a map with destination ranks as keys, the values are pickleable objects.
output: a map with source ranks as keys, each value is an object.
Transfer the map value (pickleable python object) to the proper destination rank,
A value associated with the destination rank will appear on the
destination rank as a value associated with the source rank.
Alternatively, if the input is a list of nhost objects (some can be None),
then so is the return value.
'''

from common import *
import util

ptime = False

def all2all(data, size=0):
  enter = h.startsw()
  r = _all2all(data, size)
  if ptime and rank == 0: print 'all2all elapsed time = %g'% (h.startsw()-enter)
  return  r

def _all2all(data, size=0):
  if nhost == 1:
    if size == -1:
      return (0, 0)
    return data
  if type(data) is list:
    return pc.py_alltoall(data, size)
  elif type(data) is dict:
    d = []
    for i in range(nhost):
      d.append(None)
    for i in data:
      d[i] = data[i]
    d = pc.py_alltoall(d, size)
    if size == -1:
      return d
    z = {}
    for i,x in enumerate(d):
      if x != None:
        z.update({i : x})
    return z
  raise ValueError

if __name__ == '__main__':
  d = []
  for i in range(nhost):
    d.append(i+10)
  sizes = all2all(d, -1)
  d = all2all(d)
  for r in util.serialize():
    print rank, sizes, d
