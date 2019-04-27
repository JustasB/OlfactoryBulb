from util import *
from loadbalutil import lb

def determine_complexity(model):
  cx = []
  for gid in model.mitrals:
    cx.append((lb.cell_complexity(model.mitrals[gid]), gid))
  for gid in model.granules:
    cx.append((lb.cell_complexity(sec=model.granules[gid].soma),gid))
  return cx

def write_cx_file(name, cx):
  for r in serialize():
    if r == 0:
      f = open(name, 'w')
    else:
      f = open(name, 'a')
    for i in cx:
      s = "%d %g\n"%(i[1], i[0])
      f.write(s)
    f.close()

def test():
  elapsed('Before determine_complexity')
  cx = determine_complexity(getmodel())
  elapsed('determine_complexity')
  cx.sort(key=lambda x:x[0], reverse=True)
  write_cx_file('cx.dat', cx)
  total_cx = 0
  for i in cx:
     total_cx += i[0]
  total_cx = pc.allreduce(total_cx, 1)
  if rank == 0: print 'total_cx=',total_cx
  elapsed('write_cx_file')

if __name__ == '__main__':
  test()
  finish()
