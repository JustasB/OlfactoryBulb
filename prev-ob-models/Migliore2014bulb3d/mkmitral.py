from neuron import h
h.load_file('mitral.hoc')
import genmitral

def mkmitral(gid):
  nrn = getmitral(gid)
  
  m = h.Mitral()
  m.createsec(len(nrn.dend), len(nrn.tuft))
  m.subsets()
  m.topol(0) # need to connect secondary dendrites explicitly

  for i, d in enumerate(nrn.dend):
    
    # <<< check my changed if
    if(d.parent == nrn.soma): # <<< changed name
      m.secden[i].connect(m.soma(.5))
    else:
      m.secden[i].connect(m.secden[d.parent.index](1)) # <<< changed name
  
  m.geometry()
  m.segments() # depends on geometry
  m.geometry() # again to get the hillock stylized shape

  fillall(nrn, m)
  
  m.segments() # again to get the proper number of segments for tuft and secden
  m.soma.push()
  m.x = h.x3d(0)
  m.y = h.y3d(0)
  m.z = h.z3d(0)
  h.pop_section()
  m.memb()
  return m

def fillall(n, m):
  fillshape(n.soma, m.soma)
  fillshape(n.apic, m.priden)
  for i,s in enumerate(n.dend):
    fillshape(s, m.secden[i])
  for i,s in enumerate(n.tuft):
    fillshape(s, m.tuftden[i])
  
def fillshape(s1, s2):
    s2.push()
    h.pt3dclear()
    for x in s1.points:
      h.pt3dadd(x[0], x[1], x[2], x[3])
    h.pop_section()

if __name__ == "__main__":
  for mgid in range(635):
    print mgid
    mkmitral(mgid)
  
