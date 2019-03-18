# -*- coding: cp1252 -*-
from copy import copy
import realSoma
from math import sqrt, sin, cos, pi
from misc import *
from params import *

import params

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
class Extreme:
        SOMA = -1; DENDRITE = 0; APICAL = 1; TUFT = 2
        def __init__(self):
                self.sec = None
                self.phi = 0.
                self.theta = 0.
                self.dist = 0.
                self.limit = 0.
                self.extr_type = Extreme.SOMA
                self.can_bifurc = False
                self.bif_dist = 0.
                self.bif_limit = 0.
                self.depth = 0
                self.basePhi = 0.
                self.baseTheta = 0.
                



# find intersection with Ellipsoid
def EllipsoidLineIntersec(u, p, axis):
  return Ellipsoid(params.bulbCenter,axis).intersection(p,u)

# Ellipsoidd with intersection
def EllipsoidIntersec(p, axis):
  e = Ellipsoid(params.bulbCenter, axis)
  p = e.project(p)
  lamb, phi = e.toElliptical(p)[1:]
  return p, pi*0.5-lamb, phi
  



new_version = False

if new_version:

        
  ''' init the soma position '''
  def gen_soma_pos(r, glompos, mgid):
    nmg=params.Nmitral_per_glom
    glomproj, phi_base, theta_base = EllipsoidIntersec(glompos, params.somaAxis[0])
    phi_base = 0.5*pi - phi_base
    phi_phase = 2*pi/nmg*(mgid%nmg)
    
    # generate a new position
    def gen_pos():
      phi = mgid%nmg * 2*pi / nmg + r.uniform(-0.4*pi/nmg, 0.4*pi/nmg) + phi_phase
      theta = r.normal(0.5*pi, 0.05*(pi*0.5)**2)
      rho = r.uniform(0.2, 1) * params.GLOM_DIST
      phi, theta = convert_direction(phi, theta, phi_base, theta_base)
      return Spherical.xyz(rho, phi, theta, glomproj)

    # check if the position is good
    def in_soma_zone(pos):
      soma_inf = Ellipsoid(bulbCenter, somaAxis[0])
      soma_sup = Ellipsoid(bulbCenter, somaAxis[1])
      return soma_inf.normalRadius(pos) > 1. and soma_sup.normalRadius(pos) < 1.

    while True:
      pos = gen_pos()
      if in_soma_zone(pos):
        break

    return pos

  ''' init the soma '''
  def mk_soma(mgid, nrn):
    sec = Section()
    glompos = params.glomRealCoords[int(mgid / params.Nmitral_per_glom)]
    r = params.ranstream(mgid, params.stream_soma)
    somapos = gen_soma_pos(r, glompos, mgid)

    import realSoma
    index = int(r.discunif(0, realSoma.N_SOMA-1))
    sec.points = realSoma.realSoma(index, somapos)
    nrn.soma = sec

    return somapos


  def can_change_type(extr, glompos):
    if extr.extr_type == Extreme.APICAL:
      return distance(extr.sec.points[-1], glompos) < params.GLOM_RADIUS
    return False

  ''' check the barrier '''
  def feasible(p, extr, nrn, glompos):
    if extr.extr_type == Extreme.APICAL:
      return distance(p, glompos) < distance(extr.sec.points[-1], glompos)
    elif extr.extr_type == Extreme.TUFT:
      d = distance(p, glompos)/params.GLOM_RADIUS
      return d <= 1 and d >= 0.75
    return True
    
  # noise
  def noise(r, b,c=0,d=0):
    return rLaplace(r, 0., b)

  # generate a walk, contains the growing rules.
  def genWalk(r, extr, glomPos):


          
          r = r[extr.extr_type]
          rho = WALK_RHO
          phi = extr.phi
          theta = extr.theta
          diam = 1.                        
          pts = extr.sec.points
          
          if extr.extr_type == Extreme.APICAL:
                  diam = APIC_DIAM
                  dphi = noise(r, params.NS_PHI_B)
                  dtheta = noise(r, params.NS_THETA_B)
                  phi += dphi
                  theta += dtheta
                  absphi, abstheta = convert_direction(phi, theta, extr.basePhi, extr.baseTheta)
                  p = Spherical.xyz(rho, absphi, abstheta, extr.sec.points[-1])
                  p += [ diam ]
                  return p, phi, theta
          elif extr.extr_type == Extreme.TUFT:
                  #def correctTuft(norm, p):
                  #        _rho, phi, theta = Spherical.to(getP(norm, versor(glomPos, p), p), extr.sec.points[-1])
                  #        return phi, theta
                  #
                  #
                  #_p = Spherical.xyz(rho, phi, theta, pts[-1])
                  #dglom = distance(_p, glomPos)
                  #if dglom > 0.9 * GLOM_RADIUS:
                  #        phi, theta = correctTuft(dglom - .9 * GLOM_RADIUS, _p)
                  #elif dglom < 0.6 * GLOM_RADIUS:
                  #        phi, theta = correctTuft(dglom - .6 * GLOM_RADIUS, _p)
                  ##
                  diam = TUFT_DIAM
                  dphi = noise(r, params.NS_PHI_B)
                  dtheta = noise(r, params.NS_THETA_B)
                  phi += dphi
                  theta += dtheta
                  absphi, abstheta = convert_direction(phi, theta, extr.basePhi, extr.baseTheta)
                  p = Spherical.xyz(rho, absphi, abstheta, extr.sec.points[-1])
                  p += [ diam ]
          elif extr.extr_type == Extreme.DENDRITE:
                  ## diam
                  diam = gen_dend_diam(extr.dist)
                  pts = extr.sec.points

                  # first curvature
                  #if len(pts) < 5 and extr.depth == 0:
                    #      theta += (pi / 3. - pi / 15) / 4
                  #        _phi, _theta = ConvertDirection(phi, theta, extr.basePhi, extr.baseTheta, True)
                   #       p = Spherical.xyz(rho, _phi, _theta, pts[-1])

                  
                  # simulate Ellipsoid pression
                  def EllipsoidPression(p, axis, up, k):
                          e = Ellipsoid(bulbCenter, axis)
                          h = e.normalRadius(p)

                          # check if the pression if from up or down to the surface
                          if up:
                                  h_check = h > 1.
                                  h_diff = h
                          else:
                                  h_check = h < 1.
                                  h_diff = 1. - h
                                  
                          q = None
                          if h_check:
                                  q, _lamb, _phi = EllipsoidIntersec(p, axis)
                          else:
                                  _h, _lamb, _phi = e.toElliptical(p)
                                  F = h_diff * k
                                  q = [ -F * sin(_lamb) * cos(_phi) + p[0], -F * cos(_lamb) * cos(_phi) + p[1], -F * sin(_phi) + p[2] ]
                                  
                          vNew = versor(q, pts[-1])
                          _p = getP(rho, vNew, pts[-1])
                          _rho, _phi, _theta = Spherical.to(_p, pts[-1])
                          return _p, _phi, _theta

                  # check
                  _phi, _theta = convert_direction(phi, theta, extr.basePhi, extr.baseTheta)
                  p = Spherical.xyz(rho, _phi, _theta, pts[-1])
                  if not inSomaZone(pts[-1]):
                          e = Ellipsoid(bulbCenter, bulbAxis)
                          if inSomaZone(p):
                                  p, _phi, _theta = EllipsoidPression(p, somaAxis[1], True, 0.)
                                  phi, theta = convert_direction(_phi, _theta, extr.basePhi, extr.baseTheta, True)
                          elif e.normalRadius(pts[-1]) < e.normalRadius(p):
                                  #p, _phi, _theta = EllipsoidPression(p, bulbAxis, True, .6)
                                  p, _phi, _theta = EllipsoidPression(p, bulbAxis, True, GROW_RESISTANCE)
                                  phi, theta = convert_direction(_phi, _theta, extr.basePhi, extr.baseTheta, True)

                  _phi += noise(r, NS_PHI_B, NS_PHI_MIN, NS_PHI_MAX)
                  _theta += noise(r, NS_THETA_B, NS_THETA_MIN, NS_THETA_MAX)
                  p = Spherical.xyz(rho, _phi, _theta, pts[-1])                
                  p.append(diam)
                  return p, phi, theta



                  
          # add the noise        
          phi += noise(r, NS_PHI_B, NS_PHI_MIN, NS_PHI_MAX)
          theta += noise(r, NS_THETA_B, NS_THETA_MIN, NS_THETA_MAX)
          p = Spherical.xyz(rho, phi, theta, extr.sec.points[-1])
          p.append(diam)                               
          return p, phi, theta

  def canDelete(extr, flag_fail):
    if extr.extr_type != Extreme.APICAL:
      return (((extr.dist >= extr.limit or abs(extr.dist - extr.limit) < WALK_RHO) and not extr.can_bifurc) or flag_fail)
    
    return False

  def canBifurcate(extr):
    if extr.extr_type == Extreme.DENDRITE:
      return extr.can_bifurc and (extr.bif_dist >= extr.bif_limit or abs(extr.bif_dist - extr.bif_limit) < WALK_RHO)
    
    return False


                                                                                       
  def updateExtr(extr, p, phi, theta):
          d = distance(p, extr.sec.points[-1])
          extr.dist += d
          extr.bif_dist += d
          extr.sec.points += [ p ]                
          extr.phi = phi
          extr.theta = theta
          
  def mkDendrite(r, depth, parent):
          sec = Section()
          sec.parent = parent
          if parent: parent.sons.append(sec)
          ##
          extr = Extreme()
          extr.extr_type = Extreme.DENDRITE
          extr.sec = sec
          
          # bifurcation decision
          extr.can_bifurc = depth < len(BIFURC_PROB) and r.uniform(0, 1) <= BIFURC_PROB[depth] and extr.dist < (DEND_LEN_MAX - DEND_LEN_MIN)
          extr.depth = depth
         #extr.can_bifurc = False
          if extr.can_bifurc:
                  if depth <= 0:
                          min_len = BIFURC_LEN_MIN_1
                          max_len = BIFURC_LEN_MAX_1
                          mu = BIFURC_LEN_MU_1
                  else:
                          min_len = BIFURC_LEN_MIN_2
                          max_len = BIFURC_LEN_MAX_2
                          mu = BIFURC_LEN_MU_2
                  bl = r.negexp(mu)
                  while bl < min_len or bl > max_len: bl = r.repick()
                  extr.bif_limit = bl
          else:
                  minLen = DEND_LEN_MIN
                  if extr.dist > 0: minLen = extr.dist + DEND_LEN_MIN
                  
                  dendLen = r.normal(DEND_LEN_MU, DEND_LEN_VAR)
                  while dendLen < minLen or dendLen > DEND_LEN_MAX: dendLen = r.repick()
                  extr.limit = dendLen
                  
          return sec, extr




  def mkApical(mid, somaPos, glomPos, nrn, extrLs):
          sec = Section()
          sec.points = [ copy(somaPos) ]
          sec.points[-1].append(APIC_DIAM)
          sec.parent = nrn.soma
          nrn.soma.sons.append(sec)
          nrn.apic = sec
          extr = Extreme()
          extr.sec = sec
          #rho, phi, theta = Spherical.to(glomPos, somaPos)
          extr.phi = 0; extr.theta = 0
          extr.extr_type = Extreme.APICAL
          extr.limit = APIC_LEN_MAX
          extrLs.append(extr)



  def mkTuft(r, extrLs, nrn, pos, glomPos,basePhi,baseTheta):
          # orientation respect glomerulus
          _rho, gl_phi, gl_theta = Spherical.to(glomPos, pos)
          
          # make tuft        
          TUFTS = int(r.discunif(N_MIN_TUFT, N_MAX_TUFT))
          for i in range(TUFTS):
                  sec = Section()
                  sec.index = i
                  sec.parent = nrn.apic
                  nrn.apic.sons.append(sec)
                  nrn.tuft.append(sec)
                  sec.points = [ copy(nrn.apic.points[-1]) ]
                  sec.points[-1][-1] = TUFT_DIAM
                  extr = Extreme()
                  extr.sec = sec
                  extr.phi, extr.theta = convert_direction(i * 2 * pi / TUFTS * (1 + 0.75 * r.uniform(-0.5 / TUFTS, 0.5 / TUFTS)), pi / 4, gl_phi, gl_theta)
                  extr.limit = r.uniform(TUFT_MIN_LEN, TUFT_MAX_LEN)
                  extr.extr_type = Extreme.TUFT
                  extrLs.append(extr)
                  extrLs[-1].basePhi=basePhi;extrLs[-1].baseTheta=baseTheta
          


  # change extremity, especially at the end of apical                                                             
  def change_type(r, extrLs, i, nrn, glomPos):
                                   
          mkTuft(r[Extreme.TUFT], extrLs, nrn, extrLs[i].sec.points[-1], glomPos,extrLs[i].basePhi,extrLs[i].baseTheta)
          


  # make a bifurcation
  def Bifurcate(r, extrLs, i, nrn):
          r = r[Extreme.DENDRITE]
          def get_phi():
                  phi = r.negexp(BIFURC_PHI_MU)
                  while phi < BIFURC_PHI_MIN or phi > BIFURC_PHI_MAX: phi = r.repick()
                  return phi

          
          sec1, extr1 = mkDendrite(r, extrLs[i].depth + 1, extrLs[i].sec)
          sec1.points = [ copy(extrLs[i].sec.points[-1]) ]
          sec1.index = len(nrn.dend)
          nrn.dend.append(sec1)
          
          extr1.phi = extrLs[i].phi + get_phi(); extr1.theta = extrLs[i].theta
          extr1.dist = extrLs[i].dist
          extr1.basePhi = extrLs[i].basePhi; extr1.baseTheta = extrLs[i].baseTheta
          extrLs.append(extr1)

          sec2, extr2 = mkDendrite(r, extrLs[i].depth + 1, extrLs[i].sec)
          sec2.points = [ copy(extrLs[i].sec.points[-1]) ]
          sec2.index = len(nrn.dend)
          nrn.dend.append(sec2)
          extr2.basePhi = extrLs[i].basePhi; extr2.baseTheta = extrLs[i].baseTheta
          extr2.phi = extrLs[i].phi - get_phi(); extr2.theta = extrLs[i].theta
          extr2.dist = extrLs[i].dist
          extrLs.append(extr2)
          

  def inSomaZone(p):
          somaInf = Ellipsoid(bulbCenter, somaAxis[0]); somaSup = Ellipsoid(bulbCenter, somaAxis[1])
          return somaInf.normalRadius(p) >= 1. and somaSup.normalRadius(p) < 1.


  # initialization of algorithm
  def initGrow(mid, noDend):
          
          r = [ ranstream(mid, stream_dend), ranstream(mid, stream_apic), ranstream(mid, stream_tuft) ]
          extrLs = [ ]
          nrn = Neuron()
          somaPos, glomPos = mkSoma(mid, nrn)
          
          mkApical(r[Extreme.APICAL], somaPos, glomPos, nrn, extrLs)

          # make dendrites.
          if not noDend:
                  # get basic inclination
                  somaLvl = Ellipsoid(bulbCenter, somaAxis[0])
                  h, phi_base, theta_base = somaLvl.toElliptical(somaPos)
                  theta_base = pi / 2 - theta_base; phi_base = pi / 2 - phi_base
                  _r = r[Extreme.DENDRITE]
                  
                  ###
                  DENDRITES = int(_r.negexp(N_MEAN_DEND - N_MIN_DEND)) + N_MIN_DEND
                  while DENDRITES > N_MAX_DEND: DENDRITES = int(_r.negexp(N_MEAN_DEND - N_MIN_DEND)) + N_MIN_DEND
                  nmg=params.Nmitral_per_glom
                  dphi=2*pi/DENDRITES
                  for i in range(DENDRITES):
                          phiPhase = i*dphi + 2*pi/nmg*(mid%nmg)
                          
                          sec, extr = mkDendrite(_r, 0, nrn.soma)
                          p = copy(somaPos); p.append(gen_dend_diam(0.))
                          sec.points = [ p ]
                          sec.index = i              
                          nrn.dend.append(sec)
                          
                          ## initial direction
                          phi = phiPhase + _r.uniform(-0.4,0.4)*dphi
                          
                          theta = pi / 3. # - _r.uniform(pi / 15, pi / 12)
                          newPhi, newTheta = convert_direction(phi, theta, phi_base, theta_base)
                          extr.phi = phi; extr.theta = theta
                          extr.basePhi = phi_base; extr.baseTheta = theta_base
                          extrLs.append(extr)
          return r, nrn, extrLs, glomPos


  
  def Grow(mid, noDend):
          r = nrn = extrLs = glomPos = None
          r, nrn, extrLs, glomPos = initGrow(mid, noDend) # initialization        
          
          it = 0
          while it < GROW_MAX_ITERATIONS and len(extrLs) > 0:
                  i = 0
                  while i < len(extrLs):
                          j = 0
                          while j < GROW_MAX_ATTEMPTS:
                                  p, phi, theta = genWalk(r, extrLs[i], glomPos)
                                  
                                  if feasible(p, extrLs[i], nrn, glomPos):
                                          updateExtr(extrLs[i], p, phi, theta)
                                          break
                                  j += 1
                                           
                          if canDelete(extrLs[i], j >= GROW_MAX_ATTEMPTS):
                                  if len(extrLs[i].sec.points) <= 1 and extrLs[i].extr_type == Extreme.DENDRITE:
                                          if extrLs[i].sec.index < (len(nrn.dend) - 1):
                                                  for q in range(extrLs[i].sec.index + 1, len(nrn.dend)):
                                                          nrn.dend[q].index -= 1
                                          del nrn.dend[extrLs[i].sec.index] 
                                  
                                  del extrLs[i]        
                          elif can_change_type(extrLs[i], glomPos):
                                  change_type(r, extrLs, i, nrn, glomPos)
                                  del extrLs[i]
                          elif canBifurcate(extrLs[i]):
                                 Bifurcate(r, extrLs, i, nrn)
                                 del extrLs[i]
                          else:
                                  i += 1
                  it += 1
                  
          return nrn

  def genMitral(mid): return Grow(mid, False)
  def genSomaApicalTuft(mid): return Grow(mid, True)

else:

  def inSomaZone(p):
          somaInf = Ellipsoid(bulbCenter, somaAxis[0]); somaSup = Ellipsoid(bulbCenter, somaAxis[1])
          return somaInf.normalRadius(p) >= 1. and somaSup.normalRadius(p) < 1.
      
  ''' init the soma position '''
  def gen_soma_pos(r, glompos, mgid):
    nmg=params.Nmitral_per_glom
    glomproj1, phi_base, theta_base = EllipsoidIntersec(glompos, params.somaAxis[0])
    glomproj2 = EllipsoidIntersec(glompos, params.somaAxis[1])[0]
    glomproj = centroid([ glomproj1, glomproj2 ])
    
    phi_base = 0.5*pi - phi_base
    
    # generate a new position
    def gen_pos():
      phi = 2*pi/nmg*(mgid%nmg) +2*pi/nmg*r.uniform(-0.4,0.4)
      
      theta = r.normal(0.5*pi,0.05*(0.5*pi)**2)
      rho = r.uniform(0.5, 1) * params.GLOM_DIST
      phi, theta = convert_direction(phi, theta, phi_base, theta_base)
      return Spherical.xyz(rho, phi, theta, glomproj)

    # check if the position is good
    def in_soma_zone(pos):
      soma_inf = Ellipsoid(bulbCenter, somaAxis[0])
      soma_sup = Ellipsoid(bulbCenter, somaAxis[1])
      return soma_inf.normalRadius(pos) > 1. and soma_sup.normalRadius(pos) < 1.

    while True:
      pos = gen_pos()
      if in_soma_zone(pos):
        break

    return pos

  ''' init the soma '''
  def mk_soma(mgid, nrn):
    sec = Section()
    glompos = params.glomRealCoords[int(mgid / params.Nmitral_per_glom)]
    r = params.ranstream(mgid, params.stream_soma)
    somapos = gen_soma_pos(r, glompos, mgid)

    import realSoma
    index = int(r.discunif(0, realSoma.N_SOMA-1))
    sec.points = realSoma.realSoma(index, somapos)
    nrn.soma = sec

    return somapos

  ''' check the barrier '''
  def feasible(p, extr, nrn, glompos):
    if extr.extr_type == Extreme.APICAL:
      return distance(p, glompos) < distance(extr.sec.points[-1], glompos)
    elif extr.extr_type == Extreme.TUFT:
      d = distance(p, glompos)/params.GLOM_RADIUS
      return d <= 1 and d >= 0.75
    return True

  # noise
  def noise(r, b, minval, maxval):
          ns = rLaplace(r, 0., b)
          while ns < minval or ns > maxval: ns = rLaplace(r, 0., b)
          return ns

  # generate a walk, contains the growing rules.
  def genWalk(r, extr, glomPos):


          
          r = r[extr.extr_type]
          rho = WALK_RHO; phi = extr.phi; theta = extr.theta; diam = 1.                        
          pts = extr.sec.points
          if extr.extr_type == Extreme.APICAL:
                  _rho, phi, theta = Spherical.to(glomPos, extr.sec.points[-1])
                  diam = APIC_DIAM                
          elif extr.extr_type == Extreme.TUFT:
                  def correctTuft(norm, p):
                          _rho, phi, theta = Spherical.to(getP(norm, versor(glomPos, p), p), extr.sec.points[-1])
                          return phi, theta
                  #
                  
                  _p = Spherical.xyz(rho, phi, theta, pts[-1])
                  dglom = distance(_p, glomPos)
                  if dglom > 0.9 * GLOM_RADIUS:
                          phi, theta = correctTuft(dglom - .9 * GLOM_RADIUS, _p)
                  elif dglom < 0.6 * GLOM_RADIUS:
                          phi, theta = correctTuft(dglom - .6 * GLOM_RADIUS, _p)
                  ##
                  diam = TUFT_DIAM
          elif extr.extr_type == Extreme.DENDRITE:
                  ## diam
                  diam = gen_dend_diam(extr.dist)
                  pts = extr.sec.points

                  # first curvature
                  #if len(pts) < 5 and extr.depth == 0:
                    #      theta += (pi / 3. - pi / 15) / 4
                  #        _phi, _theta = ConvertDirection(phi, theta, extr.basePhi, extr.baseTheta, True)
                   #       p = Spherical.xyz(rho, _phi, _theta, pts[-1])

                  
                  # simulate Ellipsoid pression
                  def EllipsoidPression(p, axis, up, k):
                          e = Ellipsoid(bulbCenter, axis)
                          h = e.normalRadius(p)

                          # check if the pression if from up or down to the surface
                          if up:
                                  h_check = h > 1.
                                  h_diff = h
                          else:
                                  h_check = h < 1.
                                  h_diff = 1. - h
                                  
                          q = None
                          if h_check:
                                  q, _lamb, _phi = EllipsoidIntersec(p, axis)
                          else:
                                  _h, _lamb, _phi = e.toElliptical(p)
                                  F = h_diff * k
                                  q = [ -F * sin(_lamb) * cos(_phi) + p[0], -F * cos(_lamb) * cos(_phi) + p[1], -F * sin(_phi) + p[2] ]
                                  
                          vNew = versor(q, pts[-1])
                          _p = getP(rho, vNew, pts[-1])
                          _rho, _phi, _theta = Spherical.to(_p, pts[-1])
                          return _p, _phi, _theta

                  # check
                  _phi, _theta = convert_direction(phi, theta, extr.basePhi, extr.baseTheta)
                  p = Spherical.xyz(rho, _phi, _theta, pts[-1])
                  if not inSomaZone(pts[-1]):
                          e = Ellipsoid(bulbCenter, bulbAxis)
                          if inSomaZone(p):
                                  p, _phi, _theta = EllipsoidPression(p, somaAxis[1], True, 0.)
                                  phi, theta = convert_direction(_phi, _theta, extr.basePhi, extr.baseTheta, True)
                          elif e.normalRadius(pts[-1]) < e.normalRadius(p):
                                  #p, _phi, _theta = EllipsoidPression(p, bulbAxis, True, .6)
                                  p, _phi, _theta = EllipsoidPression(p, bulbAxis, True, GROW_RESISTANCE)
                                  phi, theta = convert_direction(_phi, _theta, extr.basePhi, extr.baseTheta, True)

                  _phi += noise(r, NS_PHI_B, NS_PHI_MIN, NS_PHI_MAX)
                  _theta += noise(r, NS_THETA_B, NS_THETA_MIN, NS_THETA_MAX)
                  p = Spherical.xyz(rho, _phi, _theta, pts[-1])                
                  p.append(diam)
                  return p, phi, theta



                  
          # add the noise        
          phi += noise(r, NS_PHI_B, NS_PHI_MIN, NS_PHI_MAX)
          theta += noise(r, NS_THETA_B, NS_THETA_MIN, NS_THETA_MAX)
          p = Spherical.xyz(rho, phi, theta, extr.sec.points[-1])
          p.append(diam)                               
          return p, phi, theta

  def canDelete(extr, flag_fail):
    if extr.extr_type != Extreme.APICAL:
      return (((extr.dist >= extr.limit or abs(extr.dist - extr.limit) < WALK_RHO) and not extr.can_bifurc) or flag_fail)
    
    return False

  def canBifurcate(extr):
    if extr.extr_type == Extreme.DENDRITE:
      return extr.can_bifurc and (extr.bif_dist >= extr.bif_limit or abs(extr.bif_dist - extr.bif_limit) < WALK_RHO)
    
    return False

  def canChangeTypeOfExtr(extr, glomPos): return extr.extr_type == Extreme.APICAL and (extr.dist >= extr.limit or distance(extr.sec.points[-1], glomPos) < GLOM_RADIUS)
                                                                                       
  def updateExtr(extr, p, phi, theta):
          d = distance(p, extr.sec.points[-1])
          extr.dist += d
          extr.bif_dist += d
          extr.sec.points += [ p ]                
          extr.phi = phi
          extr.theta = theta
          
  def mkDendrite(r, depth, parent):
          sec = Section()
          sec.parent = parent
          if parent: parent.sons.append(sec)
          ##
          extr = Extreme()
          extr.extr_type = Extreme.DENDRITE
          extr.sec = sec
          
          # bifurcation decision
          extr.can_bifurc = depth < len(BIFURC_PROB) and r.uniform(0, 1) <= BIFURC_PROB[depth] and extr.dist < (DEND_LEN_MAX - DEND_LEN_MIN)
          extr.depth = depth
         #extr.can_bifurc = False
          if extr.can_bifurc:
                  if depth <= 0:
                          min_len = BIFURC_LEN_MIN_1
                          max_len = BIFURC_LEN_MAX_1
                          mu = BIFURC_LEN_MU_1
                  else:
                          min_len = BIFURC_LEN_MIN_2
                          max_len = BIFURC_LEN_MAX_2
                          mu = BIFURC_LEN_MU_2
                  bl = r.negexp(mu)
                  while bl < min_len or bl > max_len: bl = r.repick()
                  extr.bif_limit = bl
          else:
                  minLen = DEND_LEN_MIN
                  if extr.dist > 0: minLen = extr.dist + DEND_LEN_MIN
                  
                  dendLen = r.normal(DEND_LEN_MU, DEND_LEN_VAR)
                  while dendLen < minLen or dendLen > DEND_LEN_MAX: dendLen = r.repick()
                  extr.limit = dendLen
                  
          return sec, extr



  def mkApical(mid, somaPos, glomPos, nrn, extrLs):
          sec = Section()
          sec.points = [ copy(somaPos) ]
          sec.points[-1].append(APIC_DIAM)
          sec.parent = nrn.soma
          nrn.soma.sons.append(sec)
          nrn.apic = sec
          extr = Extreme()
          extr.sec = sec
          rho, phi, theta = Spherical.to(glomPos, somaPos)
          extr.phi = phi; extr.theta = theta
          extr.extr_type = Extreme.APICAL
          extr.limit = APIC_LEN_MAX
          extrLs.append(extr)



  def mkTuft(r, extrLs, nrn, pos, glomPos):
          # orientation respect glomerulus
          _rho, gl_phi, gl_theta = Spherical.to(glomPos, pos)
          
          # make tuft        
          TUFTS = int(r.discunif(N_MIN_TUFT, N_MAX_TUFT))
          for i in range(TUFTS):
                  sec = Section()
                  sec.index = i
                  sec.parent = nrn.apic
                  nrn.apic.sons.append(sec)
                  nrn.tuft.append(sec)
                  sec.points = [ copy(nrn.apic.points[-1]) ]
                  sec.points[-1][-1] = TUFT_DIAM
                  extr = Extreme()
                  extr.sec = sec
                  extr.phi, extr.theta = convert_direction(i * 2 * pi / TUFTS * (1 + 0.75 * r.uniform(-0.5 / TUFTS, 0.5 / TUFTS)), pi / 4, gl_phi, gl_theta)
                  extr.limit = r.uniform(TUFT_MIN_LEN, TUFT_MAX_LEN)
                  extr.extr_type = Extreme.TUFT
                  extrLs.append(extr)
          
          

  # change extremity, especially at the end of apical                                                             
  def changeTypeExtr(r, extrLs, i, nrn, glomPos):
          
          if distance(extrLs[i].sec.points[-1], glomPos) != GLOM_RADIUS:
                  pos = getP(distance(extrLs[i].sec.points[-1], glomPos) - GLOM_RADIUS,
                             versor(glomPos, extrLs[i].sec.points[-1]),
                             extrLs[i].sec.points[-1])
                  stretchSection(extrLs[i].sec.points, pos)
                  
                  
          mkTuft(r[Extreme.TUFT], extrLs, nrn, pos, glomPos)
          
                  
  # make a bifurcation
  def Bifurcate(r, extrLs, i, nrn):
          r = r[Extreme.DENDRITE]
          def get_phi():
                  phi = r.negexp(BIFURC_PHI_MU)
                  while phi < BIFURC_PHI_MIN or phi > BIFURC_PHI_MAX: phi = r.repick()
                  return phi
          
          sec1, extr1 = mkDendrite(r, extrLs[i].depth + 1, extrLs[i].sec)
          sec1.points = [ copy(extrLs[i].sec.points[-1]) ]
          sec1.index = len(nrn.dend)
          nrn.dend.append(sec1)
          
          extr1.phi = extrLs[i].phi + get_phi(); extr1.theta = extrLs[i].theta
          extr1.dist = extrLs[i].dist
          extr1.basePhi = extrLs[i].basePhi; extr1.baseTheta = extrLs[i].baseTheta
          extrLs.append(extr1)

          sec2, extr2 = mkDendrite(r, extrLs[i].depth + 1, extrLs[i].sec)
          sec2.points = [ copy(extrLs[i].sec.points[-1]) ]
          sec2.index = len(nrn.dend)
          nrn.dend.append(sec2)
          extr2.basePhi = extrLs[i].basePhi; extr2.baseTheta = extrLs[i].baseTheta
          extr2.phi = extrLs[i].phi - get_phi(); extr2.theta = extrLs[i].theta
          extr2.dist = extrLs[i].dist
          extrLs.append(extr2)
          



  # initialization of algorithm
  def initGrow(mid, noDend):
          
          r = [ ranstream(mid, stream_dend), ranstream(mid, stream_apic), ranstream(mid, stream_tuft) ]
          extrLs = [ ]
          nrn = Neuron()

          glomPos = params.glomRealCoords[int(mid/params.Nmitral_per_glom)]
          somaPos = mk_soma(mid, nrn)
          
          mkApical(r[Extreme.APICAL], somaPos, glomPos, nrn, extrLs)

          # make dendrites.
          if not noDend:
                  # get basic inclination
                  somaLvl = Ellipsoid(bulbCenter, somaAxis[0])
                  h, phi_base, theta_base = somaLvl.toElliptical(somaPos)
                  theta_base = pi / 2 - theta_base; phi_base = pi / 2 - phi_base
                  _r = r[Extreme.DENDRITE]
                  
                  ###
                  DENDRITES = int(_r.negexp(N_MEAN_DEND - N_MIN_DEND)) + N_MIN_DEND
                  while DENDRITES > N_MAX_DEND: DENDRITES = int(_r.negexp(N_MEAN_DEND - N_MIN_DEND)) + N_MIN_DEND
                  
                  nmg=params.Nmitral_per_glom
                  dphi=2*pi/DENDRITES
                  phi_phase=_r.uniform(0,2*pi)
                  for i in range(DENDRITES):
                          phiPhase = i*dphi + phi_phase #2*pi/nmg*(mid%nmg)
                          
                          sec, extr = mkDendrite(_r, 0, nrn.soma)
                          p = copy(somaPos); p.append(gen_dend_diam(0.))
                          sec.points = [ p ]
                          sec.index = i              
                          nrn.dend.append(sec)
                          
                          ## initial direction
                          phi = phiPhase
                          
                          theta = pi / 3.
                          newPhi, newTheta = convert_direction(phi, theta, phi_base, theta_base)
                          extr.phi = phi; extr.theta = theta
                          extr.basePhi = phi_base; extr.baseTheta = theta_base
                          extrLs.append(extr)
          return r, nrn, extrLs, glomPos

  def Grow(mid, noDend):
          r = nrn = extrLs = glomPos = None
          r, nrn, extrLs, glomPos = initGrow(mid, noDend) # initialization        
          
          it = 0
          while it < GROW_MAX_ITERATIONS and len(extrLs) > 0:
                  i = 0
                  while i < len(extrLs):
                          j = 0
                          while j < GROW_MAX_ATTEMPTS:
                                  p, phi, theta = genWalk(r, extrLs[i], glomPos)
                                  
                                  if feasible(p, extrLs[i], nrn, glomPos):
                                          updateExtr(extrLs[i], p, phi, theta)
                                          break
                                  j += 1
                                           
                          if canDelete(extrLs[i], j >= GROW_MAX_ATTEMPTS):
                                  if len(extrLs[i].sec.points) <= 1 and extrLs[i].extr_type == Extreme.DENDRITE:
                                          if extrLs[i].sec.index < (len(nrn.dend) - 1):
                                                  for q in range(extrLs[i].sec.index + 1, len(nrn.dend)):
                                                          nrn.dend[q].index -= 1
                                          del nrn.dend[extrLs[i].sec.index] 
                                  
                                  del extrLs[i]        
                          elif canChangeTypeOfExtr(extrLs[i], glomPos):
                                  changeTypeExtr(r, extrLs, i, nrn, glomPos)
                                  del extrLs[i]
                          elif canBifurcate(extrLs[i]):
                                 Bifurcate(r, extrLs, i, nrn)
                                 del extrLs[i]
                          else:
                                  i += 1
                  it += 1
                  
          return nrn

  def genMitral(mid): return Grow(mid, False)
  def genSomaApicalTuft(mid): return Grow(mid, True)

