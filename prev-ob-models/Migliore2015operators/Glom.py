#from params import bulbCenter, glomAxis, GLOM_RADIUS
from misc import *

#glomEll = Ellipsoid(bulbCenter, glomAxis)

N_TOTAL_GLOMS = 1800.
REAL_GLOMS_XY_FNAME = 'glomxy.txt'
REALGLOMS_FNAME = 'realgloms.txt'
FALSEGLOMS_FNAME = 'falsegloms.txt'
glomRealCoords = []
glomFalseCoords = []
cuttingY = 2200

def loadGloms():
    def _loadGloms(fname, vect):
        # load gloms positions
        f = open(fname)
        line = f.readline()
        while line:
            token = line.split()
            # every line has N glom X Y of i-glomerolous
            x = float(token[0]); y = float(token[1]); z = float(token[2])
            vect.append([ x, y, z ])
            line = f.readline()
        ###
    _loadGloms(REALGLOMS_FNAME, glomRealCoords)
    _loadGloms(FALSEGLOMS_FNAME, glomFalseCoords)
        
def bulbHalfAxisZ(theta):
    import params
    if theta > pi: theta %= pi
    if theta >= 0. and theta <= pi / 2:
        return params.bulbAxis[2] / 2
    y0 = params.bulbAxis[2] / 2
    x0 = pi / 2
    b = 0.52
    a = 300. / (abs(pi / 2) ** b)
    return y0 + a * abs(theta - x0) ** b


    

                    
                    
            
        
    
