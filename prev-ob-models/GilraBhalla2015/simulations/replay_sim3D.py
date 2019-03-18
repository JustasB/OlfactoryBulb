import sys, math
from PyQt4 import QtGui,QtCore
from PyQt4.Qt import Qt
import time#, sched
import pickle
from updatepaintGL_aditya2 import newGLWindow

from oglfunc.objects import *
from oglfunc.group import *

## USAGE: python replay_sim3D.py <data file>
## TYPICAL: python replay_sim3D.py sim_saved.pickle

#s = sched.scheduler(time.time, time.sleep)

def rotate_y(x,y,z,theta):
    return (x*math.cos(theta)-z*math.sin(theta), y, x*math.sin(theta)+z*math.cos(theta))

f = open(sys.argv[1],'r')
q = pickle.load(f)
f.close()

def updateframe(w,framenum):
    if framenum>=2:
        frameoldold = framenum-2
        frameold = framenum-1
    elif framenum>=1:
        frameoldold = framenum-1
        frameold = framenum-1
    else:
        frameoldold = framenum
        frameold = framenum
    dmin = w.minVal
    for i,vizObject in enumerate(w.vizObjects):
        ## decay the spike slowly over next three frames
        dbin = [q[framenum][i],q[frameold][i],q[frameoldold][i]]
        dmaxnow = max(dbin)
        d = (dmaxnow-dmin)/(dbin.index(dmaxnow)*0.3+1.0)
        cell_identity = vizObject.cell_identity
        if cell_identity == 0:
            idx = int(d*w.factor1)
            vizObject.r,vizObject.g,vizObject.b=w.colorMap1[idx]            
        elif cell_identity == 1:
            idx = int(d*w.factor2)
            vizObject.r,vizObject.g,vizObject.b=w.colorMap2[idx]
        elif cell_identity == 2:
            idx = int(d*w.factor3)
            vizObject.r,vizObject.g,vizObject.b=w.colorMap3[idx]
    w.updateGL()
    
def separate_cells(w):
    for i,vizObject in enumerate(w.vizObjects):
        ## moose path is stored in obj.l_coords[-1]
        obj_moosepath = vizObject.l_coords[-1]
        if 'mitral' in obj_moosepath:
            vizObject.cell_identity = 0
        elif 'granule' in obj_moosepath:
            vizObject.cell_identity = 1
        elif 'PG' in obj_moosepath:
            vizObject.cell_identity = 2

app = QtGui.QApplication(sys.argv)
newWin = newGLWindow()
#self.newWin.setWindowState(Qt.WindowMaximized)
#self.newWin.setWindowState(Qt.WindowFullScreen)
w =  newWin.mgl ## instance of updatePaintGL
w.resizeGL(1024,768)
w.setMinimumSize(1024,768)
w.setMaximumSize(1024,768)
w.adjustSize()
newWin.show()
## don't have lights - colors should be the same from any direction
w.lights = False
## turn on visualization, use with self.qgl.updateViz()	
w.viz = 1 # after this all cells drawn will get color updated
## set the color map for visualization
#w.setColorMap(vizMinVal=-100e-3,vizMaxVal=80e-3,cMap='jet')
w.setColorMaps(vizMinVal=-90e-3,vizMaxVal=70e-3)

## last entry in q is configuration of objects in 3D
config = q[-1]
for i in range(len(config)):
    a = locals()[config[i][0]](w,config[i][1],config[i][2])
    a.setCellParentProps(config[i][3],config[i][4],1,0,0)
    w.vizObjects.append(a)

tiltangle = 50.0 # degr
tot_rot = 20.0
## rotate during simulation from -tot_rot/2.0 to +tot_rot/2.0
w.rotate([0,0,1],-tot_rot/2.0) # rotate about axis-vector (z-axis) by xx (10) degrees
## view at tilt angle
w.rotate([1,0,0],-50.0) # rotate about axis-vector (y-axis) by xx (-50) degrees
#w.rotate([0,1,0],-10) # rotate about axis-vector (y-axis) by xx (-10) degrees
x0,y0,z0 = (-1.5,-4,-130)
#x0,y0,z0 = (0,-4.25,-130)
w.translate([x0,y0,z0])
w.updateGL()

## categorize cells upfront as mit/gran/PG based on moose path,
## else too expensive to repeat string check for each frame
separate_cells(w)

framemin = 600
framemax = 900
theta_dt = tot_rot/float(framemax-framemin) # degrees

## leave last entry in q which is configuration, rest are Vm-s at different times
for framenum,vals in enumerate(q[0:-1]):
    if framenum<600 or framenum>900: continue
    ## translate and rotate back to original config
    w.translate([-x0,-y0,-z0])
    w.rotate([1,0,0],50.0)
    ## rotate about y
    w.rotate([0,0,1],theta_dt)
    ## again tilt it and translate for viewing
    w.rotate([1,0,0],-50.0)
    w.translate([x0,y0,z0])
    print "at framenum =",framenum
    updateframe(w,framenum)
    ## save pictures
    ####### grabFrameBuffer() has problems: a frame in 300 or so gets skewed / corrupted!
    pic = w.grabFrameBuffer()
    pic.save('movie/sim_'+str(framenum)+'.png','PNG')

#app.exec_()
sys.exit(0)
