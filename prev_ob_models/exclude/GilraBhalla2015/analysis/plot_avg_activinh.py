import pickle
import os,sys,glob
from pylab import *

sys.path.extend([".."])
from data_utils import *

arraysoff = []
arrayson = []
filelist = glob.glob(sys.argv[1]+'*.pickle')
for filename in filelist:
    f = open(filename,'r')
    (Ainjectarray, dual_firingratearray) = pickle.load(f)
    arraysoff.append(dual_firingratearray[0])
    arrayson.append(dual_firingratearray[1])
    f.close()
arraysoff = array(arraysoff)
arrayson = array(arrayson)

fig = figure()
ax = fig.add_subplot(111)
## am plotting standard deviation of each trial
## and not standard error, since a plot of one trial is being represented.
errorbar(Ainjectarray*1e9,arraysoff.mean(axis=0),\
    yerr=arraysoff.std(axis=0),color='g',marker=',',\
    label="mitral B: 0 nA", linewidth=2.0)
errorbar(Ainjectarray*1e9,arrayson.mean(axis=0),\
    yerr=arrayson.std(axis=0),color='r',marker=',',\
    label="mitral B: 1.5 nA", linewidth=2.0)
#title("Activity dependent inhibition",fontsize=36)
axes_labels(ax, "current in mitral A (nA)",\
    "mitral A soma firing rate (Hz)")
biglegend('lower right')

show()
