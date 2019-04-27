import sys, pickle

sys.path.extend(["..","../networks","../simulations","../generators","../analysis"])
from networkConstants import *
from stimuliConstants import *
from simset_odor import *
from average_odor_morphs import get_filename
from sim_utils import rebin
from pylab import * # part of matplotlib that depends on numpy but not scipy

from data_utils import * # has axes_labels()

## USAGE: python movie_overlap_plots.py

RUNTIME = REALRUNTIME + SETTLETIME

## time points for the firing rate which is read from a pickled file
firingtsteps = arange(0,RUNTIME+1e-10,FIRINGFILLDT)# include the last RUNTIME point also.
lastrespstart = -int(RESPIRATION/FIRINGFILLDT)

stimseed = 844.0
odornum = 0 # 5 for odorA, 0 for odorB
numgloms = 3
## inh_options = [ (no_singles,no_joints,no_lat,no_PGs,varyRMP), ... ]
inh = (False,False,False,False,False)
resultsdir = '../results/odor_morphs'

## load in the ORN input firing rates
fname = '../generators/firerates/firerates_2sgm_'+str(stimseed)+'.pickle'
f = open(fname,'r')
frateOdorList,fratePulseList,randomPulseList, \
randomPulseStepsList,randomResponseList,kernels \
= pickle.load(f)
f.close()

## load in the mit responses
netseed = stimseed
filename, switch_strs \
    = get_filename(netseed,stimseed,inh,numgloms,\
        None,None,None,resultsdir)
f = open(filename,'r')
(mitral_responses_list,mitral_responses_binned_list) = pickle.load(f)
f.close()

fps = 30
## The simulation is played at 1/10x i.e. 0.5 s (RESPIRATION) over 5.0 s.
## It's actually played at 1/50x, but I'll put each plot for 5 frames.
NUMFRAMES = 5*fps
## We need our bindt to be: (=1/300 s)
BINDT = RESPIRATION/float(NUMFRAMES)
#BINDT = FIRINGFILLDT # 1 ms bins, same as FIRINGFILLDT
NUMBINS = int(RESPIRATION/BINDT) # num bins for 1 resp cycle of 0.5 s
## RASTER WON'T WORK -- for single trial I used mpirank=0,
## no correspondence with multi-trial data used here.
## Also should not use trial-averaging below for multi-trial when RASTER = True.
RASTER = False
## moving average window of 0.1 s or BINDT if RASTER is False/True
if RASTER: BINAVGWIDTH = BINDT
else: BINAVGWIDTH = 0.1
mitral_responses_binned_list = \
    rebin(mitral_responses_list, numbins=NUMBINS,\
        bin_width_time=BINAVGWIDTH, numresps=1)
numavgs = len(mitral_responses_list)
mitral_responses_avg = mean(mitral_responses_binned_list, axis=0)
mitral_responses_std = std(mitral_responses_binned_list, axis=0)

orntimes = arange(FIRINGFILLDT/2.0,RESPIRATION,FIRINGFILLDT)
plottimes = arange(BINDT/2.0,RESPIRATION,BINDT)

def plot_ORN_frates(axlist,lasttime):
    """ 5 for odor A, 0 for odor B.
    Take only the last respiration cycle.
    """
    for glomnum in [1,0,2]:
        ax = axlist[glomnum]
        frate = frateOdorList[glomnum][odornum]
        ornlasttimeidx = int(time/FIRINGFILLDT)+1
        ax.plot(orntimes[:ornlasttimeidx],\
            frate[lastrespstart:lastrespstart+ornlasttimeidx],\
            ['b','r','g'][glomnum],linewidth=linewidth)
        beautify_plot(ax,drawxaxis=False,drawyaxis=False,xticks=[],yticks=[])
        if glomnum==2:
            add_scalebar(ax,matchx=False,matchy=False,hidex=True,hidey=True,\
                sizex=0.2,labelx='0.2 s',sizey=2,labely='2 Hz',\
                bbox_to_anchor=[0.8,-0.2],bbox_transform=ax.transAxes)
        ax.set_xlim(0,plottimes[-1])
        ymin,ymax=ax.get_ylim()
        ax.set_ylim(0,10.5) # 10.5 Hz for all plots, to get common scale-bar

def plot_mitresponses(axlist,lasttimeidx):
    for miti,mitnum in enumerate([2,1,0,4]):
        ax = axlist[miti]
        simresponse = mitral_responses_avg[odornum,mitnum]
        ax.errorbar(x=plottimes[:lasttimeidx],y=simresponse[:lasttimeidx],\
            color=['c','m','r','g'][miti],linewidth=linewidth)
        if mitnum==1:
            add_scalebar(ax,matchx=False,matchy=False,hidex=True,hidey=True,\
                sizex=0.2,labelx='0.2 s',sizey=8,labely='8 Hz',color='w',\
                bbox_to_anchor=[1.0,-0.4],bbox_transform=ax.transAxes)
        beautify_plot(ax,x0min=True,y0min=True,\
            drawxaxis=False,drawyaxis=False,xticks=[],yticks=[])
        ax.set_xlim(0,plottimes[-1])
        if RASTER: ax.set_ylim(0,1/BINAVGWIDTH) # spike raster - 1 spike in BINAVGWIDTH
        else: ax.set_ylim(10,40) # same ylim to set common scale bar

if __name__ == "__main__":
    for timeidx,time in enumerate(plottimes):
        ## MOVIE figures: plots to overlay on simulation movie.
        ## for 1280x720 image, use 300 dpi to convert to inches needed for figsize
        fig = figure(figsize=(4.667,2.4),\
            dpi=300,facecolor='none') # none is transparent
        #inaxlist = [fig.add_subplot(4,4,i) for i in [1,2,4]]
        #outaxlist = [fig.add_subplot(4,4,i) for i in [13,14,15,16]]
        ## axisbg doesn't work, savefig overrides it
        #outaxlist = [fig.add_subplot(4,6,i) for i in [20,15,16,23]]
        ## (left,right,top,bottom) in figure coordinates for the 4 subplots for mits 2,1,0,4 above.
        ## no points of any subplot should "collide" with that of any other subplot, hence the 0.2499.
        coordslist = [(0.225,0.375,0.2499,0.0),(0.35,0.5,0.5,0.25),(0.6,0.75,0.5,0.25),(0.75,0.9,0.2499,0.0)]
        outaxlist = []
        for i in range(4):
            gs = GridSpec(1,1)
            plotcoords = coordslist[i]
            gs.update(left=plotcoords[0],right=plotcoords[1],top=plotcoords[2],bottom=plotcoords[3])
            outaxlist.append( plt.subplot(gs[:,:]) )
        ## Upi suggested that I don't plot the ORN input.
        #plot_ORN_frates(inaxlist,time)
        plot_mitresponses(outaxlist,timeidx)
        fig_clip_off(fig)
        #fig.tight_layout() # doesn't work if using GridSpec
        fig.savefig('../figures/movie/ORN_mitresponses'+str(timeidx).rjust(10,'0')+'.png',\
            dpi=fig.dpi,transparent=True)
        print "Saved plot figure number",timeidx
