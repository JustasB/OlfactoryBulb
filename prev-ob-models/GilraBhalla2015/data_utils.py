import string
## Setting GTK backend does not work on the gj cluster,
## as I don't have a local pyGTK for my local python2.6.
## So ignore when on gj or its nodes.
import socket
hname = socket.gethostname()
if 'gulabjamun' not in hname and 'node' not in hname:
    import matplotlib
    #matplotlib.use('Agg')
    #matplotlib.use('GTK')
from pylab import *
from matplotlib import collections
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from poisson_utils import *

## choose figure or poster defaults
poster = False

if not poster:
    ####### figure defaults
    label_fontsize = 8 # pt
    plot_linewidth = 0.5 # pt
    linewidth = 1.0#0.5
    axes_linewidth = 0.5
    marker_size = 3.0 # markersize=<...>
    cap_size = 2.0 # for errorbar caps, capsize=<...>
    columnwidth = 85/25.4 # inches
    twocolumnwidth = 174/25.4 # inches
    linfig_height = columnwidth*2.0/3.0
    fig_dpi = 300
else:
    ####### poster defaults
    label_fontsize = 12 # pt
    plot_linewidth = 1.0 # pt
    linewidth = 1.0
    axes_linewidth = 1.0
    marker_size = 3.0
    cap_size = 2.0 # for errorbar caps
    columnwidth = 4 # inches
    linfig_height = columnwidth*2.0/3.0

#######################

########## You need to run:
## From gj:
## ./restart_mpd_static
## The 0th (boss process) will always be node000 as it is the first node in ~/hostfile.
## Hence from node000: cd to the working directory simulations/
## (so that sys.path has accurate relative paths)
## mpiexec -machinefile ~/hostfile -n <numprocs> ~/Python-2.6.4/bin/python2.6 <script_name> [args]
## 0 rank process is for collating all jobs. (rank starts from 0)
## I assume rank 0 process always runs on the machine whose
## X window system has a Display connected and can show the graphs!!!!!!
## The rank 0 stdout is always directed to the terminal from which mpiexec was run.
## I hope X output also works the same way.
## For long simulations save results in a text file
## for replotting later and avoid above ambiguity.
from mpi4py import MPI

mpicomm = MPI.COMM_WORLD
mpisize = mpicomm.Get_size() # Total number of processes
mpirank = mpicomm.Get_rank() # Number of my process
mpiname = MPI.Get_processor_name() # Name of my node
# The 0th process is the boss who collates/receives all data from workers
boss = 0
print 'Process '+str(mpirank)+' on '+mpiname+'.'

def calc_STA( inputlist, spiketrain, dt, STAtime):
    """ spiketrain is a list of output spiketimes,
    inputlist is the input timeseries with dt sample time,
    STAtime is time for which STA must be computed.
    User must ensure that inputlist is at least
    as long as max spiketime in spiketrain.
    endidx = int(spiketime/dt) not round(spiketime/dt)
    as index 0 represents input from time 0 to time dt.
    returns number of relevant spikes and
    sum of Spike Triggered Averages as a list of length STAtime/dt.
    """
    lenSTA = int(STAtime/dt)
    STAsum = zeros(lenSTA)
    numspikes = 0
    for spiketime in spiketrain:
        if spiketime<STAtime: continue
        endidx = int(spiketime/dt)
        startidx = endidx - lenSTA
        STAsum += inputlist[startidx:endidx]
        numspikes += 1
    return numspikes,STAsum

def get_phaseimage( phaselist, phasemax, dt,\
    overlay = False, rasterwidth = 10, rasterheight = 1 ):
    phaseim = None
    for resplist in phaselist:
        ## each spike 'tick' is horizontal with dimensions rasterwidth x rasterheight
        phaseim_line = zeros( (int(phasemax/dt)*rasterheight, rasterwidth) )
        for phase in resplist:
            ## a horizontal line of rasterwidth for every spike
            ## of a given trial, respcycle and phase.
            row = int(phase/dt)*rasterheight
            phaseim_line[row:row+rasterheight,:] = 1.0
        if phaseim is None: phaseim = phaseim_line
        else:
            if overlay: phaseim += phaseim_line
            ## on numpy array(), axis=1, keep rows same, add cols.
            else: phaseim = append( phaseim, phaseim_line, axis=1)
    return phaseim

def plot_rasters(listof_rasterlists, runtime,\
    colorlist=['r','g','b'], labellist=['v1','v2','v3'], labels=True):
    fig = figure(facecolor='none')
    ax = fig.add_subplot(111)
    numrasterlists = len(listof_rasterlists)
    for rlistnum,rasterlist in enumerate(listof_rasterlists):
        numrasters = float(len(rasterlist))
        seglist = []
        for rnum,raster in enumerate(rasterlist):
            for t in raster:
                ## append a segment for a spike
                seglist.append(((t,rlistnum+rnum/numrasters),\
                    (t,rlistnum+(rnum+1)/numrasters)))
        ## plot the raster
        if labels:
            segs = collections.LineCollection(seglist,\
                 color=colorlist[rlistnum%len(colorlist)],\
                 label=labellist[rlistnum%len(labellist)])
        else:
            segs = collections.LineCollection(seglist,\
                 color=colorlist[rlistnum%len(colorlist)])
        ax.add_collection(segs)
    ax.set_xlim(0.0,runtime)
    if labels:
        ax.set_ylim(0,numrasterlists*1.3) # extra 0.3 height for legend
        biglegend()
    else:
        ax.set_ylim(0,numrasterlists)
    axes_labels(ax,'time (s)','spike raster trial#')
    title('Spike rasters', fontsize=24)

def crosscorr(x,y):
    """ pass numpy arrays x and y, so that element by element division & multiplication works.
    The older version of correlate() in numpy 1.4 gives only a scalar for tau=0. """
    return correlate(x,y)/(sqrt(correlate(x,x))*correlate(y,y))
    
def crosscorrgram(x, y, dt, halfwindow, starttime, endtime, norm='none'):
    """ pass arrays x and y of numtrials arrays of spike times.
    x[trialnum][tnum], y[trialnum][tnum]
    dt is the binsize, T = endtime-starttime
    Valid time length of the correlogram is from
    -halfwindow to +halfwindow.
    Analysis as per http://mulab.physiol.upenn.edu/crosscorrelation.html
    I further normalize by total number of spikes (Dhawale et al 2010 fig S5).
    I restrict the reference spike train x, between starttime+halfwindow to endtime-halfwindow;
    I restrict the compared spike train y, between starttime to endtime.
    
    For each spike in x, there is a sliding window of spikes in y.
    norm = 'overall': a la Ashesh et al 2010, divide by (total #spikes in all sliding windows over y).
    Above is same as dividing by (#spikesx * (mean #spikes in a sliding window of y)).
    norm = 'analogous': divide by (sqrt(#spikesx) * sqrt(#spikesy))
    Similar to dividing by sqrt(autocorrelationx)*sqrt(autocorrelationy)
    norm = 'ref': divide by (#spikesx)
    i.e. use the number of spikes in the reference spiketrain as the norm factor.
    Normalizes such that tau=0 value of auto-corr i.e. crosscorrgram(x,x,...) = 1
    norm = 'none': no division
    
    NOTE: mean is not subtracted from the two spike trains.
    To do that, first convert the list of spiketimes to a spike raster of 0s and 1s.
    Then subtract the respective means, and sum after element-wise multiplication.
    
    Finally, this function returns the average correlogram over all the trials.
    
    With a reference spiketrain x and normalization by #spikesy,
    the crosscorrgram becomes somewhat asymmetrical wrt x and y? """

    T = endtime-starttime
    xstarttime = starttime+halfwindow
    xendtime = endtime-halfwindow
    ##  div 2 and +1 to make number of bins odd
    bins = int(4*halfwindow/dt)/2 + 1
    centralbinnum = bins/2 ## integer division
    corrgramavg = array([0.0]*bins)
    ## x[trialnum][tnum]
    numtrials = len(x)
    corrnums = 0
    for trialnum in range(numtrials):
        xtrialnum = x[trialnum]
        if len(xtrialnum) == 0: continue
        spikenumx = 0
        spikenumy_allwindows = 0
        corrgram = array([0.0]*bins)
        for tx in xtrialnum:
            ## be careful, MOOSE inserts 0.0-s at the end of the fire times list!!!
            if tx<=xstarttime or tx>=xendtime: continue
            ## central bin is centered around t=0
            ## tx=ty falls in the center of the central bin.
            ystarttime = tx-halfwindow
            yendtime = tx+halfwindow
            spikenumx += 1
            for ty in y[trialnum]:
                if ty<=ystarttime or ty>=yendtime: continue
                binnum = round((ty-tx)/dt)+centralbinnum
                corrgram[binnum] += 1.0
                spikenumy_allwindows += 1
        ## if variable spikenumy_thistrial exists, add it.
        #if 'spikenumy_thistrial' in locals():
        #    spikenumy += spikenumy_thistrial

        ## Normalization:
        ## Divide by (total #spikes in all sliding windows over y).
        ## This is same as dividing by (#spikesx * (mean #spikes in sliding windows of y)).
        if norm=='overall':
            if spikenumy_allwindows>0:
                corrgram /= float(spikenumy_allwindows)
                corrgramavg += corrgram
                corrnums += 1
            #else: corrgram = [float('nan')]
        ## Divide by (sqrt(#spikesx) * sqrt(#spikesy))
        ## Similar to dividing by sqrt(autocorrelationx)*sqrt(autocorrelationy)
        elif norm=='analogous':
            spikenumy = 0
            for ty in y[trialnum]:
                if ty<starttime or ty>endtime: continue
                spikenumy += 1
            if spikenumx>0 and spikenumy>0:
                corrgram /= (sqrt(spikenumx)*sqrt(spikenumy))
                corrgramavg += corrgram
                corrnums += 1
        ## Divide by (#spikesx)
        ## Normalizes such that tau=0 value of auto-corr i.e. crosscorrgram(x,x,...) = 1
        elif norm=='ref':
            if spikenumx>0:
                corrgram /= spikenumx
                corrgramavg += corrgram
                corrnums += 1
        else:
            corrgramavg += corrgram
            corrnums += 1
    if corrnums==0: return array([nan]*bins)
    else: return corrgramavg/float(corrnums)

## --------------------------------------
## matplotlib stuff

def axes_off(ax,x=True,y=True):
    if x:
        for xlabel_i in ax.get_xticklabels():
            xlabel_i.set_visible(False)
            xlabel_i.set_fontsize(0.0)
    if y:
        for xlabel_i in ax.get_yticklabels():
            xlabel_i.set_fontsize(0.0)
            xlabel_i.set_visible(False)
    if x:
        for tick in ax.get_xticklines():
            tick.set_visible(False)
    if y:
        for tick in ax.get_yticklines():
            tick.set_visible(False)

def set_tick_widths(ax,tick_width):
    for tick in ax.xaxis.get_major_ticks():
        tick.tick1line.set_markeredgewidth(tick_width)
        tick.tick2line.set_markeredgewidth(tick_width)
    for tick in ax.xaxis.get_minor_ticks():
        tick.tick1line.set_markeredgewidth(tick_width)
        tick.tick2line.set_markeredgewidth(tick_width)
    for tick in ax.yaxis.get_major_ticks():
        tick.tick1line.set_markeredgewidth(tick_width)
        tick.tick2line.set_markeredgewidth(tick_width)
    for tick in ax.yaxis.get_minor_ticks():
        tick.tick1line.set_markeredgewidth(tick_width)
        tick.tick2line.set_markeredgewidth(tick_width)

def axes_labels(ax,xtext,ytext,adjustpos=False,fontsize=label_fontsize,xpad=None,ypad=None):
    ax.set_xlabel(xtext,fontsize=fontsize,labelpad=xpad)
    # increase xticks text sizes
    for label in ax.get_xticklabels():
        label.set_fontsize(fontsize)
    ax.set_ylabel(ytext,fontsize=fontsize,labelpad=ypad)
    # increase yticks text sizes
    for label in ax.get_yticklabels():
        label.set_fontsize(fontsize)
    if adjustpos:
        ## [left,bottom,width,height]
        ax.set_position([0.135,0.125,0.84,0.75])
    set_tick_widths(ax,axes_linewidth)

def biglegend(legendlocation='upper right',ax=None,fontsize=label_fontsize, **kwargs):
    if ax is not None:
        leg=ax.legend(loc=legendlocation, **kwargs)
    else:
        leg=legend(loc=legendlocation, **kwargs)
    # increase legend text sizes
    for t in leg.get_texts():
        t.set_fontsize(fontsize)

def beautify_plot(ax,x0min=True,y0min=True,
        xticksposn='bottom',yticksposn='left',xticks=None,yticks=None,
        drawxaxis=True,drawyaxis=True):
    """
    x0min,y0min control whether to set min of axis at 0.
    xticksposn,yticksposn governs whether ticks are at
    'both', 'top', 'bottom', 'left', 'right', or 'none'.
    xtickx/yticks is a list of ticks, else [min,max] is taken.
    Due to rendering issues,
    axes do not overlap exactly with the ticks, dunno why.
    """
    ax.get_yaxis().set_ticks_position(yticksposn)
    ax.get_xaxis().set_ticks_position(xticksposn)
    xmin, xmax = ax.get_xaxis().get_view_interval()
    ymin, ymax = ax.get_yaxis().get_view_interval()
    if x0min: xmin=0
    if y0min: ymin=0
    if xticks is None: ax.set_xticks([xmin,xmax])
    else: ax.set_xticks(xticks)
    if yticks is None: ax.set_yticks([ymin,ymax])
    else: ax.set_yticks(yticks)
    ### do not set width and color of axes by below method
    ### axhline and axvline are not influenced by spine below.
    #ax.axhline(linewidth=axes_linewidth, color="k")
    #ax.axvline(linewidth=axes_linewidth, color="k")
    ## spine method of hiding axes is cleaner,
    ## but alignment problem with ticks in TkAgg backend remains.
    for loc, spine in ax.spines.items(): # items() returns [(key,value),...]
        spine.set_linewidth(axes_linewidth)
        if loc == 'left' and not drawyaxis:
            spine.set_color('none') # don't draw spine
        elif loc == 'bottom' and not drawxaxis:
            spine.set_color('none') # don't draw spine
        elif loc in ['right','top']:
            spine.set_color('none') # don't draw spine
    ### alternate method of drawing axes, but for it,
    ### need to set frameon=False in add_subplot(), etc.
    #if drawxaxis:
    #    ax.add_artist(Line2D((xmin, xmax), (ymin, ymin),\
    #        color='black', linewidth=axes_linewidth))
    #if drawyaxis:
    #    ax.add_artist(Line2D((xmin, xmin), (ymin, ymax),\
    #        color='black', linewidth=axes_linewidth))
    ## axes_labels() sets sizes of tick labels too.
    axes_labels(ax,'','',adjustpos=False)
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)
    return xmin,xmax,ymin,ymax

def fig_clip_off(fig):
    ## clipping off for all objects in this fig
    for o in fig.findobj():
        o.set_clip_on(False)

## ------
## from https://gist.github.com/dmeliza/3251476#file-scalebars-py

# Adapted from mpl_toolkits.axes_grid2
# LICENSE: Python Software Foundation (http://docs.python.org/license.html)

from matplotlib.offsetbox import AnchoredOffsetbox
class AnchoredScaleBar(AnchoredOffsetbox):
    def __init__(self, transform, sizex=0, sizey=0, labelx=None, labely=None, loc=4,
                 pad=0.1, borderpad=0.1, sep=2, prop=None, label_fontsize=label_fontsize, color='k', **kwargs):
        """
        Draw a horizontal and/or vertical  bar with the size in data coordinate
        of the give axes. A label will be drawn underneath (center-aligned).

        - transform : the coordinate frame (typically axes.transData)
        - sizex,sizey : width of x,y bar, in data units. 0 to omit
        - labelx,labely : labels for x,y bars; None to omit
        - loc : position in containing axes
        - pad, borderpad : padding, in fraction of the legend font size (or prop)
        - sep : separation between labels and bars in points.
        - **kwargs : additional arguments passed to base class constructor
        """
        from matplotlib.patches import Rectangle
        from matplotlib.offsetbox import AuxTransformBox, VPacker, HPacker, TextArea, DrawingArea
        bars = AuxTransformBox(transform)
        if sizex:
            bars.add_artist(Rectangle((0,0), sizex, 0, fc="none", linewidth=axes_linewidth, color=color))
        if sizey:
            bars.add_artist(Rectangle((0,0), 0, sizey, fc="none", linewidth=axes_linewidth, color=color))

        if sizex and labelx:
            textareax = TextArea(labelx,minimumdescent=False,textprops=dict(size=label_fontsize,color=color))
            bars = VPacker(children=[bars, textareax], align="center", pad=0, sep=sep)
        if sizey and labely:
            ## VPack a padstr below the rotated labely, else label y goes below the scale bar
            ## Just adding spaces before labely doesn't work!
            padstr = '\n '*len(labely)
            textareafiller = TextArea(padstr,textprops=dict(size=label_fontsize/3.0))
            textareay = TextArea(labely,textprops=dict(size=label_fontsize,rotation='vertical',color=color))
            ## filler / pad string VPack-ed below labely
            textareayoffset = VPacker(children=[textareay, textareafiller], align="center", pad=0, sep=sep)
            ## now HPack this padded labely to the bars
            bars = HPacker(children=[textareayoffset, bars], align="top", pad=0, sep=sep)

        AnchoredOffsetbox.__init__(self, loc, pad=pad, borderpad=borderpad,
                                   child=bars, prop=prop, frameon=False, **kwargs)

def add_scalebar(ax, matchx=True, matchy=True, hidex=True, hidey=True, \
    label_fontsize=label_fontsize, color='k', **kwargs):
    """ Add scalebars to axes

    Adds a set of scale bars to *ax*, matching the size to the ticks of the plot
    and optionally hiding the x and y axes

    - ax : the axis to attach ticks to
    - matchx,matchy : if True, set size of scale bars to spacing between ticks
                    if False, size should be set using sizex and sizey params
    - hidex,hidey : if True, hide x-axis and y-axis of parent
    - **kwargs : additional arguments passed to AnchoredScaleBars

    Returns created scalebar object
    """
    def f(axis):
        l = axis.get_majorticklocs()
        return len(l)>1 and (l[1] - l[0])
    
    if matchx:
        kwargs['sizex'] = f(ax.xaxis)
        kwargs['labelx'] = str(kwargs['sizex'])
    if matchy:
        kwargs['sizey'] = f(ax.yaxis)
        kwargs['labely'] = str(kwargs['sizey'])
        
    sb = AnchoredScaleBar(ax.transData, label_fontsize=label_fontsize, color=color, **kwargs)
    ax.add_artist(sb)

    if hidex : ax.xaxis.set_visible(False)
    if hidey : ax.yaxis.set_visible(False)

    return sb

## from https://gist.github.com/dmeliza/3251476#file-scalebars-py -- ends
## ------

## matplotlib stuff ends
## -----------------------------------------

def plotSpikes(firetimes, runtime, plotdt):
    firetimes = array(firetimes)
    # MOOSE often inserts one or two spiketime = 0.0 entries when storing spikes, so discount those:
    firetimes = firetimes[ where(firetimes>0.0)[0] ]
    firetimes = firetimes[ where(diff(firetimes)>2*plotdt)[0] ] # Take the falling edge of every threshold crossing.
    firearray = zeros(int(round(runtime/plotdt)),dtype=int8) # 1D array of type int8
    firelen = len(firearray)
    for firetime in firetimes:
        firearray[int(round(firelen*float(firetime))/runtime)] = 1
    return firearray

def plotBins(firetimes, numbins, runtime, settletime):
    binlist = [0]*numbins
    firetimes = array(firetimes)
    ## MOOSE often inserts one or two spiketime = 0.0 entries
    ## when storing spikes, so discount those:
    firetimes = firetimes[ where(firetimes>0.0)[0] ]
    for firetime in firetimes:
        if firetime>=settletime:
            ## The small number has been added to the Dr to ensure no index out of range errors
            ## Nothing to do about causality here:
            ## while plotting, keep bintime to right edge to make it causal
            binnum = int((firetime-settletime)/(runtime-settletime+0.0001)*numbins)
            binlist[binnum] += 1
    return [binspikes/((runtime-settletime)/float(numbins)) for binspikes in binlist] # return firing rate in Hz

def plotOverlappingBins(firetimes, numbins, time_period, settletime, bin_width_time):
    """
    Firing rate in overlapping bins (moving average).
    numbins # of bins in the time (settletime) to (time_period+settletime)
    Assumes periodic/wrapped boundary conditions with period=time_period.
    This way the end bins are accurate,
    else they will not have data to one end and show lower firing rates.
    Typically, adjust settletime to bin
    only the first or second respiratory cycle.
    """
    CAUSAL = True
    binlist = [0]*numbins
    firetimes = array(firetimes)
    ## MOOSE often inserts one or two spiketime = 0.0 entries
    ## when storing spikes, so discount those:
    firetimes = firetimes[ where(firetimes>0.0)[0] ]
    bindt = time_period/float(numbins)
    ## if CAUSAL, take spikes only to the left of bin centre_times.
    if CAUSAL: centre_times = arange(bindt, time_period+bindt/2.0, bindt)
    else: centre_times = arange(bindt/2, time_period, bindt)
    bin_half_t = bin_width_time/2.0
    rightmost_t = time_period
    for firetime in firetimes:
        ## The end bins will not show correct firing rate!
        if firetime>=settletime and firetime<(settletime+time_period):
            firetime -= settletime
            ## Each firetime is in multiple bins depending on bin_width_time
            for binnum,bin_centre_t in enumerate(centre_times):
                ## if CAUSAL, take spikes only to the left of bin centre_times.
                if CAUSAL:
                    bin_left = bin_centre_t - bin_width_time
                    bin_right = bin_centre_t
                else:
                    bin_left = bin_centre_t - bin_half_t
                    bin_right = bin_centre_t + bin_half_t
                if firetime >= bin_left and firetime < bin_right:
                    binlist[binnum] += 1
                ## Next lines implement circularity of firetimes
                if bin_left < 0 and firetime >= (bin_left+rightmost_t):
                    binlist[binnum] += 1
                if bin_right > rightmost_t and firetime < (bin_right-rightmost_t):
                    binlist[binnum] += 1
    return [float(binspikes)/bin_width_time for binspikes in binlist] # return firing rate in Hz

def calcFreq(timeTable, runtime, settletime, plotdt, threshold, spiketable):
    # input: if spiketable is True: timeTable has spike times: i.e. a MOOSE table which has stepMode = TAB_SPIKE
    # input: if spiketable is False: timeTable has Vm-s: i.e. a MOOSE table which has stepMode = TAB_BUF
    # output: (meanrate, meanrate2, events)
    # output: events is a list of times of falling edges of 'spikes' separated by at least 2*eventdt.
    # output: meanrate2 is just #spikes/time
    # output: meanrate is mean of 1/inter-spike-interval (removing very short ISIs)
    tablenumpy = array(timeTable) # convert the MOOSE table into a numpy array

    if spiketable: # timeTable has spike times
        # only those spike times which are after settle time.
        # it is important to do this even is settletime == 0.0,
        # since MOOSE inserts spurious t=0.0 spitketime entries in a spike table.
        events = tablenumpy[ where(tablenumpy>settletime)[0] ]
    else: # timeTable has Vm-s
        cutout = tablenumpy[ int(settletime/plotdt): ] # cutout only those spike times which are after settle time.
        if len(cutout) <= 0:
            events = []
        else:
            thresholded = where(cutout>threshold)[0] # gives indices whereever cutout > THRESHOLD
            # where difference between two adjacent indices in thresholded > 2
            # i.e. takes falling edge of every threshold crossing
            # THIS IS UNLIKE SPIKETABLE ABOVE WHICH TAKES RISING EDGE!
            take = where(diff(thresholded)>2)[0] # numpy's where and diff 
            indices = thresholded[ take ] # indexed by a list! works for ndarray only, not for usual python lists.
            # numpy multiplication of array by scalar -- very different from python list multiplication by integer!!!
            events = indices*plotdt + settletime

    # calculate mean firing rate as 1/inter-spike-interval (removing very short ISIs)
    if len(events)>1: # at least two events needed!
        firingRateList = array([])
        for i in range(len(events)-1):
            firingtimespan = (events[i+1]-events[i])
            firingRateList = append(firingRateList,1.0/firingtimespan)
        ############ Have to filter out the APs which have very closely spaced double firings :
        ############ typically happens for close to zero currents like 0.15nA etc.
        ## Keep removing all firing rate entries which are greater than twice the min.
        while firingRateList.max() > 2*firingRateList.min():
            firingRateList = delete(firingRateList,firingRateList.argmax())
        ################### Finally calculate the actual mean
        meanrate = array(firingRateList).mean() # 1/s = Hz
    else:
        meanrate = 0

    # mean firing rate as #spikes/time
    meanrate2 = len(events)/(runtime-settletime) # Hz

    return (meanrate, meanrate2, events)

def minimum_distance(a,b,p):
    """ a,b,p are vectors.
    Given two end-points a and b of a line segment,
    find the minimum distance from p to the the line segment.
    points could be 2D or 3D."""
    a = array(a)
    b = array(b)
    p = array(p)
    ## length sq of line segment
    len_sq = norm(b-a)**2
    if len_sq==0: return norm(p-a)
    ## take infinte line as c = a + t*(b-a)
    ## find t where the point p drops a normal to line
    t = dot(b-a,p-a)/len_sq
    ## point before a
    if t<0: return norm(p-a)
    ## point after b
    elif t>1.0: return norm(p-b)
    else: return norm( p - (a+t*(b-a)) )

def outcode(x,y,xmin,ymin,xmax,ymax):
    """ utility function for the Cohen-Sutherland clipper below. """
    outcode = 0x0 # inside
    if y>ymax: outcode |= 0x1 # top
    elif y<ymin: outcode |= 0x2 # bottom
    if x>xmax: outcode |= 0x4 # right
    elif x<xmin: outcode |= 0x8 # left
    
    return outcode

def clip_line_to_rectangle(x1,y1,x2,y2,xmin,ymin,xmax,ymax):
    """ clip line segment from (x1,y1) to (x2,y2) inside a rectange. 
    2D points only. Cohen-Sutherland algorithm (Wikipedia)"""
    outcode1 = outcode(x1,y1,xmin,ymin,xmax,ymax)
    outcode2 = outcode(x2,y2,xmin,ymin,xmax,ymax)
    accept = False
    done = False
    while not done:
        if not (outcode1 | outcode2): # segment fully inside
            accept = True
            done = True
        elif outcode1 & outcode2:   # segment fully outside i.e.
                                    # both points are to the left/right/top/bottom of rectangle
            done = True
        else:
            if outcode1>0: outcode_ex = outcode1
            else: outcode_ex = outcode2
            if outcode_ex & 0x1 :  # top
                x = x1 + (x2 - x1) * (ymax - y1) / (y2 - y1)
                y = ymax
            elif outcode_ex & 0x2: # bottom
                x = x1 + (x2 - x1) * (ymin - y1) / (y2 - y1)
                y = ymin
            elif outcode_ex & 0x4: # right
                y = y1 + (y2 - y1) * (xmax - x1) / (x2 - x1)
                x = xmax
            else: # left
                y = y1 + (y2 - y1) * (xmin - x1) / (x2 - x1)
                x = xmin
 
            ## get the new co-ordinates of the line
            if (outcode_ex == outcode1):
                x1 = x
                y1 = y
                outcode1 = outcode(x1, y1, xmin, ymin, xmax, ymax)
            else:
                x2 = x
                y2 = y
                outcode2 = outcode(x2, y2, xmin, ymin, xmax, ymax)

    return accept,x1,y1,x2,y2

## copied savitsky_golay from scipy cookbook: http://www.scipy.org/Cookbook/SavitzkyGolay
def savitzky_golay(y, window_size, order, deriv=0):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techhniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv]
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m, y, mode='valid')

def circular_convolve(x,y,n):
    """
    From: http://www.dspguru.com/dsp/tutorials/a-little-mls-tutorial
    I need to compute the circular cross-correlation of y and x
                       N-2
    Ryx[n] = 1/(N-1) * SUM{ y[i] * x[i-n] }
                       i=0
    python indexing (negative indices) automatically gives you circularity!!!
    """
    N = len(y)
    return array( 1.0/(N-1)*sum( [ y[i]*x[i-n] for i in range(N-1) ] ) )
    
def circular_correlate(x,y):
    N = len(y)
    return array( [ circular_convolve(x,y,n) for n in range(len(y)) ] )

############################ Information theory functions ###############################

def binary_entropy(p0,p1):
    H = 0.0
    if p0>0.0:
        H -= p0*log2(p0)
    if p1>0.0:
        H -= p1*log2(p1)
    return H

def makestring(intlist):
    s = ''
    for val in intlist:
        s += str(val)
    return s

## yield defines a generator
def find_substr_endchars(mainstr,substr,delay=0):
    """
    Returns a list (generator) of characters
    each following a 'substr' occurence in 'mainstr' after 'delay'.
    """
    ## don't use regular expressions re module, which finds only non-overlapping matches
    ## we want to find overlapping matches too.
    substrlen = len(substr)
    while True:
        idx = mainstr.find(substr)
        ## find returns -1 if substr not found
        if idx != -1:
            endcharidx = idx+substrlen+delay
            if endcharidx<len(mainstr):
                yield mainstr[endcharidx]
            else: # reached end of string
                break
            ## chop the mainstr just after the start of substr,
            ## not after the end, as we want overlapping strings also
            mainstr = mainstr[idx+1:]
        else: # substr not found
            break

## yield defines a generator
def find_substrs12_endchars(sidestr,mainstr,substr1,substr2,delay1=0,delay2=0):
    """
    Returns a list (generator) of characters from mainstr Y,
    each following a substr1 occurence in sidestr X after delay1,
    and following a substr2 occurence in mainstr Y after delay2.
    """
    ## don't use regular expressions re module, which finds only non-overlapping matches
    ## we want to find overlapping matches too.
    substr2len = len(substr2)
    substr1len = len(substr1)
    abs_idx1 = 0 ## mainstr is getting chopped, but we maintain abs index on sidestr
    while True:
        idx2 = mainstr.find(substr2)
        ## find returns -1 if substr2 not found
        if idx2 != -1:
            endcharidx2 = idx2+substr2len+delay2
            ### NOTE: abs_startidx1 is one earlier than definition!!! I think necessary for causality.
            ## put +1 below to switch to definition in Quinn et al 2010
            abs_startidx1 = abs_idx1 + endcharidx2 - substr1len-delay1
            if endcharidx2<len(mainstr): # mainstr Y has characters left?
                if abs_startidx1 >= 0: # sidestr X has sufficient chars before?
                    ## sidestr has substr1 before the char to be returned? and mainstr is not over
                    ## IMP: below if's first term is the only place directed info enters.
                    ## Remove first term below and you get just the entropy of mainstr Y: VERIFIED.
                    #print sidestr[abs_startidx1:abs_startidx1+substr1len], substr1, abs_startidx1
                    if sidestr[abs_startidx1:abs_startidx1+substr1len]==substr1:
                        yield mainstr[endcharidx2]
            else: # reached end of string
                break
            ## chop the mainstr just after the start of substr2,
            ## not after the end, as we want overlapping strings also
            mainstr = mainstr[idx2+1:]
            ## don't chop sidestr as substr1len may be greater than substr2len
            ## in the next iteration, idx2 will be relative, but for sidestr we maintain abs_idx1
            abs_idx1 += idx2+1
        else: # substr2 not found
            break

def calc_entropyrate(spiketrains,markovorder,delay=0):
    """
    spiketrains is a list of spiketrain = [<0|1>,...]
    should be int-s, else float 1.0-s make the str-s go crazy!
    J = markovorder >= 1. Cannot handle non-Markov i.e. J=0 presently!!
    Returns entropy rate, assuming markov chain of order J :
    H(X_{J+1}|X_J..X_1).
    Assumes spike train bins are binary i.e. 0/1 value in each timebin.
    delay is self delay of effect of X on X.
    NOTE: one time step i.e. causal delay is permanently present, 'delay' is extra.
    """
    Hrate = 0.0
    N = 0
    if markovorder>0:
        ## create all possible binary sequences priorstr=X_1...X_J
        ## convert integer to binary repr str of length markovorder padded with zeros (=0)
        reprstr = '{:=0'+str(markovorder)+'b}'
        priorstrs = [ reprstr.format(i) for i in range(int(2**markovorder)) ]
    else:
        ## return numpy nan if markovorder <= 0
        return nan
    ## Convert the list of timebins to a string of 0s and 1s.
    ## Don't do it in loops below, else the same op is repeated len(priorstrs) times.
    ## Below conversion is quite computationally expensive.
    mcs = []
    for spiketrain in spiketrains:
        ## A generator expression is given as argument to makestring
        mcs.append(makestring(val for val in spiketrain))

    ## Calculate entropy for each priorstr, and sum weighted by probability of each priorstr
    for priorstr in priorstrs:
        num1s = 0
        num0s = 0
        for mc in mcs:
            for postchar in find_substr_endchars(mc,priorstr,delay):
                if int(postchar): # if the character just after priorstr is nonzero i.e. 1
                    num1s += 1
                else:
                    num0s += 1
        N_givenprior = float(num1s + num0s)
        ## H(X|Y) = \sum p(Y=y)*H(X|Y=y) ; the normalization by N is done at the end
        ## p(Y=y) = N_givenprior/N where N is total after all loops
        if N_givenprior>0:
            Hrate += N_givenprior * binary_entropy(num0s/N_givenprior,num1s/N_givenprior)
            N += N_givenprior
    if N!=0: Hrate = Hrate/N
    return Hrate

def calc_dirtinforate(spiketrains1,spiketrains2,markovorder1,markovorder2,delay1=0,delay2=0):
    """
    Returns directed information rate from spiketrains1 X to spiketrains2 Y.
    Returns directed information rate (lim_{n->/inf} 1/n ...),
    assuming train2 as markov chain of order K,
    and train1 affecting it with markov order J.
    I(X^n->Y^n) = H( Y_{J+1} | Y_J..Y_1 ) - H( Y_{L} | Y^{L-1}_{L-J} X^{L-1}_{L-K} ),
    where L = max(J,K).
    NOTE: I have changed X^{L}_{L-K-1} in definition above to X^{L-1}_{L-K} for causality!
    Assumes spike train bins are binary i.e. integer 0/1 value in each timebin.
    spiketrains1 and 2 are each a list of spiketrain = [<0|1>,...]
    should be int-s, else float 1.0-s make the str-s go crazy!
    dimensions of both must be the same.
    J = markovorder1 >= 1. Cannot handle non-Markov i.e. J=0 presently!!
    K = markovorder2 >= 1. Cannot handle non-Markov i.e. K=0 presently!!
    Keep J,K<5, else too computationally intensive.
    The prior substrings are searched delay1 and delay2 before Y_n in trains 1 and 2.
    delay1 is lateral/side delay of effect of X on Y,
    delay2 is self/main delay of effect of Y on Y.
    NOTE: one time step i.e. causal delay is permanently present, delay1 and 2 are extra.
    """
    dirtIrate_term2 = 0.0
    N = 0
    ## for the 'cause' spike train
    if markovorder1>0:
        ## create all possible binary sequences priorstr=X_1...X_J
        ## convert integer to binary repr str of length markovorder padded with zeros (=0)
        reprstr = '{:=0'+str(markovorder1)+'b}'
        priorstrs1 = [ reprstr.format(i) for i in range(int(2**markovorder1)) ]
    else:
        ## return numpy nan if markovorder <= 0
        return nan
    ## for the 'effect' spike train
    if markovorder2>0:
        ## create all possible binary sequences priorstr=X_1...X_K
        ## convert integer to binary repr str of length markovorder padded with zeros (=0)
        reprstr = '{:=0'+str(markovorder2)+'b}'
        priorstrs2 = [ reprstr.format(i) for i in range(int(2**markovorder2)) ]
    else:
        ## return numpy nan if markovorder <= 0
        return nan

    ## Convert the list of timebins to a string of 0s and 1s.
    ## Don't do it in loops below, else the same op is repeated len(priorstrs) times.
    ## Below conversion is quite computationally expensive.
    mcs1 = []
    for spiketrain in spiketrains1:
        ## A generator expression is given as argument to makestring
        mcs1.append(makestring(val for val in spiketrain))
    mcs2 = []
    for spiketrain in spiketrains2:
        ## A generator expression is given as argument to makestring
        mcs2.append(makestring(val for val in spiketrain))

    ## Calculate entropy for each combo of priorstr 1 & 2,
    ## and sum weighted by probability of each combo
    for priorstr1 in priorstrs1:
        for priorstr2 in priorstrs2:
            num1s = 0
            num0s = 0
            for chaini,mc1 in enumerate(mcs1):
                mc2 = mcs2[chaini]
                for postchar in find_substrs12_endchars(mc1,mc2,priorstr1,priorstr2,delay1,delay2):
                    ## if the character just after priorstr1 & priorstr2, is nonzero i.e. 1
                    if int(postchar):
                        num1s += 1
                    else:
                        num0s += 1
            N_givenpriors = float(num1s + num0s)
            ## H(Y|Y^X^) = \sum p(Y^=y^)*H(Y|Y^=y^,X^=x^) ;
            ## the normalization by N is done at the end
            ## p(Y^=y^,X^=x^) = N_givenpriors/N where N is total after all loops
            if N_givenpriors>0:
                dirtIrate_term2 += N_givenpriors * \
                    binary_entropy(num0s/N_givenpriors,num1s/N_givenpriors)
                N += N_givenpriors
    if N!=0: dirtIrate_term2 = dirtIrate_term2/N

    ## H( Y_{J+1} | Y_J..Y_1 )
    dirtIrate_term1 = calc_entropyrate(spiketrains2,markovorder2,delay2)
    ## I(X^n->Y^n) = H( Y_{J+1} | Y_J..Y_1 ) - H( Y_{L} | Y^{L-1}_{L-J} X^{L-1}_{L-K} )
    dirtIrate = dirtIrate_term1 - dirtIrate_term2
    return dirtIrate

def get_spiketrain_from_spiketimes(\
        spiketimes,starttime,timerange,numbins,warnmultiple=True,forcebinary=True):
    """ bin number of spikes from starttime to endtime, into dt bins.
    if warnmultiple, warn if multiple spikes are binned into a single bin.
    if forcebinary, set multiple spikes in a bin to 1.
    """
    ## important to make these int, else spikestrs in entropy calculations go haywire!
    spiketrain = zeros(numbins,dtype=int)
    for spiketime in spiketimes:
        spiketime_reinit = spiketime-starttime
        if 0.0 < spiketime_reinit < timerange:
            binnum = int(spiketime_reinit/timerange*numbins)
            spiketrain[binnum] += 1
    if forcebinary or warnmultiple:
        multiplespike_indices = where(spiketrain>1)[0] # spiketrain must be a numpy array()
        if len(multiplespike_indices)>0:
            ## if non-empty number of multiple spikes indices, set to 1, print warning
            if forcebinary:
                spiketrain[multiplespike_indices] = 1
            if warnmultiple:
                ## do not print warnings if user has turned them off
                print "There are more than 1 spikes in", \
                    len(multiplespike_indices), "number of bins."
                if forcebinary: print "Have forced them all to be 1."
    return spiketrain

##############################################################################
