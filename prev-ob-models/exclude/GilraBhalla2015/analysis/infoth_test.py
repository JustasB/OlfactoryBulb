
import sys
sys.path.extend([".."])

from data_utils import *

def check_substr():
    """check find_substr_endchars"""
    a = '101010100'
    print a
    print [ch for ch in find_substr_endchars(a,'1010')]

## test entropy rate

def IID_entropyrate():
    """equal prob of 1s and 0s IID. Entropy rate should be 1 bit per use."""
    spiketrains = [[int(uniform(0,1)+0.5) for i in range(1000)]]
    print "Equal prob of 1s and 0s IID. Hrate =",calc_entropyrate(spiketrains,1)

def markov_order1_entropyrate():
    """markov process of order 1
    p(0)=0.6, p(1)=0.4"""
    spiketrains = []
    for j in range(10):
        mc = [1]
        for i in range(10000):
            if mc[-1]:
                p = uniform(0,1)
                if p<0.25: mc.append(1)
                else: mc.append(0)
            else:
                p = uniform(0,1)
                if p<0.5: mc.append(1)
                else: mc.append(0)
        spiketrains.append(mc)
    entropyrates = []
    ## careful: don't go to higher orders.. machine stalls badly..
    for i in range(10):
        entropyrates.append(calc_entropyrate(spiketrains,i))
    figure(facecolor='w')
    plot(entropyrates,'r-,')

def plot_table(rasters,rowlabels,collabels,data,cellcolours,titlestr,figfilename):
    ## 'plot' a table
    fig = figure(figsize=(8, 6), dpi=100)
    ax = fig.add_axes([0.14, 0.75, 0.85, 0.1])
    axes_off(ax)
    ## loop over rasters in reverse order, as they are plotted from bottom upwards
    for rasteri,raster in enumerate(rasters[::-1]):
        raster = array(raster)
        ## find out indices of 1-s and plot them:
        rasterindices = where(raster==1)[0]
        ax.plot(rasterindices,zeros(len(rasterindices))+rasteri,'|k',\
            markersize=20, markeredgewidth='2') # | is the marker
    ax.set_ylim(-0.5,rasteri+0.5)
    dirtItable = ax.table(cellText=data, cellColours=cellcolours, rowLoc='right',\
        rowLabels=rowlabels, colLabels=collabels, colLoc='center', loc='bottom')
    table_props = dirtItable.properties()
    table_cells = table_props['child_artists']
    for cell in table_cells:
        cell.set_height(1.35)
        cell.set_fontsize(18)
    ax.set_title(titlestr,fontsize=14)
    ## tight_layout() doesn't seem to work with table
    #fig.tight_layout()
    #fig.savefig(figfilename,dpi=300)

def plot_dirtIrates(spiketrains1,spiketrains2,delay1,delay2,\
                    dirtIcutoff,filename='none.svg'):
    collabels = ['Order 1','2','3']
    rowlabels = ['Delay 0','1','2','3']
    dirtIs = []
    cellcolours = []
    for measure_delay in range(4):
        print "delay =",measure_delay
        dirtIorders = []
        cellcoloursorders = []
        for i in range(1,4):
            print "markovorder =",i
            causal_dirtI = calc_dirtinforate(spiketrains1,spiketrains2,\
                i,i,measure_delay,measure_delay)
            oppcausal_dirtI = calc_dirtinforate(spiketrains2,spiketrains1,\
                i,i,measure_delay,measure_delay)
            dirtIstr = '    causal = {:1.3f}\nopp causal = {:1.3f}'\
                .format(causal_dirtI,oppcausal_dirtI)
            print dirtIstr
            dirtIorders.append(dirtIstr)
            if causal_dirtI>dirtIcutoff: cellcoloursorders.append('r')
            else: cellcoloursorders.append('w')
        dirtIs.append(dirtIorders)
        cellcolours.append(cellcoloursorders)
    titlestr = "Copy IID spike train 1 to 2 with delays = "\
            +str(delay1)+" & "+str(delay2)+\
            "\n Measure with markov order & common delay as below"+\
            "\n causal is 1->2, opp causal is 2->1"
    plot_table([spiketrains1[0][0:100],spiketrains2[0][0:100]],\
        rowlabels,collabels,dirtIs,cellcolours,titlestr,filename)


def copycat_dirtIrate(delay=0):
    """Test directed information: trains1 has full causal exc effect on trains2.
    Copies spiketrain1 to 2 with causal delay i.e. delay+1.
    spiketrains and spiketrains1 is IID, p(0)=0.5."""
    spiketrains = []
    spiketrains1 = []
    spiketrains2 = []
    for j in range(10):
        ## need delay+1 num of copies of 1, to access mc1[-delay-1] below
        ## maintain same length of spike trains, hence common start length
        mc = [1]*(delay+1)
        mc1 = [1]*(delay+1)
        mc2 = [1]*(delay+1)
        for i in range(10000):
            ## mc2 fully depends on mc1 with delay
            if mc1[-delay-1]:
                mc2.append(1)
            else:
                mc2.append(0)
            ## mc1 is IID
            p = uniform()
            if p<0.5: mc1.append(1)
            else: mc1.append(0)
            ## mc is IID
            p = uniform()
            if p<0.5: mc.append(1)
            else: mc.append(0)
        spiketrains.append(mc)
        spiketrains1.append(mc1)
        spiketrains2.append(mc2)

    collabels = ['Order 1','2','3']
    rowlabels = ['Delay 0','1','2','3']
    dirtIs = []
    cellcolours = []
    for measure_delay in range(4):
        print "delay =",measure_delay
        dirtIorders = []
        cellcoloursorders = []
        for i in range(1,4):
            print "order =",i
            causal_dirtI = calc_dirtinforate(spiketrains1,spiketrains2,\
                i,i,measure_delay,measure_delay)
            oppcausal_dirtI = calc_dirtinforate(spiketrains2,spiketrains1,\
                i,i,measure_delay,measure_delay)
            acausal_dirtI = calc_dirtinforate(spiketrains,spiketrains2,\
                i,i,measure_delay,measure_delay)
            dirtIstr = '    causal = {:1.3f}\nopp causal = {:1.3f}\n    acausal = {:1.3f}'\
                .format(causal_dirtI,oppcausal_dirtI,acausal_dirtI)
            print dirtIstr
            dirtIorders.append(dirtIstr)
            if causal_dirtI>0.9: cellcoloursorders.append('r')
            else: cellcoloursorders.append('w')
        dirtIs.append(dirtIorders)
        cellcolours.append(cellcoloursorders)
    titlestr = "Copy IID spike train 1 to 2 with causal delay = "+str(delay)+\
            "\n Measure with markov order & common delay as below"+\
            "\n causal is 1->2, opp causal is 2->1, acausal is IID->2"
    plot_table([spiketrains1[0][0:100],spiketrains2[0][0:100]],\
        rowlabels,collabels,dirtIs,cellcolours,titlestr,'copycat_mydefn.svg')

def partialcopycat_dirtIrate(delay1=0,delay2=0,inh=True):
    """Test directed information: trains1 has partial causal exc effect on trains2.
    Copies partially spiketrain1 to 2, with 1 delayed by delay1 and 2 delayed by delay2,
    (these are causal delays, hence delay1+1 and delay2+1)
    and flipping if inh."""
    spiketrains1 = []
    spiketrains2 = []
    if inh:
        one=0
        zer=1
    else:
        one=1
        zer=0
    for j in range(10):
        ## need delay+1 num of copies of 1, to access mc1[-delay1-1] below
        ## maintain same length of spike trains, hence common start length
        delay = max(delay1,delay2)
        mc1 = [1]*(delay+1)
        mc2 = [1]*(delay+1)
        for i in range(10000):
            ## mc2 depends on mc1 and mc2,
            if mc1[-delay1-1] and mc2[-delay2-1]: # both are 1
                mc2.append(one)
            elif mc1[-delay1-1] or mc2[-delay2-1]: # elif either is 1
                p = uniform(0,1)
                if p<0.75: mc2.append(one)
                else: mc2.append(zer)
            else: # both are zero
                mc2.append(zer)
            ## mc1 is IID
            p = uniform()
            if p<0.5: mc1.append(1)
            else: mc1.append(0)
        spiketrains1.append(mc1)
        spiketrains2.append(mc2)

    plot_dirtIrates(spiketrains1,spiketrains2,delay1,delay2,0.1,'partialcopycat_mydefn.svg')

def bicausal_dirtIrate(delay1=0,delay2=0,inh=True):
    """Test directed information: bidirectional exc/inh effect."""
    spiketrains1 = []
    spiketrains2 = []
    if inh:
        one=0
        zer=1
    else:
        one=1
        zer=0
    for j in range(10):
        ## need delay+1 num of copies of 1, to access mc1[-delay1-1] below
        ## maintain same length of spike trains, hence common start length
        delay = max(delay1,delay2)
        mc1 = [1]*(delay+1)
        mc2 = [1]*(delay+1)
        for i in range(10000):
            ## mc1 and mc2 depend on mc1 and mc2, but asymmetrically
            ## delay is on both here!
            if mc1[-delay1-1] and mc2[-delay2-1]: # both are 1
                p = uniform(0,1)
                if p<0.85: mc2.append(one)
                else: mc2.append(zer)
                p = uniform(0,1)
                if p<0.65: mc1.append(one)
                else: mc1.append(zer)
            elif mc1[-delay1-1] or mc2[-delay2-1]: # elif either is 1
                p = uniform(0,1)
                if p<0.6: mc2.append(one)
                else: mc2.append(zer)
                p = uniform(0,1)
                if p<0.5: mc1.append(one)
                else: mc1.append(zer)
            else: # both are zero
                p = uniform(0,1)
                if p<0.8: mc2.append(zer)
                else: mc2.append(one)
                p = uniform(0,1)
                if p<0.6: mc1.append(zer)
                else: mc1.append(one)
        spiketrains1.append(mc1)
        spiketrains2.append(mc2)

    plot_dirtIrates(spiketrains1,spiketrains2,delay1,delay2,0.01,'bicausal_mydefn.svg')


if __name__ == "__main__":
    #copycat_dirtIrate()
    partialcopycat_dirtIrate(delay1=1,delay2=2,inh=False)
    partialcopycat_dirtIrate(delay1=1,delay2=2,inh=True)
    #bicausal_dirtIrate(delay1=1,delay2=2,inh=True)
    show()
