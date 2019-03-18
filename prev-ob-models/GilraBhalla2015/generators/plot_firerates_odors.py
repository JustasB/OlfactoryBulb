import sys, pickle

sys.path.extend(["..","../networks","../simulations"])
from networkConstants import *
from stimuliConstants import *
from simset_odor import *
import generate_firerates_odors as odor_gen

from pylab import * # part of matplotlib that depends on numpy but not scipy

from data_utils import * # has axes_labels()

## USAGE: python2.6 plot_firerates_odors.py <firerates_filename>

#### We have dual exponential functions as impulse response / kernels
#### for respiration, odor A and odor B (each is different for every glomerulus)
#### These kernels are convolved with the respiration pulse
#### to generate the firing rate which is fed to a Poissonian generator.

RUNTIME = REALRUNTIME + SETTLETIME

## time points for the firing rate which is read from a pickled file
firingtsteps = arange(0,RUNTIME+1e-10,FIRINGFILLDT)# include the last RUNTIME point also.
numt = len(firingtsteps)
extratime = arange(0,3*RESPIRATION+1e-10,FIRINGFILLDT)
randompulsetime = arange(0,PULSE_RUNTIME+1e-10,FIRINGFILLDT)

fname = sys.argv[1]
f = open(fname,'r')
frateOdorList,fratePulseList,randomPulseList, \
randomPulseStepsList,randomResponseList,kernels \
= pickle.load(f)
f.close()

def odor_and_air_stimuli():
    ### mitral and PG odor ORNs firing files
    #for glomnum in range(NUM_GLOMS):
    #    fig = figure(facecolor='w')
    #    ax = fig.add_subplot(111)
    #    title('Glom '+str(glomnum)+': odor morphs', fontsize=30)
    #    frateperglomList = []
    #    for odoridx,(odorA, odorB) in enumerate(inputList):
    #        frate = frateOdorList[glomnum][odoridx]
    #        plot(firingtsteps, frate, linewidth=2, color=(odorA,odorB,0), marker=',')
    #    axes_labels(ax,'time (s)','ORN firing rate (Hz)')
    ## show ORN inputs for lat inh 
    ## odor A
    fig = figure(facecolor='w')
    ax = fig.add_subplot(111)
    title('Odor A: glom0 (r); glom1 (g); glom2 (b)', fontsize=30)
    for glomnum in [0,1,2]:
        frateperglomList = []
        odoridx = 5
        frate = frateOdorList[glomnum][odoridx]
        plot(firingtsteps, frate, ['r','g','b'][glomnum], linewidth=2, marker=',')
    axes_labels(ax,'time (s)','ORN firing rate (Hz)',fontsize=30)
    ## odor B
    fig = figure(facecolor='w')
    ax = fig.add_subplot(111)
    title('Odor B: glom0 (r); glom1 (g); glom2 (b)', fontsize=30)
    for glomnum in [0,1,2]:
        frateperglomList = []
        odoridx = 0
        frate = frateOdorList[glomnum][odoridx]
        plot(firingtsteps, frate, ['r','g','b'][glomnum], linewidth=2, marker=',')
    axes_labels(ax,'time (s)','ORN firing rate (Hz)',fontsize=30)
    # air
    fig = figure(facecolor='w')
    ax = fig.add_subplot(111)
    title('air: glom0 (r); glom1 (g); glom2 (b)', fontsize=30)
    for glomnum in [0,1,2]:
        frateperglomList = []
        odoridx = 6
        frate = frateOdorList[glomnum][odoridx]
        plot(firingtsteps, frate, ['r','g','b'][glomnum], linewidth=2, marker=',')
    axes_labels(ax,'time (s)','ORN firing rate (Hz)',fontsize=30)

def pulse_stimuli(numgloms=NUM_GLOMS):
    pulsetime = arange(0.0,PULSE_DURATION+1e-10,FIRINGFILLDT)
    for glomnum in range(numgloms):
        fig1 = figure()
        title('Glomerulus '+str(glomnum)+' : Pulse Responses')
        ax1 = fig1.add_subplot(111)
        lenpulselist = len(pulseList)
        for pulseidx,(pulsedelayA,pulsedurationA,pulsedelayB,pulsedurationB) in enumerate(pulseList):
            frate = fratePulseList[glomnum][0][pulseidx]
            colorfrac = pulseidx/float(lenpulselist)
            pulseA = zeros([len(pulsetime)])+pulseidx*3
            pulseA[int(pulsedelayA/FIRINGFILLDT):int((pulsedelayA+pulsedurationA)/FIRINGFILLDT)]=pulseidx*3+2
            ax1.plot(pulsetime, pulseA, color=(colorfrac,1-colorfrac,1), marker=',')
            ax1.plot(pulsetime, frate, color=(colorfrac,1-colorfrac,0), marker=',')

def kernels(numgloms=NUM_GLOMS):
    for glomnum in range(numgloms):
        #### plot the kernels
        fig = figure(facecolor='w')
        ## A super axes to set a common ylabel
        bigAxes = axes(frameon=False) # hide frame
        xticks([]) # don't want to see any ticks on this axis
        yticks([])
        text(-0.14,0.3,'arb units', fontsize=24, rotation='vertical')
        
        kernelA = fratePulseList[glomnum][1]
        kAmin = min(kernelA)-0.01
        kAmax = max(kernelA)+0.01
        kernelB = fratePulseList[glomnum][2]
        kBmin = min(kernelB)-0.01
        kBmax = max(kernelB)+0.01
        kernelR = fratePulseList[glomnum][3]
        kRmin = min(kernelR)-0.01
        kRmax = max(kernelR)+0.01

        ax = fig.add_subplot(3,1,1)
        ax.set_yticks([kAmin,0,kAmax])
        ax.set_yticklabels(["%.2f"%kAmin,'0',"%.2f"%kAmax])
        for label in ax.get_yticklabels():
            label.set_fontsize(20)
        axes_off(ax,x=True,y=False)
        ax.plot(extratime, kernelA, color=(1,0,0), marker=',',\
            linestyle='solid',linewidth=2,label='OdorA kernel')
        ax.set_xlim(0,2)
        ax.set_ylim(kAmin,kAmax)
        biglegend()

        ax = fig.add_subplot(3,1,2)
        ax.set_yticks([kBmin,0,kBmax])
        ax.set_yticklabels(["%.2f"%kBmin,'0',"%.2f"%kBmax])
        for label in ax.get_yticklabels():
            label.set_fontsize(20)
        axes_off(ax,x=True,y=False)
        ax.plot(extratime, kernelB, color=(0,1,0), marker=',',\
            linestyle='solid',linewidth=2,label='OdorB kernel')
        ax.set_xlim(0,2)
        ax.set_ylim(kBmin,kBmax)
        biglegend()

        ax = fig.add_subplot(3,1,3)
        ax.set_yticks([kRmin,0,kRmax])
        ax.set_yticklabels(["%.2f"%kRmin,'0',"%.2f"%kRmax])
        for label in ax.get_yticklabels():
            label.set_fontsize(20)
        ax.set_xticks([0,2])
        ax.set_xticklabels(['0','2'])    
        ax.plot(extratime, kernelR, color=(0,0,0), marker=',',\
            linestyle='solid',linewidth=2,label='Air kernel')
        axes_labels(ax,'time (s)','',adjustpos=False)
        ax.set_xlim(0,2)
        ax.set_ylim(kRmin,kRmax)
        biglegend()
        fig.suptitle('Glomerulus '+str(glomnum)+': odor and air kernels',fontsize=24)

def random_pulse_stimuli(numgloms=NUM_GLOMS):
    numpulses = float(len(randomPulseList))
    pulse_air = randomPulseList[0]
    for glomnum in range(numgloms):
        frate_air = randomResponseList[glomnum][0]
        # the first 'frate' in randomResponseList[glomnum] is for air (rectangle pulse),
        # the second is a random pulsed air pulse.
        # then alternately odorA and odorB.
        # the last odorA and odorB frates are combined and given simultaneously.
        for pulse_i,frate in enumerate(randomResponseList[glomnum][1:-1]):
            # randomPulseList[glomnum][frate,...]
            if pulse_i == 0:
                # incompatible units addition! only schematic:
                # no rectangular pulse added, already pure air.
                pulse = array(randomPulseStepsList[1])*10.0 + ORN_BGND_MEAN
                frate = frate + ORN_BGND_MEAN
            elif pulse_i < (numpulses-3):
                # incompatible units addition! only schematic:
                pulse = array(randomPulseStepsList[pulse_i+1])*10.0 + pulse_air + ORN_BGND_MEAN
                frate = frate + frate_air + ORN_BGND_MEAN
            else:
                # incompatible units addition! only schematic:
                pulse = array(randomPulseStepsList[pulse_i+1])*10.0 + \
                        array(randomPulseStepsList[pulse_i+2])*10.0 + pulse_air + ORN_BGND_MEAN
                frate = frate + randomResponseList[glomnum][pulse_i+2] + frate_air + ORN_BGND_MEAN
            if pulse_i in [1,2,5]:
                fig1 = figure()
                nump = [1,2,5].index(pulse_i)
                title('Glomerulus '+str(glomnum)+' : Random Pulse Responses, odor '\
                    +['A','B','A+B'][nump])
                ax1 = fig1.add_subplot(111)
                plot(randompulsetime, pulse, color=['r','b','g'][nump],\
                    marker=',',linestyle='solid',label='pulse')
                #twinx()
                plot(randompulsetime, frate, color=['m','c','y'][nump],\
                    marker=',',linestyle='solid',label='ORN response')

def interglom_stimuli():
    # before taking mean activity over all glomeruli, we should clip negative values
    # because SA cells only receive _excitatory_ input from ET cells.
    frateOdorListClipped = clip(frateOdorList,0.0,1e6)
    # mean activity over gloms (axis=0) as a function of time
    fratemean = mean(frateOdorListClipped,axis=0)
    # Short axon network delays and neuronal integration is modelled as 25ms integrated excitation to PGs
    # the delay should be put in as part of the synapse
    integrate_window_length = int(SA_integrate_time/FIRINGFILLDT)
    # Make a normalized integration window
    # ones() makes array of floats by default, so no integer division problems
    rect = ones((integrate_window_length,))/integrate_window_length
    fig = figure()
    ax = fig.add_subplot(111)
    ax.set_title('SAs')
    ax.set_ylabel('Hz')
    ax.set_xlabel('time (s)')
    for odoridx,(odorA,odorB) in enumerate(inputList):
        ## convolution takes moving average,
        ## I must use mode='full' and take [0:numt], else the SA->(ET->)PG excitation is acausal.
        ## of course, the initial part from t = 0 till 'SA_integrate_time' is invalid,
        ## but that is within SETTLETIME.
        interglom_rate = convolve(fratemean[odoridx], rect, mode='full')[0:numt]
        #plot(firingtsteps, fratemean[odoridx], color=(odorA,odorB,0), marker=',')
        ax.plot(firingtsteps, interglom_rate, color=(odorA,odorB,0), marker=',')
        #### Extra plots of individual Odor responses
        #fig2 = figure()
        #ax2 = fig2.add_subplot(111)
        #for i in range(len(frateOdorListClipped)):
        #    ax2.plot(firingtsteps, frateOdorList[i][odoridx])

def paperORNinputfig(odornum):
    if odornum==0: # odor A
        kernelOdor = fratePulseList[central_glom][1]
        odorpulsenum = 2
    else:
        kernelOdor = fratePulseList[central_glom][2]
        odorpulsenum = 3
    kOmin = min(kernelOdor)-0.01
    kOmax = max(kernelOdor)+0.01
    kernelR = fratePulseList[central_glom][3]
    kRmin = min(kernelR)-0.01
    kRmax = max(kernelR)+0.01
    
    ######################## PULSE INPUT
    fig = figure(figsize=(columnwidth,linfig_height/2.0),dpi=300,facecolor='w') # 'none' is transparent

    ## odor kernel
    ax = plt.subplot2grid((2,8),(1,0),rowspan=1,colspan=2,frameon=False)
    #text(-0.1,0.9,'b', weight='bold', fontsize=label_fontsize, transform = ax.transAxes)
    text(0.3,0.75,'odor kernel', fontsize=label_fontsize, transform = ax.transAxes)
    ax.plot(extratime, kernelOdor, color='k', marker=',',\
        linestyle='solid',linewidth=linewidth,label='odor kernel')
    _,_,_,ymax = beautify_plot(ax,drawxaxis=False,drawyaxis=False,xticks=[],yticks=[])

    add_scalebar(ax,matchx=False,matchy=False,hidex=True,hidey=False,\
        sizex=0.5,labelx='0.5 s',sizey=20,labely='arb',\
        bbox_to_anchor=[1.2,0.0],bbox_transform=ax.transAxes)

    ## air kernel
    ax = plt.subplot2grid((2,8),(0,0),rowspan=1,colspan=2,frameon=False)
    #text(-0.1,0.9,'a', weight='bold', fontsize=label_fontsize, transform = ax.transAxes)
    text(0.35,0.75,'air kernel', fontsize=label_fontsize, transform = ax.transAxes)
    ax.plot(extratime, kernelR, color='k', marker=',',\
        linestyle='solid',linewidth=linewidth,label='air kernel')
    ## xtickxposn='none' not working, hence set_xticks() below!
    beautify_plot(ax,xticksposn='none',drawxaxis=False,drawyaxis=False,xticks=[],yticks=[])
    ## set the air ymax to that of the odor, so that they can be compared
    ax.set_ylim([0,ymax])

    pulse_air = randomPulseList[0]
    ## the first 'frate' in randomResponseList[glomnum] is for air (rectangle pulse),
    ## the second is a random pulsed air pulse.
    ## then alternately odorA and odorB.
    ## the last odorA and odorB frates are combined and given simultaneously.
    frate_air = randomResponseList[central_glom][0]
    ## randomPulseList[glomnum][frate,...]
    pulse = array(randomPulseStepsList[odorpulsenum])
    frate = randomResponseList[central_glom][odorpulsenum] + frate_air + ORN_BGND_MEAN
    ## suction pulse
    ax = plt.subplot2grid((2,8),(0,3),rowspan=1,colspan=2,frameon=False)
    #text(-0.1,0.9,'c', weight='bold', fontsize=label_fontsize, transform = ax.transAxes)
    plot(randompulsetime, pulse_air, color='k',linewidth=linewidth,\
        marker=',',linestyle='solid',label='suction')
    beautify_plot(ax,drawxaxis=False,drawyaxis=False,xticks=[],yticks=[])
    add_scalebar(ax,matchx=False,matchy=False,hidex=True,hidey=False,\
        sizex=2,labelx='2 s',sizey=0.1,labely='0.1',\
        bbox_to_anchor=[0.7,-0.1],bbox_transform=ax.transAxes)

    ## random pulse
    ax = plt.subplot2grid((2,8),(1,3),rowspan=1,colspan=2,frameon=False)
    #text(-0.1,0.9,'d', weight='bold', fontsize=label_fontsize, transform = ax.transAxes)
    plot(randompulsetime, pulse, color='k',linewidth=plot_linewidth,\
        marker=',',linestyle='solid',label='pulse')
    beautify_plot(ax,yticksposn='none',drawxaxis=False,drawyaxis=False,xticks=[],yticks=[])
    ## ORN frate
    ax = plt.subplot2grid((2,8),(0,6),rowspan=2,colspan=2,frameon=False)
    #text(-0.1,0.9,'e', weight='bold', fontsize=label_fontsize, transform = ax.transAxes)
    plot(randompulsetime, frate, color='k',linewidth=linewidth,\
        marker=',',linestyle='solid',label='ORN response')
    beautify_plot(ax,drawxaxis=False,drawyaxis=False,xticks=[],yticks=[])
    ax.set_ylim([0,10]) # IMP: needed to have same scalebar as resp response below!
    add_scalebar(ax,matchx=False,matchy=False,hidex=True,hidey=False,\
        sizex=2,labelx='2 s',sizey=2,labely='2 Hz',\
        bbox_to_anchor=[1.0,-0.05],bbox_transform=ax.transAxes)
    
    ## pulse input pictorial math
    ax.arrow( 0.3, 0.7, 0.06, 0, fc="k", ec="k", head_width=0.01, head_length=0.02,\
        linewidth=linewidth, transform=fig.transFigure )
    ax.arrow( 0.3, 0.3, 0.06, 0, fc="k", ec="k", head_width=0.01, head_length=0.02,\
        linewidth=linewidth, transform=fig.transFigure )
    ax.text( 0.32, 0.7, '*', fontsize=label_fontsize*1.2, transform=fig.transFigure )
    ax.text( 0.32, 0.3, '*', fontsize=label_fontsize*1.2, transform=fig.transFigure )
    ax.arrow( 0.62, 0.7, 0.08, -0.15, fc="k", ec="k", head_width=0.01, head_length=0.02,\
        linewidth=linewidth, transform=fig.transFigure )
    ax.arrow( 0.62, 0.3, 0.08, 0.15, fc="k", ec="k", head_width=0.01, head_length=0.02,\
        linewidth=linewidth, transform=fig.transFigure )
    ax.text( 0.722, 0.47, '+', fontsize=label_fontsize, transform=fig.transFigure )
    ax.arrow( 0.75, 0.5, 0.03, 0, fc="k", ec="k", head_width=0.01, head_length=0.02,\
        linewidth=linewidth, transform=fig.transFigure )

    fig.tight_layout()
    #fig.canvas.draw() ## redraws the figure
    subplots_adjust(left=0.03)
    fig_clip_off(fig)
    fig.savefig('../figures/input_example_pulse.svg',dpi=fig.dpi)
    fig.savefig('../figures/input_example_pulse.png',dpi=fig.dpi)

    #################### RESPIRATORY INPUT
    fig = figure(figsize=(columnwidth,linfig_height/2.0),dpi=300,facecolor='w') # 'none' is transparent

    ## odor kernel
    ax = plt.subplot2grid((2,8),(1,0),rowspan=1,colspan=2,frameon=False)
    #text(-0.1,0.9,'g', weight='bold', fontsize=label_fontsize, transform = ax.transAxes)
    text(0.3,0.75,'odor kernel', fontsize=label_fontsize, transform = ax.transAxes)
    ax.plot(extratime, kernelOdor, color='k', marker=',',\
        linestyle='solid',linewidth=linewidth,label='odor kernel')
    _,_,_,ymax = beautify_plot(ax,drawxaxis=False,drawyaxis=False,xticks=[],yticks=[])

    add_scalebar(ax,matchx=False,matchy=False,hidex=True,hidey=False,\
        sizex=0.5,labelx='0.5 s',sizey=20,labely='arb',\
        bbox_to_anchor=[1.2,0.0],bbox_transform=ax.transAxes)

    ## air kernel
    ax = plt.subplot2grid((2,8),(0,0),rowspan=1,colspan=2,frameon=False)
    #text(-0.1,0.9,'f', weight='bold', fontsize=label_fontsize, transform = ax.transAxes)
    text(0.35,0.75,'air kernel', fontsize=label_fontsize, transform = ax.transAxes)
    ax.plot(extratime, kernelR, color='k', marker=',',\
        linestyle='solid',linewidth=linewidth,label='air kernel')
    ## xtickxposn='none' not working, hence set_xticks() below!
    beautify_plot(ax,xticksposn='none',drawxaxis=False,drawyaxis=False,xticks=[],yticks=[])
    ## set the air ymax to that of the odor, so that they can be compared
    ax.set_ylim([0,ymax])

    ## resp waveform
    for resprow in [0,1]: ## plot two resp waveforms
        ax = plt.subplot2grid((2,8),(resprow,3),rowspan=1,colspan=2,frameon=False)
        #text(-0.1,0.9,['h','i'][resprow], weight='bold', fontsize=label_fontsize, transform = ax.transAxes)
        unclippedRespirationPulses = \
            odor_gen.normresponse(array(\
                [odor_gen.periodic_respiration_fn(t,negclip=False) for t in odor_gen.extratime]))
        clippedRespirationPulses = \
            odor_gen.normresponse(array(\
                [odor_gen.periodic_respiration_fn(t,negclip=True) for t in odor_gen.extratime]))
        #plot(odor_gen.extratime[odor_gen.startidx:], unclippedRespirationPulses[odor_gen.startidx:],\
        #    color='k',linewidth=linewidth,marker='',linestyle='dotted',label='respiration cycles')
        plot(odor_gen.extratime[odor_gen.startidx:], clippedRespirationPulses[odor_gen.startidx:],\
            color='k',linewidth=linewidth,marker=',',linestyle='solid',label='respiration cycles')
        beautify_plot(ax,drawxaxis=False,drawyaxis=False,xticks=[],yticks=[])    
    add_scalebar(ax,matchx=False,matchy=False,hidex=True,hidey=False,\
        sizey=0.3,labely='0.3',sizex=0.5,labelx='0.5 s',\
        bbox_to_anchor=[1.5,-0.1],bbox_transform=ax.transAxes)

    ## ORN resp frate
    ax = plt.subplot2grid((2,8),(0,6),rowspan=2,colspan=2,frameon=False)
    #text(-0.1,0.9,'j', weight='bold', fontsize=label_fontsize, transform = ax.transAxes)
    frate = frateOdorList[central_glom][odornum]
    plot(firingtsteps, frate, 'k', linewidth=linewidth, marker=',')
    beautify_plot(ax,drawxaxis=False,drawyaxis=False,xticks=[],yticks=[])
    ax.set_ylim([0,10]) # IMP: needed to have same scalebar as traech-ed response below!
    add_scalebar(ax,matchx=False,matchy=False,hidex=True,hidey=False,\
        sizex=0.5,labelx='0.5 s',sizey=2,labely='2 Hz',\
        bbox_to_anchor=[1.0,-0.05],bbox_transform=ax.transAxes)

    ## resp input pictorial math
    ax.arrow( 0.3, 0.7, 0.06, 0, fc="k", ec="k", head_width=0.01, head_length=0.02,\
        linewidth=linewidth, transform=fig.transFigure )
    ax.arrow( 0.3, 0.3, 0.06, 0, fc="k", ec="k", head_width=0.01, head_length=0.02,\
        linewidth=linewidth, transform=fig.transFigure )
    ax.text( 0.32, 0.7, '*', fontsize=label_fontsize*1.2, transform=fig.transFigure )
    ax.text( 0.32, 0.3, '*', fontsize=label_fontsize*1.2, transform=fig.transFigure )
    ax.arrow( 0.62, 0.7, 0.08, -0.15, fc="k", ec="k", head_width=0.01, head_length=0.02,\
        linewidth=linewidth, transform=fig.transFigure )
    ax.arrow( 0.62, 0.3, 0.08, 0.15, fc="k", ec="k", head_width=0.01, head_length=0.02,\
        linewidth=linewidth, transform=fig.transFigure )
    ax.text( 0.722, 0.47, '+', fontsize=label_fontsize, transform=fig.transFigure )
    ax.arrow( 0.75, 0.5, 0.03, 0, fc="k", ec="k", head_width=0.01, head_length=0.02,\
        linewidth=linewidth, transform=fig.transFigure )

    fig.tight_layout()
    #fig.canvas.draw() ## redraws the figure
    subplots_adjust(left=0.03)
    fig_clip_off(fig)
    fig.savefig('../figures/input_example_resp.svg',dpi=fig.dpi)
    fig.savefig('../figures/input_example_resp.png',dpi=fig.dpi)

def paper_decorr_schematic(odoridx):
    """ 5 for odor A, 0 for odor B.
    Take only the last respiration cycle.
    """
    lastrespstart = -int(RESPIRATION/FIRINGFILLDT)
    for glomnum in [0,1,2]:
        fig = figure(figsize=(columnwidth/5.0,linfig_height/4.0),\
            dpi=300,facecolor='none') # none is transparent
        ax = fig.add_subplot(111)
        frate = frateOdorList[glomnum][odoridx]
        plot(firingtsteps[lastrespstart:], frate[lastrespstart:],\
            ['r','g','b'][glomnum], linestyle=['solid','solid','dashed'][glomnum],\
            linewidth=linewidth)
        beautify_plot(ax,drawxaxis=False,drawyaxis=False,xticks=[],yticks=[])
        if glomnum==2:
            add_scalebar(ax,matchx=False,matchy=False,hidex=True,hidey=True,\
                sizex=0.2,labelx='0.2 s',sizey=2,labely='2 Hz',\
                bbox_to_anchor=[0.8,-0.2],bbox_transform=ax.transAxes)
        ax.set_xlim(RUNTIME-RESPIRATION,RUNTIME)
        ymin,ymax=ax.get_ylim()
        ax.set_ylim(0,10.5) # 10.5 Hz for all plots, to get common scale-bar
        fig_clip_off(fig)
        fig.tight_layout()
        fig.savefig('../figures/decorr/input_glom'+str(glomnum)+'.svg',\
            dpi=fig.dpi,transparent=True)

if __name__ == "__main__":
    #interglom_stimuli()
    #odor_and_air_stimuli()
    #pulse_stimuli()
    #random_pulse_stimuli(3)
    #kernels(3)
    #paperORNinputfig(0) ## 0 for odor A, 1 for odor B
    paper_decorr_schematic(0) ## 5 for odor A, 0 for odor B

    show()
