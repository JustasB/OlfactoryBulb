# -*- coding: utf-8 -*-
"""
Created on Sun Mar  6 18:22:04 2011

@author: -
"""
import os
import numpy
from matplotlib import pyplot
from neuronpy.graphics import spikeplot
from bulbspikes import *
from neuronpy.util import spiketrain
from params import sim_var

homedir = os.path.join(os.path.relpath('..'))
analysis_path = homedir

def format_axes(ax, dt=1, ylim=(0.,4.)):
    #ax.set_xticks(numpy.arange(0,num_intervals,(num_intervals-1)/4.))
    #ax.set_xticklabels(['$-\pi$','$-\pi/2$','$0$','$\pi/2$','$\pi$'], fontsize=18)
    xlim = ax.get_xlim()
    timesteps=int((xlim[1]*dt-xlim[0]*dt)/2.)
    ax.set_xticks(numpy.linspace(xlim[0],xlim[1],5))
    ax.set_xticklabels(numpy.asarray(numpy.linspace(-timesteps,timesteps,5), dtype=int))
    ax.set_xlabel('lag (ms)')
    ax.set_ylim(ylim)
    ax.set_ylabel('Synchronization magnitude')

def draw_cell(cellid, ax, color='black'):
    xloc = 10+cellid*20
    # Lateral dends
    y = numpy.abs(numpy.subtract(range(101), xloc))
    yvec = numpy.log(numpy.add(y,1))
    ax.plot(range(101), yvec, color=color)
    # Soma
    ax.fill_between(range(101), numpy.ones(101), yvec, \
        where=numpy.ma.masked_where(yvec < 1., yvec).mask, \
        color=color, linewidth=0.)
    # Glom
    ax.plot([xloc], [9], color=color, marker='o', markersize=10, markerfacecolor='white', markeredgecolor=color)
    ax.plot([xloc], [9], color=color, marker='o', markersize=9, alpha=0.25)
    ax.plot([xloc], [9], color=color, marker='1', markersize=7, markeredgewidth=2)
    # Primary dendrite
    ax.plot([xloc, xloc], [0,8], color=color, linewidth=2)
    format_schematic_axis(ax)
    
def draw_weights(cellids, ax, color='black',scale=1.):
    """Draw granule cells"""
    import synweightsnapshot
    sws = synweightsnapshot.SynWeightSnapshot( \
            nummit=sim_var['num_mitral'], \
            numgran=sim_var['num_granule'])
    
    raw=sws.read_file(sim_var['wt_input_file'], 
            os.path.join(homedir, sim_var['weight_dir']))
    sws.parse_data(raw)
    for cellid in cellids:
        wts = sws.m2g[cellid,:,0]
        wts = wts/numpy.max(wts)
    
        for i in range(len(wts)):
            if wts[i] > 0.0001:
                cellloc = 10+cellid*20
                y = numpy.abs(i - cellloc)
                yloc = numpy.log(numpy.add(y,1))
                gloc = -3.5+((i%2)*1.5)
                ax.plot([i],[yloc], marker='o', markerfacecolor=color, markersize=4.*scale, markeredgecolor=color)
                ax.plot([i,i],[yloc, gloc], color=color)
                ax.plot([i],[gloc], marker='^', markerfacecolor=color, markersize=6.*scale, markeredgecolor=color)
    format_schematic_axis(ax)
    
def format_schematic_axis(ax):
    ax.set_xlim((0,100))
    xticks = [10,30,50,70,90]
    ax.set_xticks(xticks)
    ax.set_xticklabels(numpy.multiply(xticks,10))
    ax.set_xlabel('distance in microns')
    ax.set_ylim((-5,11))
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.set_yticks([])
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('black')
    ax.xaxis.set_ticks_position('bottom')


def read_weightevents():
    M = numpy.loadtxt(os.path.join(analysis_path, 'stimweightevents.txt'))
    data = []
    for i in range(5):
        data.append([])
    for m in M:
        data[int(m[0])].append(m[1])
    return data

def read_delayevents():
    M = numpy.loadtxt(os.path.join(analysis_path, 'stimdelayevents.txt'))
    data = []
    for i in range(5):
        data.append([])
    for m in M:
        data[int(m[0])].append(m[1])
    return data


def raster(pair=[0,4], cluster_width=5, fi=.005, xlim=(1000,2000)):
#    pos1 = (10+pair[0]*20, cluster_width, 1, pair)
#    pos2 = (10+pair[1]*20, cluster_width, 1, pair)
#    stim_odor_mags = numpy.ones(5)*.55

    fig = pyplot.figure(figsize=(9.5,5.7))
    raster_ax = fig.add_axes([.1,.1,.8,.27])
    schematic_ax = fig.add_axes([.1,.85,.8,.1])
    syn_ax = fig.add_axes([.1,.45,.8,.225])

    draw_cell(pair[0], schematic_ax, color='red')
    draw_cell(pair[1], schematic_ax, color='blue')
    draw_weights(pair, schematic_ax, color='black')

    # Analyze an output file in some_dir
    bulb_spikes = BulbSpikes(sim_time=sim_var['tstop'])
    bulb_spikes.read_file(os.path.join(homedir,'spikeout.spk'))
    breath_events = numpy.loadtxt(os.path.join(homedir, 'breathevents.txt'))

    wts = read_weightevents()
    delays = read_delayevents()
    
    dt = 1
    tstop = xlim[1]
    x = numpy.arange(0,tstop,dt)
    y0 = numpy.zeros(tstop/dt)
    y1 = numpy.zeros(tstop/dt)
    EXP = numpy.exp(numpy.multiply(x,-1./200.))-numpy.exp( \
            numpy.multiply(x,-1./20.))
    
    idx = 0
    for b in breath_events:
        if b >= tstop:
            break
        else:
            dtidx = int((b+delays[pair[0]][idx])/dt)
            y0[dtidx:] += EXP[:-dtidx]*wts[pair[0]][idx]
            dtidx = int((b+delays[pair[1]][idx])/dt)
            y1[dtidx:] += EXP[:-dtidx]*wts[pair[1]][idx]
        idx += 1
    redplt = syn_ax.plot(x,y0, color='red')
    blueplt = syn_ax.plot(x,y1, color='blue')
    for breath in breath_events:
        breathplt = syn_ax.plot([breath, breath], [0,2], linestyle='--', \
                color='gray', linewidth=2)
    syn_ax.set_xlim(xlim)
    syn_ax.set_ylim(0,1.6)
    syn_ax.set_yticks([])
    syn_ax.set_xticks([])
    syn_ax.set_ylabel('EPSC onto tuft')
    leg = syn_ax.legend([breathplt, redplt, blueplt], \
            ['sniff event', 'input onto red', 'input onto blue'], \
            bbox_to_anchor=(0, 1.15, 1., .102), loc=1, ncol=3, mode="expand", \
            borderaxespad=0., handletextpad=.2)
    # Mark sniff interval
    for i in range(len(breath_events)):
        if breath_events[i] > xlim[0]:
            span = syn_ax.annotate('', xy=(breath_events[i], .28),  xycoords='data',
                               xytext=(breath_events[i+1], .28), \
                               textcoords='data', \
                               arrowprops=dict(arrowstyle="|-|", linewidth=2)
                               )
    
            syn_ax.text((breath_events[i]+breath_events[i+1])/2., .53, \
                    'sniff every\n150 - 250 ms', \
                    horizontalalignment='center', verticalalignment='top', \
                    backgroundcolor='white')
            break

    # Mark amplitude interval
    span = syn_ax.annotate('', xy=(1190, 1.28),  xycoords='data',
                       xytext=(1190, 1.12), \
                       textcoords='data', \
                       arrowprops=dict(arrowstyle="|-|", linewidth=2)
                       )

    syn_ax.text(1215, 1.21, \
            '+/- 5%', \
            horizontalalignment='left', verticalalignment='center')
    
    # Mark delay interval
    for i in range(len(breath_events)):
        if breath_events[i] > 1400:
            span = syn_ax.annotate('', xy=(breath_events[i]-2, .5),  xycoords='data',
                               xytext=(breath_events[i]+17, .5), \
                               textcoords='data', \
                               arrowprops=dict(arrowstyle="|-|", linewidth=2)
                               )
    
            syn_ax.text(breath_events[i]+7.5, .28, \
                    'delay 0-15 ms', \
                    horizontalalignment='center', verticalalignment='top', \
                    backgroundcolor='white')
            break
    

    spikes = bulb_spikes.get_mitral_spikes()
    ref=spikes[pair[0]]
    comp=spikes[pair[1]]
    gcspikes = bulb_spikes.get_granule_spikes()
    mididx = 10+pair[0]*20
    gcleft = gcspikes[mididx-int(cluster_width/2.):mididx+int(cluster_width/2.)+1]
    mididx = 10+pair[1]*20
    gcright = gcspikes[mididx-int(cluster_width/2.):mididx+int(cluster_width/2.)+1]

    sp = spikeplot.SpikePlot(fig=fig, savefig=False)
    sp.set_markercolor('blue')
    sp.set_markeredgewidth(2.)
    sp.set_markerscale(4)
    sp.plot_spikes([comp], label='comp', cell_offset=cluster_width*2+5, \
            draw=False )
    sp.set_markercolor('red')
    sp.plot_spikes([ref], label='ref', cell_offset=cluster_width*2+2, \
            draw=False)
    sp.set_markerscale(1.3)

    sp.set_markeredgewidth(1.5)
    sp.set_markercolor('blue')
    sp.plot_spikes(gcright, label='gcright', cell_offset=cluster_width, \
            draw=False)
    sp.set_markercolor('red')
    sp.plot_spikes(gcleft, label='gcleft', cell_offset=0, \
            draw=False)

    coincidences, mask_a, mask_b, ratio = \
            spiketrain.get_sync_traits(ref, comp, window=5)
#        idx = 0
#        for i in mask_a:
#            if i == 1:
#                raster_ax.plot([ref[idx]],[cluster_width*2+1.9], marker='o', color='red')
#            idx += 1
    idx = 0
    for i in mask_b:
        if i == 1:
            if comp[idx] >= xlim[0] and comp[idx] < xlim[1]:
                raster_ax.text(comp[idx],cluster_width*2+8.5, '*', \
                    color='purple', fontweight='bold', \
                    horizontalalignment='center', verticalalignment='center')
            #raster_ax.plot([comp[idx]],[cluster_width*2+7], marker='o', color='blue')
        idx += 1

    raster_ax.text(2000,cluster_width*2+8.5, '(synchronized)', color='purple', \
            horizontalalignment='center', verticalalignment='center',
            fontsize=11)

    raster_ax.set_yticks([])
    ylim = (0.5, cluster_width*2+7.5)
    for breath in breath_events:
        raster_ax.plot([breath, breath], [ylim[0], ylim[1]], linestyle='--', color='gray', linewidth=2)

    sp.update_xlim(xlim)
    raster_ax.set_ylim(ylim)
    raster_ax.set_xlabel('time (ms)')
    raster_ax.set_ylabel('spike output\n  granule      mitral\n\n', horizontalalignment='center')

    pos = schematic_ax.get_position()
    schematic_ax.text(.025, pos.ymax+.02, 'A)', transform=fig.transFigure, 
          verticalalignment='baseline')
    pos = syn_ax.get_position()
    syn_ax.text(.025, pos.ymax+.07, 'B)', transform=fig.transFigure, 
          verticalalignment='baseline')            
    pos = raster_ax.get_position()
    raster_ax.text(.025, pos.ymax+.02, 'C)', transform=fig.transFigure, 
          verticalalignment='baseline')            

#    fig.savefig(os.path.join(analysis_path, 'raster_w%d_(%d-%d)_%.3f.pdf') %(cluster_width, pair[0], pair[1], fi))
    fig.savefig(os.path.join(analysis_path, 'fig1.pdf'))
raster()
