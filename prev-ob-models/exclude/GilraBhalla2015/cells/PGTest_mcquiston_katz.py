#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import math

#### python2.6 PGTest_mcquiston_katz.py <PG2010|PG2013>

sys.path.extend(["..","../channels","neuroml"])

from moose.utils import *
from load_channels import *
from moose.neuroml.MorphML import *
from data_utils import *

from pylab import * # part of matplotlib that depends on numpy but not scipy

SIMDT = 1e-5 # seconds
PLOTDT = 1e-5 # seconds

class PGTest:

    def __init__(self,figname):
        load_channels()
        MML = MorphML({'temperature':CELSIUS})
        ## Figure 2G,H has 2013 (LTS) first, then 2010 (plateauing)
        ## 50pA injected if PG2010, else 100pA injected -- see in run() below
        if '2010' in figname:
            MML.readMorphMLFromFile('PG_aditya2010unified_neuroML_L1_L2_L3.xml',{})
            self.figname = 'PG2010'
        elif '2013' in figname:
            MML.readMorphMLFromFile('PG_aditya2013unified_neuroML_L1_L2_L3.xml',{})
            self.figname = 'PG2013'
        else:
            print "Give PG2010 or PG2013 as argument"
            sys.exit(1)
        self.libname = 'PG'
        #MML.readMorphMLFromFile('PG_aditya2013_neuroML_L1_L2_L3.xml',{})
        #self.libname = 'PG_LTS'
        self.context = moose.PyMooseBase.getContext()
        self.PGcell = self.context.deepCopy(self.context.pathToId('/library/'+self.libname),self.context.pathToId('/'),"PG")
        self.PGsoma = moose.Compartment('/PG/soma_0')
        self.soma_vm = self.setupTable('soma_vm', self.PGsoma,'Vm')
        self.caconc_conc = self.setupChannelTable('Ca_mit_conc_Ca', moose.CaConc(self.PGsoma.path+'/Ca_mit_conc'),'Ca')
        self.tca_Ik = self.setupChannelTable('TCa_Ik', moose.HHChannel(self.PGsoma.path+'/TCa_d'),'Ik')
        self.na_rat_ms_X = self.setupChannelTable('Na_rat_ms_X',moose.HHChannel(self.PGsoma.path+"/Na_rat_ms"),"X")
        self.na_rat_ms_Y = self.setupChannelTable('Na_rat_ms_Y',moose.HHChannel(self.PGsoma.path+"/Na_rat_ms"),"Y")
        self.kdr_ms_X = self.setupChannelTable('KDR_ms_X',moose.HHChannel(self.PGsoma.path+"/KDR_ms"),"X")
        self.ka_ms_X = self.setupChannelTable('KA_ms_X',moose.HHChannel(self.PGsoma.path+"/KA_ms"),"X")
        self.tca_d_X = self.setupChannelTable('TCa_d_X',moose.HHChannel(self.PGsoma.path+"/TCa_d"),"X")
        self.tca_d_Y = self.setupChannelTable('TCa_d_Y',moose.HHChannel(self.PGsoma.path+"/TCa_d"),"Y")
        self.ih_cb_X = self.setupChannelTable('Ih_cb_X',moose.HHChannel(self.PGsoma.path+"/Ih_cb"),"X")

    def run(self):
        self.context.setClock(0, SIMDT, 0)
        self.context.setClock(1, SIMDT, 0) #### The hsolve and ee methods use clock 1
        self.context.setClock(2, SIMDT, 0) #### hsolve uses clock 2 for mg_block, nmdachan and others.
        self.context.setClock(PLOTCLOCK, PLOTDT, 0)
        self.context.reset()
        
        ########### Settle
        self.PGsoma.inject = 0.0
        self.context.step(0.5) # must be float to interpret as runtime, integer is interpreted as number of steps
        oldlen = len(pg.soma_vm)
        self.inject = ones(oldlen)*self.PGsoma.inject

        ########### Hyperpolarize
        ## I am using -50pA, ideally should use -100pA as per McQuiston and Katz.
        ## But then, the hyperpolarization is -120mV for PG2010.
        ## Figure 6 has 2013 first, then 2010
        if '2010' in self.figname:
            self.PGsoma.inject = -50e-12 # 50pA for 2010 PG
        else:
            self.PGsoma.inject = -100e-12 # 100pA for 2013 PG
        self.context.step(0.6) # must be float to interpret as runtime, integer is interpreted as number of steps
        self.inject = append(self.inject,ones(len(pg.soma_vm)-oldlen)*self.PGsoma.inject)
        oldlen = len(pg.soma_vm)

        ########### No injection
        self.PGsoma.inject = 0.0
        self.context.step(0.5) # must be float to interpret as runtime, integer is interpreted as number of steps
        self.inject = append(self.inject,ones(len(pg.soma_vm)-oldlen)*self.PGsoma.inject)
        oldlen = len(pg.soma_vm)

        ########### Depolarize
        ## I am using 50pA, ideally should use 100pA as per McQuiston and Katz.
        ## But then have to put in too much K2 and KA to get pegging during depolzn,
        ## and then cell does not spike even with 10000Hz synaptic input.
        if '2010' in self.figname:
            self.PGsoma.inject = 50e-12 # 50pA for 2010 PG
        else:
            self.PGsoma.inject = 100e-12 # 100pA for 2013 PG
        self.context.step(0.6) # must be float to interpret as runtime, integer is interpreted as number of steps
        self.inject = append(self.inject,ones(len(pg.soma_vm)-oldlen)*self.PGsoma.inject)
        oldlen = len(pg.soma_vm)

        ############ No injection
        self.PGsoma.inject = 0.0
        self.context.step(0.5) # must be float to interpret as runtime, integer is interpreted as number of steps
        self.inject = append(self.inject,ones(len(pg.soma_vm)-oldlen)*self.PGsoma.inject)
        oldlen = len(pg.soma_vm)

    def setupTable(self, name, compmt, qtyname):
        # Setup the tables to pull data
        vmTable = moose.Table(name, moose.Neutral(compmt.path+"/data"))
        vmTable.stepMode = TAB_BUF #TAB_BUF: table acts as a buffer.
        vmTable.connect("inputRequest", compmt, qtyname)
        vmTable.useClock(PLOTCLOCK)
        return vmTable

    def setupChannelTable(self, name, channel, qtyname):
        # Setup the tables to pull data
        vmTable = moose.Table(name, moose.Neutral(channel.path+"/data"))
        vmTable.stepMode = TAB_BUF #TAB_BUF: table acts as a buffer.
        vmTable.connect("inputRequest", channel, qtyname)
        vmTable.useClock(PLOTCLOCK)
        return vmTable
        
if __name__ == "__main__":
    pg = PGTest(sys.argv[1])
    print "soma diameter = ",pg.PGsoma.diameter," m."
    print "soma length = ",pg.PGsoma.length," m."
    print "soma Rm = ",pg.PGsoma.Rm," Ohms."
    print "soma Cm = ",pg.PGsoma.Cm," Farads."
    print "soma Ra = ",pg.PGsoma.Ra," Ohms."
    print "soma Na gmax = ",moose.HHChannel(pg.PGsoma.path+"/Na_rat_ms").Gbar," Siemens."
    print "soma K gmax = ",moose.HHChannel(pg.PGsoma.path+"/KDR_ms").Gbar," Siemens."
    pg.run()
    ###### Paper figure 2: cells electrophysiology
    ## Figure 6 has 2013 first, then 2010
    tlist = arange(0.0,PLOTDT*len(pg.soma_vm),PLOTDT)*1000 - 450
    ## top figure is given less height as it doesn't need xticklabels and xlabel
    if pg.figname == 'PG2013': figheight = linfig_height/2.0 * 0.85
    else: figheight = linfig_height/2.0 * 1.15
    fig = figure(figsize=(columnwidth/2.0, figheight), dpi=300, facecolor='white')
    if pg.figname == 'PG2013':
        gs1 = GridSpec(4,1)
        gs1.update(left=0.28, right=0.95, top=0.9, bottom=0.05)
    else:
        gs1 = GridSpec(4,1)
        gs1.update(left=0.28, right=0.95, top=0.93, bottom=0.3)
        gs2 = GridSpec(1,1)
        ## top of gs2 should not equal bottom of gs1,
        ## else one plot at the bottom of gs2 disappears
        gs2.update(left=0.28, right=0.95, top=0.299, bottom=0.21)
    
    ## Volatge trace
    ax1 = plt.subplot(gs1[:-1,:]) # except bottom row of gs1
    ax1.plot(tlist,array(pg.soma_vm)*1e3,',-k',linewidth=plot_linewidth)

    ## current injection
    ax2 = plt.subplot(gs1[-1,:]) # bottom row of gs1
    ax2.plot(tlist,array(pg.inject)*1e12,',-k',linewidth=plot_linewidth)
    axes_labels(ax2,'','',fontsize=label_fontsize) # sets default label_fontsize

    for ax in [ax1,ax2]:
        xmin,xmax,ymin,ymax = beautify_plot(ax,x0min=True,y0min=False,\
                                drawxaxis=False,drawyaxis=True,xticks=[],yticks=[])
        ax.set_yticks([ymin,0,ymax])
        ax.set_xlim(0,2000)
        ax.set_xticks([])
        ax2.set_yticks([])

    if pg.figname == 'PG2013':
        axes_labels(ax1,'','Vm (mV)',fontsize=label_fontsize,xpad=1,ypad=-2) # sets default label_fontsize
        ax2.set_ylabel('200 pA',fontsize=label_fontsize,rotation='horizontal',labelpad=1)
    else:
        axes_labels(ax1,'','Vm (mV)',fontsize=label_fontsize,xpad=1,ypad=-6) # sets default label_fontsize        
        ax2.set_ylabel('100 pA',fontsize=label_fontsize,rotation='horizontal',labelpad=1)
        ax3 = plt.subplot(gs2[:,:]) # full gs2 for the bottom time axis
        ax3.set_xlim(0,2000)
        beautify_plot(ax3,x0min=True,y0min=False,\
            drawxaxis=True,drawyaxis=False,xticks=[],yticks=[])
        ax3.set_xticks([0,1000,2000])
        ax3.set_xticklabels(['0','1','2'])
        axes_labels(ax3,'time (s)','',fontsize=label_fontsize,xpad=1) # sets default label_fontsize

    fig_clip_off(fig)
    fig.savefig('../figures/connectivity/cells/'+pg.figname+'_ep.png', dpi=fig.dpi)
    ## Somehow this svg becomes 15MB and inkscape needs ~3GM RAM to load it!
    ## Basically this svg uses <use xlink:href=...> tags which other matplotlib svg-s don't! Surprising!
    #fig.savefig('../figures/connectivity/cells/'+pg.figname+'_ep.svg')

    #figure()
    #plot(pg.na_rat_ms_X,',-b')
    #plot(pg.na_rat_ms_Y,',-g')
    #plot(pg.kdr_ms_X,',-r')
    #plot(pg.ka_ms_X,',-m')
    #plot(pg.tca_d_X,',-y')
    #plot(pg.tca_d_Y,',-c')
    #plot(pg.ih_cb_X,',-k')
    #plot(pg.caconc_conc,',-g')
    #plot(pg.tca_Ik,',-g')

    show()
