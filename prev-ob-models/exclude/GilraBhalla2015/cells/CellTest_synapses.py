#!/usr/bin/env python
# -*- coding: utf-8 -*-

########## Need to run:
## ~/Python2.6.4/bin/python2.6 CellTest_synapses.py <cellname> <synapsename> \
##      <Iclamp/Vclamp> <numsyns_tobeplotted> <timeforeachsynapsesim> \
##      <simultaneous|staggered> <compartment_name>
## the option <simultaneous|stagerred> only needs to be given if numsyns_tobeplotted is > 1.
## synapsefactor would be say GABA_factor or AMPA_factor
## (used in the network -- see networkConstants.py), typically 1.0, etc.
## Presently only event based synapses are supported.
## Typical usage: python2.6 CellTest_synapses.py mitral ORN_mitral <Vclamp|Iclamp> 1 500e-3 simultaneous
## Eg1: python2.6 CellTest_synapses.py mitral granule_mitral <Vclamp|Iclamp> 1 500e-3 simultaneous Seg0_prim_dend_1_16
## Eg2: python2.6 CellTest_synapses.py mitral ORN_mitral Iclamp 1 500e-3 simultaneous Seg0_glom_81_102
## Eg3: python2.6 CellTest_synapses.py granule mitral_granule Iclamp 100 6e-3 staggered
## Eg4: python2.6 CellTest_synapses.py mitral granule_mitral Iclamp 50 100e-3 staggered

import os
import sys
import math

sys.path.extend(["..","../channels","../synapses","../neuroml","../simulations","../networks"])

from moose.utils import *
from synapseConstants import *
from simset_inhibition import * # has SIMDT and PLOTDT but SETTLETIME is overridden below
SETTLETIME = 500e-3 # give a large settletime for measuring post synaptic potentials

from load_channels import *
from load_synapses import *
from networkConstants import * # for GABA_factor, AMPA_factor and NMDA_factor

from pylab import * # part of matplotlib that depends on numpy but not scipy
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes # for inset
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

from data_utils import *

from moose.neuroml.MorphML import *

SATURATING_MITGRAN = False
if SATURATING_MITGRAN:
    mitgran_ampa = 'mitral_granule_saturatingAMPA'
    mitgran_nmda = 'mitral_granule_saturatingNMDA'
else:
    mitgran_ampa = 'mitral_granule_AMPA'
    mitgran_nmda = 'mitral_granule_NMDA'

#### inject a current into the tuft of mitral.
#### This to check effect of IPSP in prim dend vs sec dend.
#injtuft = 200e-12#1e-10 # A
injtuft = 0.0 # A
#injsoma = 400e-12 # A
injsoma = 0.0 # A

## nerve_shock is for injecting a sharp current into the mitral tuft after SETTLETIME
nerve_shock = False
num_tuft_comps_shock = 53 # number of tuft compartents to inject current into
tuft_shock_inject = 1e-12 # this current is injected into each of the num_tuft_comps_shock above

## set 0 mV for IPSCs, -65mV for EPSCs
VCLAMP_V = -65e-3#0e-3 # mV

class CellTest:

    def __init__(self, cellname, pseudosynname, libsynname,
        synfactor, Vclamp, numsyns, syntime, stagger, comp_name=None):
        load_channels()
        ## synchan_activation_correction is for graded synapses -- not useful for this script
        load_synapses(synchan_activation_correction)
        MML = MorphML({'temperature':CELSIUS})
        if cellname == 'PG':
            #filename = 'PG_aditya2010_neuroML_L1_L2_L3.xml'
            #filename = 'PG_aditya2012_neuroML_L1_L2_L3.xml'
            #filename = 'PG_aditya2013_neuroML_L1_L2_L3.xml'
            #cellname = 'PG_LTS'
            ## Choose one the two below
            ## set the ORN->PG and mitral->PG to 0.45 vs 1.25 nS for plateauing vs LTS
            #filename = 'PG_aditya2010unified_neuroML_L1_L2_L3.xml' # plateauing i.e. non-LTS
            filename = 'PG_aditya2013unified_neuroML_L1_L2_L3.xml' # LTS
            self.somaName = 'soma_0'
        elif cellname == 'mitral':
            #filename = 'mitral_bbmit1993davison_neuroML_L1_L2_L3_mod.xml'
            #filename = 'mitral_bbmit1993davison_neuroML_L1_L2_L3.xml'
            filename = 'mitral_bbmit1993davison_neuroML_L1_L2_L3_mod_withspikeinit.xml'
            #filename = 'mitral_bbmit1993davison_neuroML_TEST_L1_L2_L3.xml'
            #filename = 'mitral_bbmit1993davisonMS_neuroML_L1_L2_L3.xml'
            #filename = 'mitral_migliore_etal_2007.xml'
            self.somaName = 'Seg0_soma_0'
        elif cellname == 'granule':
            filename = 'granule_granadityaMS2007_neuroML_L1_L2_L3.xml'
            self.somaName = 'soma_0'
        self.cellName = cellname
        self.synName = libsynname
        self.pseudoSynName = pseudosynname
        self.synFactor = synfactor
        self.Vclamp = Vclamp
        self.numSyns = numsyns
        self.synTime = syntime
        self.stagger = stagger
        self.compName = comp_name
        if self.stagger:
            self.synRUNTIME = self.synTime*(self.numSyns)+2*SETTLETIME
        else:
            self.synRUNTIME = self.synTime+2*SETTLETIME

        self.cellsDict = MML.readMorphMLFromFile(filename,{})
        self.context = moose.PyMooseBase.getContext()
        self.cell = moose.Cell( self.context.deepCopy(\
            self.context.pathToId('/library/'+cellname),self.context.pathToId('/'),cellname) )
        self.soma = moose.Compartment('/'+self.cellName+'/'+self.somaName)
        self.somaVm = setupTable('soma_vm', self.soma,'Vm')
        if cellname=='mitral':
            ## inject a current into tuft of mitral
            tuftcomp = moose.Compartment('/'+self.cellName+'/Seg0_glom_81_102')
            tuftcomp.inject = injtuft
            ## inject a current into soma of mitral
            self.soma.inject = injsoma
            ## if nerve shock, inject a sharp current into tuft compartments
            if nerve_shock:
                ## We want a single sharp current injection at SETTLETIME, so cannot use soma.inject
                ## segment info - see MorphML_reader.py
                ## Only connect to 25 tuft compartments
                for seginfo in (self.cellsDict[self.cellName].values())[0:num_tuft_comps_shock]:
                    ## seginfo[5] is a list of potential synapses at this segment
                    if 'ORN_mitral' in seginfo[5]:
                        ## wrap this segment, seginfo[0] is segment name
                        tuftcomp = moose.Compartment('/'+self.cellName+'/'+seginfo[0])
                    ## compartment, name_extn, start_time, duration, current (all SI)
                    setup_iclamp(tuftcomp, '_nerveshock',\
                        SETTLETIME, 1e-3, tuft_shock_inject)

        self.spikeTableList = []
        self.attachSpikeTables()
        #printCellTree(self.cell)

        if cellname=='mitral':
            ## monitor Vm at the base of the tuft
            tuftbase_seg = moose.Compartment('/'+self.cellName+\
                #'/Seg0_tuftden_19_23') # for migliore and shepherd 2007 cell
                '/Seg0_prim_dend_5_20') # tuft base
            self.tuftBaseTable = setupTable('mitdendTable',tuftbase_seg,'Vm')

            ## monitor Vm at secondary dendrite
            dend_seg = moose.Compartment('/'+self.cellName+\
                #'/Seg0_sec_dendp4_0_254') # at 43.86 microns from soma
                #'/Seg0_sec_dendd3_0_204') # at 190.56 microns from soma
                '/Seg0_sec_dendd4_2_269') # at 1004.47 microns from soma
            self.secDendTable = setupTable('mitdendTable',dend_seg,'Vm')

            ## monitor Vm in the tuft compartment
            ## same as current injection above
            self.tuftTable = setupTable('mitTuftTable',tuftcomp,'Vm')

        ### testing - channels
        ## setup only one of the Tables below
        #self.var = setupTable('PG_mitral', moose.SynChan('/mitral/Seg0_glom_0_21/PG_mitral'), 'Gk')
        #self.var = setupTable('mitral_granule_NMDA', \
        #    moose.SynChan('/'+self.cellName+'/'+self.compName+'/mitral_granule_NMDA'), 'Gk')
        #self.var = setupTable('PG_mitral', moose.SpikeGen('/PG/soma_0/PG_mitral_spikegen'), 'state')
        #self.var = setupTable('soma_vm', moose.Compartment('/PG/soma_0'), 'Vm')
        #self.var = setupTable('mitral_glom', moose.Compartment('/mitral/Seg0_glom_0_21'), 'Vm')
        #self.var = setupTable('mitral_secdend', moose.Compartment('/mitral/Seg0_sec_dendd2_7_200'), 'Vm')
        #self.var = setupTable('mitral_secdend', moose.Compartment('/mitral/Seg0_sec_dendd2_0_165'), 'Vm')
        #self.var = setupTable('mitral_glom', moose.Compartment('/mitral/Seg0_glom_81_102'), 'Vm')
        #self.var = setupTable('soma_vm', moose.Compartment('/mitral/Seg0_soma_0'), 'Vm')
        #self.var = setupTable('soma_ca', moose.CaConc(self.soma.path+'/Ca_mit_conc'), 'Ca')
        
    def attachSpikeTables(self):
        ## If user did not supply a compartment_name, make a list of synapses
        ## connected to allowed compartments in the cell's xml file
        if self.compName is None:
            synapse_list = []
            synapse_list2 = []
            ## segment info - see MorphML_reader.py
            for seginfo in self.cellsDict[self.cellName].values():
                ## seginfo[5] is a list of potential synapses at this segment
                if self.pseudoSynName in seginfo[5]:
                    ## wrap this segment, seginfo[0] is segment name
                    compartment = moose.Compartment('/'+self.cellName+'/'+seginfo[0])
                    synapse = connectSynapse(self.context, compartment, self.synName, self.synFactor)
                    ## NMDA is also created if AMPA mit->gran is made.
                    if self.synName == mitgran_ampa:
                        synapse2 = connectSynapse(self.context, compartment, mitgran_nmda, NMDA_factor)
                        synapse_list2.append(synapse2)
                    synapse_list.append(synapse)
                    self.compName = seginfo[0] # the last synapse's comp becomes self.compName
        ## If user supplied a compartment_name,
        ## just have a synapse connected to the specified compartment
        else:
            #printCellTree(self.cell)
            comp_path = '/'+self.cellName+'/'+self.compName
            if self.context.exists(comp_path): compartment = moose.Compartment(comp_path)
            else:
                print "Error: non-existent compartment", comp_path
                sys.exit(1)
            synapse_list = [ connectSynapse(self.context, compartment, self.synName, self.synFactor) ]
            if self.synName == mitgran_ampa:
                synapse_list2 = [ connectSynapse(self.context, compartment, mitgran_nmda, NMDA_factor) ]
        for i in range(self.numSyns):
            ## Bad practice - This should be part of NetworkML
            ## - no random numbers in simulations, only in generators.
            randomsynindex = int(uniform(0,len(synapse_list)))
            synapse = synapse_list[randomsynindex]
            print "Connecting ", synapse.path
            ## put in i also, because randomly a synapse might get connected twice.
            spiketable = moose.TimeTable(synapse.path+'/tt_'+str(i))
            #### SynChan's synapse MsgDest takes time as its argument.
            ## Thus spiketable should contain a list of spike times.
            spiketable.connect("event", synapse,"synapse")
            if self.synName == mitgran_ampa:
                synapse2 = synapse_list2[randomsynindex]
                spiketable.connect("event", synapse2,"synapse")
            self.spikeTableList.append(spiketable)
            ## Presently only method 4 i.e. loading from file is supported for TimeTable.
            ## Hence the need to write a single value in a file.
            fn = 'temp_spikefile.txt'
            f = open(fn,'w')
            if self.stagger:
                f.write(str(i*self.synTime+SETTLETIME))
            else:
                f.write(str(SETTLETIME))
            f.close()
            spiketable.filename = fn
            os.remove(fn)
            if self.Vclamp:
                    self.synITable = setupTable('synITable',synapse,'Ik')
        if self.Vclamp:
            #### I block Ca and K channels in the cell
            #### We don't want these channels to be active when cell is held at 0mV; Na inactivates.
            blockChannels(self.cell, ['K','Ca']) # in moose_utils.py
            ## PID gain: I think for Davison 4-comp-mitral/granule: 0.5e-5 # optimal gain
            ## too high 0.5e-4 drives it to oscillate at high frequency,
            ## too low 0.5e-6 makes it have an initial overshoot (due to Na channels?)
            ## But for BBmit1993, gain of 1e-6 is optimal
            self.PIDITable = setup_vclamp(self.soma, '_somavclamp', 0,
                self.synRUNTIME, VCLAMP_V, gain=1e-6)
            print "Connected Vclamp to soma at",VCLAMP_V,"mV."

    def run(self):
        #printCellTree(self.cell)
        resetSim(self.context, SIMDT, PLOTDT)
        self.context.step(self.synRUNTIME) # each synapse is given 500ms!

        ## Paper figure 2 for granule cell electrophysiology
        fig = figure(figsize=(columnwidth/2.0,linfig_height/2.0),dpi=300,facecolor='w') # 'none' is transparent
        ax = fig.add_subplot(111)
        starttime = SETTLETIME*1000-10
        timevec = linspace(0.0,self.synRUNTIME,len(self.somaVm))*1000 - starttime
        self.somaVm = array(self.somaVm)*1000
        plot(timevec, self.somaVm,'-k', label='soma', linewidth=plot_linewidth, linestyle='-') # ms and mV
        xmin,xmax,ymin,ymax = beautify_plot(ax,x0min=True,y0min=False,drawxaxis=True,drawyaxis=True,\
            xticks=arange(0,self.synRUNTIME*1000-starttime,200),yticks=[])
        ax.set_xlim(0,800)
        ax.set_yticks([ymin,0,ymax])
        axes_labels(ax,'time (ms)','Vm(mV)',fontsize=label_fontsize,xpad=1,ypad=-5.5) # sets fontsize for ticks and labels
        #fig.text(0.04,0.7,'Vm (mV)', fontsize=label_fontsize, rotation='vertical', transform = fig.transFigure)
        
        if self.cellName == 'granule':
            ## Paper figure 2 for granule cell electrophysiology
            ## zoomed inset showing synaptic integration
            ## bbox_to_anchor can be instance of BboxBase, (x, y, width, height of the bbox), or (x, y) with width=height=0
            #axins = zoomed_inset_axes(ax, zoom=10, bbox_to_anchor=(0.43,1.0), bbox_transform=ax.transAxes)
            axins = inset_axes(ax, width=0.4, height=0.6, bbox_to_anchor=(0.5,1.2), bbox_transform=ax.transAxes)
            ## set the linewidth of the axes (spines)
            for s in axins.spines.values():
                 s.set_linewidth(axes_linewidth)
            ## Important: you need to plot the original plot again!!!
            plot(timevec, self.somaVm,'-k', label='soma', linewidth=plot_linewidth, linestyle='-') # ms and mV        
            zoom_xmin,zoom_xmax = 210, 250
            axins.set_xlim(zoom_xmin,zoom_xmax)
            ## take out the Vm data that is in the zoom region
            ## timevec may not have zoom_xmin exactly, so use index for value just below it, similar for zoom_xmax
            zoom_data = self.somaVm[where(timevec<zoom_xmin)[0][-1]:where(timevec>zoom_xmax)[0][0]]
            axins.set_ylim(min(zoom_data),max(zoom_data))
            ### autoscale will not work -- as the full data, not a subset is plotted !!
            #axins.autoscale(enable=True, axis='y', tight=True)
            #axins.autoscale_view(tight=True, scaley=True)
            #ymin, ymax = axins.get_yaxis().get_view_interval()
            #axins.set_ylim(ymin,ymax)
            axins.set_xticks([])
            axins.set_yticks([])
            mark_inset(ax, axins, loc1=3, loc2=4, fc="none", ec="0.5", linewidth=axes_linewidth)
            #add_scalebar(axins,matchx=False,matchy=False,hidex=False,hidey=False,\
            #    sizex=20,labelx='20 ms',sizey=1,labely='  1 mV',\
            #    bbox_to_anchor=[1.3,0.0],bbox_transform=axins.transAxes)
            axes_labels(axins,' {:3.0f} ms'.format(zoom_xmax-zoom_xmin),\
                '     {:3.2f} mV'.format(max(zoom_data)-min(zoom_data)),\
                xpad=-6,ypad=-6,fontsize=label_fontsize-2)
        elif self.cellName=='mitral':
            #plot(timevec, array(self.var)*1000,',-b', label='Sec dend Vm')
            plot(timevec, array(self.tuftBaseTable)*1000,'-b',\
                label='tuft base', linewidth=plot_linewidth, linestyle='--')
            plot(timevec, array(self.tuftTable)*1000,'-g',\
                label='tuft', linewidth=plot_linewidth, linestyle='-.')
            biglegend()
            
            ################ Paper figure 2
            if nerve_shock:
                fig2 = figure(figsize=(columnwidth/2.0,linfig_height/2.0),dpi=300,facecolor='w') # 'none' is transparent
                ax = fig2.add_subplot(111)
                plot(timevec, self.somaVm,'-k', label='soma', linewidth=linewidth, linestyle='-') # ms and mV
                plot(timevec, array(self.tuftBaseTable)*1000,'-b',\
                    label='tuft base', linewidth=linewidth, linestyle='--')
                beautify_plot(ax,drawxaxis=False,drawyaxis=False,xticks=[],yticks=[])
                add_scalebar(ax,matchx=False,matchy=False,hidex=True,hidey=True,\
                    sizex=3,labelx='3 ms',sizey=20,labely='20 mV',\
                    bbox_to_anchor=[0.3,0.6],bbox_transform=ax.transAxes)
                ax.set_xlim(SETTLETIME*1000-10,(SETTLETIME+15e-3)*1000) #ms
                ax.set_ylim(-70,50)
                fig2.tight_layout()
                fig2.savefig('../figures/connectivity/cells/mitral_spikeinit.png',dpi=fig2.dpi)
                fig2.savefig('../figures/connectivity/cells/mitral_spikeinit.svg',dpi=fig2.dpi)

        fig.tight_layout()
        savefig('../figures/connectivity/cells/'+self.pseudoSynName+'.png',dpi=fig.dpi)
        savefig('../figures/connectivity/cells/'+self.pseudoSynName+'.svg',dpi=fig.dpi)

        #figure(facecolor='w')
        #plot(timevec, self.var,',-b', label='Gk')
        if self.Vclamp:
            figure(facecolor='w')
            plot(timevec, self.synITable,'g-,')
            title("synapse current")
            xlabel('time (ms)',fontsize='large')
            ylabel('I (A)',fontsize='large')
            figure(facecolor='w')
            plot(timevec, self.PIDITable,'g-,')
            title("V clamp current")
            xlabel('time (ms)',fontsize='large')
            ylabel('I (A)',fontsize='large')
        show()
        
if __name__ == "__main__":
    seed([122.0])
    numsyns = int(sys.argv[4])
    stagger = False
    if numsyns > 1 and sys.argv[6] == 'staggered':
        stagger = True
    if sys.argv[3]=='Vclamp': Vclamp = True
    elif sys.argv[3]=='Iclamp': Vclamp = False
    else:
        print "You must specify Vclamp or Iclamp. Using Iclamp by default."
        Vclamp = False
    if sys.argv[2]=='granule_mitral':
        libsynname = 'granule_mitral_GABA'
        synfactor = GABA_factor*strongsynfactorinh
    elif sys.argv[2]=='mitral_granule':
        libsynname = mitgran_ampa
        synfactor = AMPA_factor
    else:
        libsynname = sys.argv[2]
        synfactor = 1.0
    if len(sys.argv)>=8: comp_name = sys.argv[7]
    else: comp_name = None
    cell = CellTest(cellname=sys.argv[1],pseudosynname=sys.argv[2],libsynname=libsynname,
        synfactor=synfactor,Vclamp=Vclamp, numsyns=numsyns,
        syntime=float(sys.argv[5]),stagger=stagger, comp_name=comp_name)
    cell.run()
