#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import math

sys.path.extend(["..","../neuroml","../channels",\
    "../synapses","../cells","../networks","../generators"])

from networkConstants import *
from stimuliConstants import * # has MAXNUMAVG_GRANS

from load_channels import *
from load_synapses import *
from load_cells import *
from moose.neuroml import *
from moose_utils import * # also 'import *'-s the moose.utils file

from pylab import * # part of matplotlib that depends on numpy but not scipy

class OBNetwork:

    ## unique str = 'morphs', etc and mpirank are passed,
    ## so that temp files of morphs, pulses, etc do not overlap
    def __init__(self, OBfile, synchan_activation_correction,
        tweaks, mpirank, uniquestr, granfilebase, resptuned=False, spiketable=False):
        self.mpirank = mpirank
        self.uniquestr = uniquestr
        self.context = moose.PyMooseBase.getContext()
        self.granfilebase = granfilebase
        self.spiketable = spiketable
        ## netseed is between 'seed' and '_' / '.xml' in OBfile
        ## ensure there is an _ before .xml, to consider both endings.
        OBfile_underscored = OBfile.replace('.xml','_.xml')
        self.netseed = float((OBfile_underscored.split('seed')[1]).split('_')[0])
        load_channels()
        ## synchan_activation_correction depends on SIMDT,
        ## hence passed all the way down to
        ## granule_mitral_GABA synapse from the top level
        load_synapses(synchan_activation_correction)
        self.cellSegmentDict = load_cells()
        
        filename = OBfile
        NML = NetworkML({'temperature':CELSIUS})
        ## Below returns populationDict = { 'populationname1':(cellname,{instanceid1:moosecell, ... }) , ... }
        ## and projectionDict = { 'projectionname1':(source,target,[(syn_name1,pre_seg_path,post_seg_path),...]) , ... }
        (self.populationDict,self.projectionDict) = \
            NML.readNetworkMLFromFile(filename,self.cellSegmentDict,params=tweaks)

        if VARY_GRANS_RMP:
            self.set_granules_rmp(gran_RMP_SD)
        if VARY_PGS_RMP:
            self.set_pgs_rmp(PG_RMP_SD) # let PGs vary also
            self.set_different_pgs() # set cells as PG_LTS or PG
            ## Vary synapse/channel Gbar-s only after setting them as PG_LTS or PG above
            self.set_pgs_input_conductances() # set input conductance based on PG type
            #self.set_pgs_tca_var() # changes TCa for all PGs
            ## for PG_LTS, varying ELeak or TCa doesn't change RMP, varying KCa does
            self.set_pgs_kca_var() # changes KCa for PG_LTSs

        #### Granule baseline synapses and connected timetables are created
        #### by the networkml reader from the network xml file.
        #### But the firefiles need to be loaded into timetables.
        self.connect_granule_baselines_to_files('granule_baseline')

        #### mitral baseline synapses and connected timetables are created
        #### by the networkml reader from the network xml file.
        #### But the firefiles need to be loaded into timetables.
        #print "Connecting 10000 files per mit"
        ## connect baseline inhibition to mitrals at lateral dendrites
        #self.connect_granule_baselines_to_files('mitral_baseline')
        ## connect baseline inhibition to mitrals at tuft dendrites
        #self.connect_granule_baselines_to_files('tuft_mitral')

        #printCellTree(moose.Cell('granules_singles_0'))
        # unordered dictionary self.populationDict['mitrals'][1] is
        # renamed to a more friendly self.mitralTable
        self.mitralTable = self.populationDict['mitrals'][1]
        self.setupMitralRecordingTables()
        
        #printNetTree()
        print "OB Network from",OBfile,"loaded!"

    ## odor_morphs, odor_pulses, activdep_inhibition.py, inhibition.py, etc.
    ## depend on granule baselines to be connected by this network loader.
    ## So do NOT remove and place below function in some other file.
    def connect_granule_baselines_to_files(self, basename):
        for projname in self.projectionDict.keys():
            if basename in projname:
                for i,proj in enumerate(self.projectionDict[projname][2]):
                    ## just wrap the existing timetable created when loading in neuroml
                    ## it has a field fileNumbers that has which lines to take from the firingfile
                    tt = moose.TimeTable(proj[2]+'/'+proj[0]+'_tt') # post_segment_path+'/'+syn_name+'_tt'
                    ## Each trial / avgnum has its own set of granule baseline spike trains
                    ## identified by (self.mpirank-1)%MAXNUMAVG_GRANS
                    attach_spikes(self.granfilebase+'_'+\
                        str((self.mpirank-1)%MAXNUMAVG_GRANS), tt, self.uniquestr+str(self.mpirank))
                    if i%1000==0: print self.granfilebase,"for proj#",i,"at rank",self.mpirank
                print basename,"connected at",self.mpirank

    def setupMitralRecordingTables(self):
        for mitral in self.populationDict['mitrals'][1].values():
            ## Setup the tables to pull data
            ## assumes at least one soma and takes the first!
            mitral.soma = moose.Compartment(get_matching_children(mitral, ['Soma','soma'])[0])
            mitral._vmTableSoma = setupTable("vmTableSoma",mitral.soma,'Vm')
            mitral.dend = moose.Compartment(get_matching_children(mitral, ['Seg0_sec_dendd4_1_264'])[0])
            ## assumes at least one sec_dend and takes the first!
            #mitral.dend = moose.Compartment(get_matching_children(mitral, ['sec_dend'])[0])
            mitral._vmTableDend = setupTable("vmTableDend",mitral.dend,'Vm')
            #mitral.glom = moose.Compartment(get_matching_children(mitral, ['Seg0_prim_dend_5'])[0])
            ## assumes at least one glom and takes the first!
            mitral.glom = moose.Compartment(get_matching_children(mitral, ['glom'])[0])
            mitral._vmTableGlom = setupTable("vmTableGlom",mitral.glom,'Vm')
            if self.spiketable:
                mitral._vmTableSoma.stepMode = TAB_SPIKE
                mitral._vmTableSoma.stepSize = THRESHOLD
                mitral._vmTableDend.stepMode = TAB_SPIKE
                mitral._vmTableDend.stepSize = THRESHOLD
                mitral._vmTableGlom.stepMode = TAB_SPIKE
                mitral._vmTableGlom.stepSize = THRESHOLD

    def set_granules_rmp(self,sd):
        ## set the seed here based on netseed,
        ## so that there is reproducible generation of random cell variations for a given network.
        seed([self.netseed])
        ## populationDict = { 'populationname1':(cellname,{instanceid1:moosecell, ... }) , ... }
        ## tweak_field() is from moose.utils
        ## important to sort cells by instance_ids for repeatability over runs (dict is unordered).
        ## some granules may have been 'tweak'-ed out, check if existent
        pop_keys = self.populationDict.keys()
        if 'granules_singles' in pop_keys:
            for gran in self.sortedbykey_dictvals(self.populationDict['granules_singles'][1]):
                ## allow 6 sd-s each side (range in Cang & Issacson 2003)
                varEm = self.trunc_2side_stdnormal(6.0)*sd
                tweak_field(gran.path+'/##[TYPE=Compartment]', 'Em', 'Em+'+str(varEm))
        if 'granules_joints' in pop_keys:
            for gran in self.sortedbykey_dictvals(self.populationDict['granules_joints'][1]):
                ## allow 6 sd-s each side (range in Cang & Issacson 2003)
                varEm = self.trunc_2side_stdnormal(6.0)*sd
                tweak_field(gran.path+'/##[TYPE=Compartment]', 'Em', 'Em+'+str(varEm))
        if 'granules_multis' in pop_keys:
            for gran in self.sortedbykey_dictvals(self.populationDict['granules_multis'][1]):
                ## allow 6 sd-s each side (range in Cang & Issacson 2003) 1.5mV SD
                varEm = self.trunc_2side_stdnormal(6.0)*sd
                tweak_field(gran.path+'/##[TYPE=Compartment]', 'Em', 'Em+'+str(varEm))

    def set_pgs_rmp(self,sd):
        ## set the seed here based on netseed,
        ## so that there is reproducible generation of random cell variations for a given network.
        seed([self.netseed+1])
        ## populationDict = { 'populationname1':(cellname,{instanceid1:moosecell, ... }) , ... }
        ## tweak_field() is from moose.utils
        ## important to sort cells by instance_ids for repeatability over runs (dict is unordered).
        ## PGs may have been 'tweak'-ed out, check if existent
        pop_keys = self.populationDict.keys()
        if 'PGs' in pop_keys:
            for gran in self.sortedbykey_dictvals(self.populationDict['PGs'][1]):
                ## allow 8 sd-s each side (range in Shao et al 2009) 1mV SD
                varEm = self.trunc_2side_stdnormal(8.0)*sd
                tweak_field(gran.path+'/##[TYPE=Compartment]', 'Em', 'Em+'+str(varEm))

    def set_different_pgs(self):
        ## 2/3 of the PG cells are the loaded LTS,
        ## but the rest are converted to plateauing cells, by changing their channel values.
        
        ## set the seed here based on netseed,
        ## so that there is reproducible generation of random cell variations for a given network.
        seed([self.netseed+2])
        ## populationDict = { 'populationname1':(cellname,{instanceid1:moosecell, ... }) , ... }
        ## tweak_field() is from moose.utils
        ## important to sort cells by instance_ids for repeatability over runs (dict is unordered).
        ## PGs may have been 'tweak'-ed out, check if existent
        pop_keys = self.populationDict.keys()
        if 'PGs' in pop_keys:
            for pg in self.sortedbykey_dictvals(self.populationDict['PGs'][1]):
                ## 2/3 of the PG cells are LTS (McQuiston & Katz,2001)
                if uniform(0.0,1.0)<0.67:
                    ## 2013 LTS PG cell -- this field is needed by later functions to query pg type
                    pg.addField('PG_type')
                    pg.setField('PG_type','PG_LTS')
                else:
                    ## 2010 plateauing PG cell -- this field is needed by later functions to query pg type
                    pg.addField('PG_type')
                    pg.setField('PG_type','PG')
                    ## Get the compartments in the cell
                    comp_id_list = moose.context.getWildcardList(pg.path+'/##[TYPE=Compartment]', True)
                    for comp_id in comp_id_list:
                        pg_comp = moose.Compartment(comp_id)
                        if 'soma' in pg_comp.name: soma_chan = True
                        else: soma_chan = False
                        ## Values in cell neuroml files are in mS/cm^2. Convert & use compartment surface area
                        comp_factor = 1e1*math.pi*pg_comp.diameter*pg_comp.length ## 1e1 is S/m^2 from mS/cm^2

                        ## Get the channels in the compartment
                        chan_id_list = moose.context.getWildcardList(pg_comp.path+'/##[TYPE=HHChannel]', True)
                        for chan_id in chan_id_list:
                            pg_chan = moose.HHChannel(chan_id)
                            if 'Na' in pg_chan.name:
                                if soma_chan: pg_chan.Gbar = 20.0*comp_factor
                                else: pg_chan.Gbar = 5.0*comp_factor
                            if 'K2' in pg_chan.name:
                                if soma_chan: pg_chan.Gbar = 40.0*comp_factor
                                else: pg_chan.Gbar = 20.0*comp_factor
                            if 'Ih' in pg_chan.name:
                                if soma_chan: pg_chan.Gbar = 0.3*comp_factor
                                else: pg_chan.Gbar = 0.1*comp_factor
                            if 'TCa' in pg_chan.name:
                                if soma_chan: pg_chan.Gbar = 3.0*comp_factor
                                else: pg_chan.Gbar = 1.0*comp_factor
                            if 'KA' in pg_chan.name:
                                if soma_chan: pg_chan.Gbar = 5.0*comp_factor
                                else: pg_chan.Gbar = 5.0*comp_factor
                        chan_id_list = moose.context.getWildcardList(pg_comp.path+'/##[TYPE=HHChannel2D]', True)
                        for chan_id in chan_id_list:
                            pg_chan = moose.HHChannel2D(chan_id)
                            if 'Kca' in pg_chan.name:
                                if soma_chan: pg_chan.Gbar = 0.0*comp_factor
                                else: pg_chan.Gbar = 0.0*comp_factor

    def set_pgs_input_conductances(self):
        ## set the seed here based on netseed,
        ## so that there is reproducible generation of random cell variations for a given network.
        seed([self.netseed+3])
        ## populationDict = { 'populationname1':(cellname,{instanceid1:moosecell, ... }) , ... }
        ## tweak_field() is from moose.utils
        ## important to sort cells by instance_ids for repeatability over runs (dict is unordered).
        ## PGs may have been 'tweak'-ed out, check if existent
        pop_keys = self.populationDict.keys()
        if 'PGs' in pop_keys:
            for pg in self.sortedbykey_dictvals(self.populationDict['PGs'][1]):
                if 'LTS' not in pg.getField('PG_type'):
                    id_list = moose.context.getWildcardList(pg.path+'/##[TYPE=SynChan]', True)
                    ## The only synapses on the PGs are ORN->PG and mitral->PG. Reduce for both.
                    for moose_id in id_list:
                        pg_syn = moose.SynChan(moose_id)
                        pg_syn.Gbar = pg_syn.Gbar*0.45/1.25 # 0.45nS for 2010 PG, instead of 1.25nS for 2013 PG.

    def set_pgs_tca_var(self):
        ## set the seed here based on netseed,
        ## so that there is reproducible generation of random cell variations for a given network.
        seed([self.netseed+4])
        ## populationDict = { 'populationname1':(cellname,{instanceid1:moosecell, ... }) , ... }
        ## tweak_field() is from moose.utils
        ## important to sort cells by instance_ids for repeatability over runs (dict is unordered).
        ## PGs may have been 'tweak'-ed out, check if existent
        pop_keys = self.populationDict.keys()
        if 'PGs' in pop_keys:
            for pg in self.sortedbykey_dictvals(self.populationDict['PGs'][1]):
                id_list = moose.context.getWildcardList(pg.path+'/##[TYPE=HHChannel]', True)
                ## Get the HHChannel TCa_d
                for moose_id in id_list:
                    pg_chan = moose.HHChannel(moose_id)
                    if 'TCa' in pg_chan.name:
                        pg_chan.Gbar = pg_chan.Gbar*uniform(0.5,1.0) # 50% to 100% variation of TCa

    def set_pgs_kca_var(self):
        ## set the seed here based on netseed,
        ## so that there is reproducible generation of random cell variations for a given network.
        seed([self.netseed+5])
        ## populationDict = { 'populationname1':(cellname,{instanceid1:moosecell, ... }) , ... }
        ## tweak_field() is from moose.utils
        ## important to sort cells by instance_ids for repeatability over runs (dict is unordered).
        ## PGs may have been 'tweak'-ed out, check if existent
        pop_keys = self.populationDict.keys()
        if 'PGs' in pop_keys:
            for pg in self.sortedbykey_dictvals(self.populationDict['PGs'][1]):
                if 'LTS' in pg.getField('PG_type'):
                    id_list = moose.context.getWildcardList(pg.path+'/##[TYPE=HHChannel2D]', True)
                    ## Get the HHChannel2D KCa
                    for moose_id in id_list:
                        pg_chan = moose.HHChannel(moose_id)
                        if 'Kca_mit_usb_pg' in pg_chan.name:
                            pg_chan.Gbar = pg_chan.Gbar*uniform(0.5,1.5) # 50% to 150% variation of KCa

    def sortedbykey_dictvals(self,instances_dict):
        dict_keys = instances_dict.keys()
        dict_keys.sort()
        return [instances_dict[key] for key in dict_keys]

    def trunc_2side_stdnormal(self,numsd):
        randomnum = standard_normal()
        if randomnum<-numsd:
            randomnum = -numsd
        elif randomnum>numsd:
            randomnum = numsd
        return randomnum
