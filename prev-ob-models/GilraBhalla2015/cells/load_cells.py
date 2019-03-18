#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import math

sys.path.extend(["..","../channels"])
from channelConstants import * # for CELSIUS, the global temperature

import moose
from moose.neuroml import *

def load_cells():
    MML = MorphML({'temperature':CELSIUS})
    ## writing ../cells/ below allows load_cells() to be called from within any sister directory

    ## original mitral cell of Upi with half RM (half tau), -58mV resting potential
    ## more Na & K in primary dend, some Na & K in tuft, electrode leak removed.
    #mitral_dict = MML.readMorphMLFromFile('../cells/mitral_bbmit1993davison_neuroML_L1_L2_L3_mod.xml',{})
    ## modified as above + for spike initiation at soma for small tuft input and at tuft for large tuft input:
    ## RA reduced; added hillock and initial segment (with Migliore & Shepherd's Na channel), tuft and prim dend more excitable
    mitral_dict = MML.readMorphMLFromFile('../cells/mitral_bbmit1993davison_neuroML_L1_L2_L3_mod_withspikeinit.xml',{})
    #mitral_dict = MML.readMorphMLFromFile('../cells/mitral_bbmit1993davison_neuroML_L1_L2_L3.xml',{})
    #mitral_dict = MML.readMorphMLFromFile('../cells/mitral_bbmit1993davison_neuroML_TEST_L1_L2_L3.xml',{})
    ## modified mitral cell with changes as above plus: RM is 1/5th,
    ## Migliore & Shepherd Na channels, special hillock MS Na channel
    #mitral_dict = MML.readMorphMLFromFile('../cells/mitral_bbmit1993davisonMS_neuroML_L1_L2_L3.xml',{})

    ## Low threshold spiking cell, LTS with shoulder and plateau potential due to extra TCa (low threshold Ca).
    PG_dict = MML.readMorphMLFromFile('../cells/PG_aditya2013unified_neuroML_L1_L2_L3.xml',{})
    ### For the 2012 cell, to get the Ca conc correct, below code was used. Obsolete with the new 2013 cell.
    ### Low threshold spiking cell, LTS with shoulder and plateau potential due to extra TCa (low threshold Ca).
    ### Be sure to uncomment the CaBasal setting below too!
    ### For just the new 2012 PG cell set CaPool.CaBasal to 0.05 mM, else the PG fires a settling burst at t=0.
    ### PG_dict['PG'] = { segid1 : [ segname,(proximalx,proximaly,proximalz),
    ###  (distalx,distaly,distalz),diameter,length,[potential_syn1, ... ] ] , ... }
    ### This is a complete hackjob, ideally, I should use a different CaPool.
    #PG_dict = MML.readMorphMLFromFile('../cells/PG_aditya2012_neuroML_L1_L2_L3.xml',{})
    #for seg in PG_dict['PG'].values():
    #    segname = seg[0]
    #    segpath = '/library/PG/'+segname
    #    if not moose.context.exists(segpath):
    #        print "Missing compartment",segpath
    #        sys.exit(1)
    #    compartment = moose.Compartment(segpath)
    #    for child in compartment.getChildren(compartment.id): # channels and synapses
    #        if moose.Neutral(child).className in ['CaConc']:
    #            CaPool = moose.CaConc(child)
    #            CaPool.CaBasal = 0.05 # mM
    ### hackjob for PG cell 2012 ends.

    granule_dict = MML.readMorphMLFromFile('../cells/granule_granadityaMS2007_neuroML_L1_L2_L3.xml',{})

    ## make one dictionary out of above dicts.
    cellSegmentList = {}
    cellSegmentList.update(mitral_dict)
    cellSegmentList.update(PG_dict)
    cellSegmentList.update(granule_dict)
    return cellSegmentList
