#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import math

sys.path.extend(["..","../channels/"])

from NaChannel import *
from KFastChannel import *
from KSlowChannel import *
from CaLChannel import *
from KAChannel import *
## use KCaChannel instead of KCaMPIChannel, in a non-MPI i.e. non-parallel run,
## to generate the KCaA.dat and KCaB.dat files.
#from KCaChannel import *
from KCaMPIChannel import *
## use KCaChannel_PG instead of KCaMPIChannel_PG, in a non-MPI i.e. non-parallel run,
## to generate the KCaA_PG.dat and KCaB_PG.dat files.
#from KCaChannel_PG import *
from KCaMPIChannel_PG import *
from CaPool import *
from KMChannel import *
from CaTChannel import *

from NaMitChannelMS import *
from KAChannelMS import *
from KDRChannelMS import *

import moose
from moose.neuroml import *

FARADAY = 96154.0 # Coulombs # from cadecay.mod : 1/(2*96154.0) = 5.2e-6 which is the Book of Genesis / readcell value
#FARADAY = 96485.3415 # Coulombs # from Wikipedia

def load_channels():
    # Na channels in /library should be called *Na* not *na*. No other channel should have Na in its name.
    # I search later for Na channels by string *Na*
    NaChannel("/library/Na_mit_usb")
    KFastChannel("/library/K2_mit_usb") ## CAUTION: K2 is Kfast
    KSlowChannel("/library/K_mit_usb") ## CAUTION: K is Kslow
    CaLChannel("/library/LCa3_mit_usb") # L-type Ca channel (high threshold)
    KAChannel("/library/KA_bsg_yka")
    KCaChannel("/library/Kca_mit_usb")
    KCaChannel_PG("/library/Kca_mit_usb_pg")
    CaPool("/library/Ca_mit_conc")
    KMChannel("/library/KM_bsg_upi")
    CaTChannel("/library/TCa_rat_ag") # T-type Ca channel (low threshold)
    # NaMitChannelMS(sh, Vna, MOOSEpathname)
    # mit Na and gran Na channels are different in sh and ENa. 
    NaMitChannelMS(10e-3, 50e-3, "/library/Na_mit_ms")
    NaMitChannelMS(0e-3, 50e-3, "/library/Na_mit_initialsegment_ms")
    # for Granule cell, maintaining the old name Na_rat_ms to avoid changing .p file (now invalid argument!),
    # even though this Na channel for granule cell is from Migliore and Shepherd 2008
    NaMitChannelMS(15e-3, 60e-3, "/library/Na_rat_ms")
    KAChannelMS("/library/KA_ms")
    KDRChannelMS("/library/KDR_ms")
    #self.context.readNeuroML(../channels/IhChannel.xml,"/library/Ih_cb")
    CML = ChannelML({'temperature':CELSIUS})
    CML.readChannelMLFromFile('../channels/Ih_cb.xml')
    CML.readChannelMLFromFile('../channels/TCa_d.xml')
    ## extras from neuroConstruct examples
    #CML.readChannelMLFromFile('../channels/CaHVA_Chan.xml')
    #CML.readChannelMLFromFile('../channels/CaL_Chan.xml')

def connect_CaConc(compartment_list):
    context = moose.PyMooseBase.getContext()
    #### Connect the Ca pools and channels
    #### Ca channels should have an extra field called 'ion' defined and set in MOOSE.
    #### Ca dependent channels like KCa should have an extra field called 'ionDependency' defined and set in MOOSE.
    #### Am connecting these at the very end so that all channels and pools have been created
    for compartment in compartment_list:
        if context.exists(compartment.path+'/Ca_mit_conc'): # Ca Pool
            caconc = moose.CaConc(compartment.path+'/Ca_mit_conc')
            for child in compartment.getChildren(compartment.id):
                neutralwrap = moose.Neutral(child)
                if neutralwrap.className == 'HHChannel':
                    channel = moose.HHChannel(child)
                    ### If 'ion' field is not present, the Shell returns '0', cribs and prints out a message but it does not throw an exception
                    if channel.getField('ion') == 'Ca':
                        channel.connect('IkSrc',caconc,'current')
                        #print 'Connected ',channel.path
                if neutralwrap.className == 'HHChannel2D':
                    channel = moose.HHChannel2D(child)
                    ### If 'ionDependency' field is not present, the Shell returns '0', cribs and prints out a message but it does not throw an exception
                    if channel.getField('ionDependency') == 'Ca':
                        caconc.connect('concSrc',channel,'concen')
                        #print 'Connected ',channel.path


