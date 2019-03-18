#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import math

# The PYTHONPATH should contain the location of moose.py and _moose.so
# files.  Putting ".." with the assumption that moose.py and _moose.so
# has been generated in ${MOOSE_SOURCE_DIRECTORY}/pymoose/ (as default
# pymoose build does) and this file is located in
# ${MOOSE_SOURCE_DIRECTORY}/pymoose/examples
# sys.path.append('..\..')
try:
    import moose
except ImportError:
    print "ERROR: Could not import moose. Please add the directory containing moose.py in your PYTHONPATH"
    import sys
    sys.exit(1)

from channelConstants import *

VKCa = -80.0e-3 # Volts # kca3.mod has a vshift=-10mV in addition to ek=-70mV

GKCa = 142.0*sarea # Siemens, from mit4.hoc

def calc_KCa_alpha_y(v,Ca):
    #return math.exp((v-65e-3)/27e-3) * 500.0e3*(0.015-Ca)/(math.exp((0.015-Ca)/0.0013)-1) # This is from Bhalla and Bower's 1993 paper: perhaps the v-65e-3 is a typo. It should be as below.
    return math.exp((v+70e-3)/27e-3)*1e3 * 500.0*(0.015-Ca)/(math.exp((0.015-Ca)/0.0013)-1) # This is from kca3.mod: different from Bhalla and Bower 1993 paper. see above

def calc_KCa_beta_y(v,Ca):
    return 0.05e3

CaMIN = 0.0
CaMAX = 1.0e-2 # mol/m^2 same as millimol/litre
CaNDIVS = 100

class KCaChannel(moose.HHChannel2D):
    """KCa channel inherits from HHChannel2D."""
    def __init__(self, *args):
        """Setup the KCa channel with defaults"""
        moose.HHChannel2D.__init__(self,*args)
        ### Since this channel is table based and this script is for creation of KCaA.dat and KCaB.dat,
        ### while KCaMPIChannel.py is for reading in these tables,
        ### we hard code the VMIN, VMAX and NDIVS so that those in globalConstants.py
        ### do not spoil the reading in of already created tables.
        VMIN = -0.1 # V
        VMAX = 0.05 # V
        NDIVS = 150
        dv = (VMAX-VMIN)/NDIVS
        ## For HHChannel2D Ek, Gbar, etc don't get set via python assignments!!!
        ## Have to use moose shell commands.
        #self.Ek = VKCa
        self.getContext().runG('setfield '+self.path+' Ek '+str(VKCa))
        #self.Gbar = GKCa
        self.getContext().runG('setfield '+self.path+' Gbar '+str(GKCa))
        self.addField('ion')
        self.setField('ion','K')
        self.addField('ionDependency')
        self.setField('ionDependency','Ca')
        ## This will create HHGate2D instance xGate inside the KCa channel.
        ## (2D since it is within HHChannel2D)
        self.Xpower = 1
        ## VOLT_C1_INDEX: VOLTAGE message specifies x variable
        ## and CONCEN1 variable specifies y variable of X_A and X_B tables;
        ## pg 323, sec 19.4.6 of Book of Genesis
        self.Xindex = "VOLT_C1_INDEX"
        ## xGate was already created and wrapped when Xpower was set non-zero
        #self.xGate = moose.HHGate2D(self.path + "/xGate")
        #self.getContext().runG("showfield "+self.path +"/xGate/A -all")
        self.xGate.A.xmin = VMIN
        self.xGate.A.xmax = VMAX
        #self.xGate.A.xdivs = NDIVS # these get overridden by the number of values in the table
        self.xGate.B.xmin = VMIN
        self.xGate.B.xmax = VMAX
        #self.xGate.B.xdivs = NDIVS # these get overridden by the number of values in the table

        ### HHGate2D is not wrapped properly in pyMOOSE.
        ### ymin, ymax and ydivs are not exposed.
        ### Setting them creates new and useless attributes within HHGate2D without warning!
        ### Hence use runG to set these via Genesis command
        self.getContext().runG("setfield "+self.path+"/xGate/A"+\
            #" ydivs "+str(CaNDIVS)+\ # these get overridden by the number of values in the table
            " ymin "+str(CaMIN)+\
            " ymax "+str(CaMAX))
        self.getContext().runG("setfield "+self.path+"/xGate/B"+\
            #" ydivs "+str(CaNDIVS)+\ # these get overridden by the number of values in the table
            " ymin "+str(CaMIN)+\
            " ymax "+str(CaMAX))

        selfdir = os.path.dirname(__file__)
        if selfdir != '': selfdir += os.sep
        ftableA = open(selfdir+"KCaA.dat","w")
        ftableB = open(selfdir+"KCaB.dat","w")
        v = VMIN
        dCa = (CaMAX-CaMIN)/CaNDIVS
        for i in range(NDIVS+1):
            Ca = CaMIN
            for j in range(CaNDIVS+1):
                alpha = calc_KCa_alpha_y(v,Ca)
                ftableA.write(str(alpha)+" ")
                ftableB.write(str(alpha+calc_KCa_beta_y(v,Ca))+" ")
                Ca += dCa
            ftableA.write("\n")
            ftableB.write("\n")
            v += dv
        ftableA.close()
        ftableB.close()

        ### PRESENTLY, Interpol2D.cpp in MOOSE only allows loading via a data file,
        ### one cannot set individual entries A[0][0] etc.
        ### Thus pyMOOSE also has not wrapped Interpol2D
        self.getContext().runG("call "+self.path+"/xGate/A load "+selfdir+"KCaA.dat 0")
        self.getContext().runG("call "+self.path+"/xGate/B load "+selfdir+"KCaB.dat 0")
        
        # Test print of table
        #self.getContext().runG("call "+self.path+"/xGate/A print "+selfdir+"KCaA_out.dat")
        #self.getContext().runG("call "+self.path+"/xGate/B print "+selfdir+"KCaB_out.dat")

        # calc_mode is LIN_INTERP i.e. 1 by default, so no need to set it as below.
        #self.getContext().runG("setfield "+self.path+" X_A->calc_mode 1")
        #self.getContext().runG("setfield "+self.path+" X_B->calc_mode 1")

