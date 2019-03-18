#!/usr/bin/env python
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

# ALL SI UNITS

class CaPool(moose.CaConc):
    def __init__(self, *args):
        moose.CaConc.__init__(self, *args)
        ## 0.05 is the value at which Ca settles in the PG cell,
        ## setting Ca_Base to 0.05 makes intial burstiness of PG cell 2012 go away,
        ## yet it makes no difference to its later firing!
        ## but setting it to 0.05 makes the mitral very low firing.
        ## From cadecay.mod but Book of Genesis says 5e-5 mol/m^3 is typical.
        self.CaBasal = 1e-5 # mol/m^3, same as mmol/litre = mMolar
        ## dunno why python wrapper has both CaBasal and Ca_Base
        ## setting Ca_Base to 0.05 makes no difference to PG cell 2012's
        ## firing or initial burstiness or Ca response.
        self.Ca_Base = 1e-5 # mol/m^3, same as mmol/litre = mMolar
        self.tau = 10e-3 # second # From cadecay.mod
        self.getContext().runG('setfield ' + self.path + ' ceiling 1e6')
        self.getContext().runG('setfield ' + self.path + ' floor 0.0')
        
    def connectCaChannels(self, channel_list):
        """Connects the Ca2+ channels in channel_list as a source of
        Ca2+ to the pool."""
        for channel in channel_list:
                if not hasattr(channel, 'connected_to_pool') or not channel.connected_to_pool:
                    channel.connect('IkSrc', self, 'current')
                    channel.connected_to_pool = True
                else:
                    print channel.path, 'already connected'
                
    def connectDepChannels(self, channel_list):
        """Connect channels in channel_list as dependent channels"""
        for channel in channel_list:
            if not hasattr(channel, 'connected_to_ca') or not channel.connected_to_ca:
                self.connect("concSrc", channel, "concen")
                #self.getContext().runG('addmsg '+self.path+' '+channel.path+' CONCEN1 Ca')
                channel.connected_to_ca = True
            else:
                print "WARNING: Ignoring non-KCaChannel", channel.path
# 
# capool.py ends here
