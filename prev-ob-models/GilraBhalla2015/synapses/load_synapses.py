#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

sys.path.extend(["..","../channels","../synapses"])

from synapseConstants import *
from channelConstants import * # has CELSIUS
import moose
from moose.neuroml import *

from mitral_granule_AMPA import *
from mitral_granule_saturatingAMPA import *
from mitral_granule_NMDA import *
from mitral_granule_saturatingNMDA import *
from granule_mitral_GABA import *
from mitral_self_GABA import *

def load_synapses(synchan_activation_correction):
    mitral_granule_AMPA("/library/mitral_granule_AMPA")
    mitral_granule_saturatingAMPA("/library/mitral_granule_saturatingAMPA")
    mitral_granule_NMDA("/library/mitral_granule_NMDA")
    mitral_granule_saturatingNMDA("/library/mitral_granule_saturatingNMDA")
    granule_mitral_GABA(synchan_activation_correction,"/library/granule_mitral_GABA")
    mitral_self_GABA("/library/mitral_self_GABA")     
    CML = ChannelML({'temperature':CELSIUS})
    CML.readChannelMLFromFile('../synapses/ORN_PG.xml')
    CML.readChannelMLFromFile('../synapses/ORN_mitral.xml')
    CML.readChannelMLFromFile('../synapses/mitral_PG.xml')
    CML.readChannelMLFromFile('../synapses/PG_mitral.xml')
    CML.readChannelMLFromFile('../synapses/SA_PG.xml')

def attachORNs(cell, ORNCellSynapse, numORNSyns):
    ##### Inhibitory Synapse NOT added as I am not sure if it is due to PG inhibition - Perhaps Cleland and Sethupathy also did not model it       
    ##### Excitatory + Inhibitory(not added) Synapse combo taken from the paper Djurisic etal 2008 JNeurosci.
    ##### Actually it's only an excitatory synapse, but they have used the inhibitory one to model the later time course.
    ##### Though this might be needed to account for PG cell inhibition?
    ##### Gbar-s have been set to ensure 16mV EPSP at the glom tuft as done in the paper.
    ##### Paper's defaults give only 8mV EPSP peak at glom tuft. So here multiplied by two. Cannot use supplemental table values as no cytoplasmic resistivity Ri in this model. Only axial resistance from glom to prim which doesn't matter much [maybe it does - loading rather than input resistance?]. Makes sense only if dendritic tree with many compartments.
    ##### Also no idea of how much to change the inhibitory part without Ri, so multiplied that also by 2

    synapse_list = []
    for compartment in cell.getChildren(cell.id): # compartments
        for child in moose.Neutral(compartment).getChildren(compartment): # channels and synapses
            if moose.Neutral(child).name in [ORNCellSynapse]:
                synapse_list.append(moose.SynChan(child))
    spikeTableList = []
    for i in range(numORNSyns):
        synapse = synapse_list[int(uniform(0,len(synapse_list)))] # Bad practice - This should be part of NetworkML - no random numbers in simulations, only in generators.
        print "Connecting ", synapse.path
        spiketable = moose.TimeTable(synapse.path+'/tt_'+str(i))
        #### SynChan's synapse MsgDest takes time as its argument. Thus spiketable should contain a list of spike times.
        spiketable.connect("event", synapse,"synapse")
        synapse.setWeight(0, 1) # 0th element in synaptic array set to weight 1
        spikeTableList.append(spiketable)
    return spikeTableList
