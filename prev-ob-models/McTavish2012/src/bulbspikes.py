# -*- coding: utf-8 -*-
#*****************************************************************************
#       Copyright (C) 2010 THOMAS MCTAVISH <Thomas.McTavish@yale.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import sys
import numpy
from neuronpy.util import spiketrain
import bisect

class BulbSpikes:
    """
    Object to store and retrieve mitral and granule cell olfactory bulb spikes.

    Works with a more generic :class: `spiketrain` object that deals with 
    all of the spikes in a simulation to separate out mitral cell and granule 
    cell somatic spikes. It is possible to get some simple diagnostics as well, 
    like a spike plot and frequencies.

    AUTHORS:

    - THOMAS MCTAVISH (2010-02-05): initial version
    """
    
    def __init__(self, sim_time=60000, num_mit=5, num_gran=100, mit_spikes=[],
    gran_spikes=[], syn_spikes={}):
        """
        Initialize member variables.
        
        :var num_mit: Number of mitral cells.
        
        :var num_gran: Number of granule cells.
        
        :var sim_time: Duration of the simulation in ms.
        
        :var mit_spikes: Mitral cell spikes as a list of lists.
        
        :var gran_spikes: Granule cell spikes as a list of lists.        
        """
        self.raw=None
        self.num_mit=num_mit
        self.num_gran=num_gran
        self.sim_time=sim_time
        self.mit_spikes=mit_spikes
        self.gran_spikes=gran_spikes

    def syn_gid(self, src_gid, tgt_gid):
        """Get the synaptic gid given the target gid and the source gids."""
        i = 0
        if (src_gid < self.num_mit): # Target is granule
            i = (tgt_gid * self.num_mit + src_gid + 1 + 
                    self.num_mit + self.num_gran) * 100 + 1
        else: # Target is mitral
            i = (src_gid * self.num_mit + tgt_gid + 1 + 
                    self.num_mit + self.num_gran) * 100
        return i
                
    def syngid2srctgt(self, syn_gid):
        m2g = syn_gid % 2
        i = syn_gid - m2g
        i /= 100
        i -= self.num_mit
        i -= self.num_gran
        i -= 1
        mgid = i % self.num_mit
        i -= mgid
        ggid = i/self.num_mit - self.num_mit
        if m2g == 1:
            return (mgid, ggid)
        else:
            return (ggid, mgid)

    def read_file(self,fname):
        """
        Create the raw data dict object read from the spike file.
        
        :param fname: Name of the spike file to read.
        """
        try:
            self.raw=spiketrain.read_file(fname)
        except Exception:
            self.raw=None
            raise

    def get_mitral_spikes(self,r=None,ids=None):
        """
        Get the spikes of the mitral cells.
        
        :param r: Optional. If specified, it must be a tuple of length 2 to 
            specify the range of values to acquire. The first parameter is the 
            lower bound timestamp and the second parameter is the upper bound
            timestamp. Spikes that fall within this window will be returned.
            
        :param ids: Optional. If specified, this must be a list to specify which
            ids will be returned. If the spike train of a single cell is
            desired, it must be encapsulated in a list with square brackets, 
            for example [<idx>]. IDs are zero-based.
        
        :returns: A list of lists containing the spike times of the mitral 
            cells.
        """
        
        # The spikes from the complete simulation are in the raw dict. If
        # the mit_spikes have not been copied from that object, do so now, and
        # store them in RAM.
        if len(self.mit_spikes)==0:
            # Try and get the spikes from the raw dict
            self.mit_spikes=[]   # Mitral cell spikes
            
            for i in range(self.num_mit):
                # The ids stored in the raw data dict are gid based.
                # Here, we are looking for numbers 0 to num_mit - 1.
                self.mit_spikes.append(spiketrain.get_spikes(self.raw,i))
            
            max_time = max(max(self.mit_spikes))
            if max_time>self.sim_time:
                self.sim_time=max_time
                
        if ids is None:
            ids=range(self.num_mit)

        if type(r) is tuple and len(r)==2:
            sub=[]
            for i in ids:
                sub.append([])
                # Use binary search for the spikes
                idxmin=bisect.bisect_left(self.mit_spikes[i],r[0])
                idxmax=bisect.bisect_right(self.mit_spikes[i],r[1])
                sub[i]=self.mit_spikes[i][idxmin:idxmax]
                
            return sub
        else:
            sub=[]
            idx = 0
            for i in ids:
                sub.append([])
                sub[idx]=self.mit_spikes[i]
                idx += 1
                
            return sub
                
        return self.mit_spikes
            
    def get_mitral_frequency(self,r=None,ids=None):
        """
        Get the spiking frequency of the mitral cells.
        
        :param r: Optional. If specified, it must be a tuple of length 2 to 
            specify the range of values to acquire. The first parameter is the 
            lower bound timestamp and the second parameter is the upper bound
            timestamp. Spikes that fall within this window will be returned.
            
        :param ids: Optional. If specified, it will return the spikes of the 
            ids specified in the list. If the spike train of a single cell is
            desired, it must be encapsulated in a list. IDs are zero-based.
            
        :returns: A vector of frequencies of each requested mitral cell.
        """
        # Simply divide the length of each spike vector by the simulation time
        # or window time.
        spikes=[] # The local copy of spikes to get the frequency on
        try:
            spikes=self.get_mitral_spikes(r,ids)
            if len(spikes)==0:
                    raise Exception('Unable to get mitral cell spikes')
        except Exception:
            raise
        freqm=numpy.zeros(self.num_mit)
        timems=self.sim_time
        if type(r) is tuple and len(r)==2:
            timems=float(r[1]-r[0])
        for i in range(len(spikes)):
            freqm[i]=float(len(spikes[i]))*1000./timems
        return freqm
    
    def get_granule_spikes(self,r=None,ids=None):
        """
        Get the spikes of the granule cells.
        
        :param r: Optional. If specified, it must be a list of length 2 to 
            specify the range of values to acquire. The first parameter is the 
            lower bound timestamp and the second parameter is the upper bound
            timestamp. Spikes that fall within this window will be returned.
            
        :param ids: Optional. If specified, it will return the spikes of the 
            ids specified in the list. If the spike train of a single cell is
            desired, it must be encapsulated in a list. IDs are zero-based.
        
        :returns: A list of lists containing the spike times of the granule 
            cells.
        """
        
        # The spikes from the complete simulation are in the raw dict object. If
        # the gran_spikes have not been copied from that object, do so now.
        if len(self.gran_spikes)==0:
            # Try and get the spikes from the raw dict
            self.gran_spikes=[]   # Granule cell spikes
            
            for i in range(self.num_mit,self.num_mit+self.num_gran):
                # The ids stored in the raw data dict are gid based.
                # Here, we are looking for numbers num_mit to 
                # num_mit + num_gran - 1.
                self.gran_spikes.append(spiketrain.get_spikes(self.raw,i))
        if type(r) is tuple and len(r)==2:
            sub=[]
            for i in range(self.num_gran):
                sub.append([])
                # Use binary search for the spikes
                idxmin=bisect.bisect_left(self.gran_spikes[i],r[0])
                idxmax=bisect.bisect_right(self.gran_spikes[i],r[1])
                sub[i]=self.gran_spikes[i][idxmin:idxmax]
            return sub
        return self.gran_spikes
                
    def get_granule_frequency(self,r=None,ids=None):
        """
        Get the spiking frequency of the granule cells.
        
        :param r: Optional. If specified, it must be a tuple of length 2 to 
            specify the range of values to acquire. The first parameter is the 
            lower bound timestamp and the second parameter is the upper bound
            timestamp. Spikes that fall within this window will be returned.
            
        :param ids: Optional. If specified, it will return the spikes of the 
            ids specified in the list. If the spike train of a single cell is
            desired, it must be encapsulated in a list. IDs are zero-based.
            
        :returns: A vector of frequencies of each requested granule cell.
        """
        # Simply divide the length of each spike vector by the simulation time
        # or window time.
        spikes=[] # The local copy of spikes to get the frequency on
        try:
            spikes=self.get_granule_spikes(r,ids)
            if len(spikes)==0:
                    raise Exception('Unable to get mitral cell spikes')
        except Exception:
            raise
        freqg=numpy.zeros(self.num_gran)
        timems=self.sim_time
        if type(r) is tuple and len(r)==2:
            timems=float(r[1]-r[0])
        for i in range(len(spikes)):
            freqg[i]=float(len(spikes[i]))*1000./timems
        return freqg
        
    def get_syn_spikes(self, src_gid, tgt_gid):
        return spiketrain.get_spikes(self.raw,self.syn_gid(src_gid, tgt_gid))

# If executed from the command line, read the file.
#if __name__ == '__main__':
#    read_file(sys.argv[1])
