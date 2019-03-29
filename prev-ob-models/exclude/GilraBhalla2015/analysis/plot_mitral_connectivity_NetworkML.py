# -*- coding: utf-8 -*-
from xml.etree import ElementTree as ET
import string
import sys
import operator # has itemgetter() for sorted()

sys.path.extend(["..","../neuroml","../cells","../channels","../networks","../simulations"])
from neuroml_utils import *
from sim_utils import * # has build_tweaks()
from networkConstants import *
## these are needed to load mitral cell and get segids and synapse locations
from load_channels import *
from load_synapses import *
from load_cells import *
from simset_odor import * # has synchan_activation_correction
load_channels()
## synchan_activation_correction depends on SIMDT,
## hence passed all the way down to
## granule_mitral_GABA synapse from the top level
load_synapses(synchan_activation_correction)
## cellSegmentDict[<cellname>] = { segid1 : [ segname,(proximalx,proximaly,proximalz),
##     (distalx,distaly,distalz),diameter,length,[potential_syn1, ... ] ] , ... }
cellSegmentDict = load_cells()

from pylab import *

## USAGE:
## python2.6 plot_mitral_connectivity_NetworkML.py <networkmlfile.xml>
## OR
## python2.6 plot_mitral_connectivity_NetworkML.py SYN_DECAY
## The plots are for weights and numbers of all exc/inh synapses on the central/lateral mitrals.
## The printing on the terminal gives weights located proximally vs distally from central mits.

TWOGLOMS = False#True # whether to check for mits 0,2 (True) or 0,1 (False)
INCLUDEMULTIS = True # whether to include multi-s in joint connectivity

class NetworkML():

    def __init__(self,numgloms,Exc_notInh=False,CentralWts_notLateral=True):
        self.numgloms=numgloms
        if TWOGLOMS: self.mitB=2
        else: self.mitB=1

        ## Whether exc or inh part of the reciprocal synapses should be considered
        ## False by default, as interested in inhibition onto central 'premits'
        self.Exc_notInh = Exc_notInh
        if Exc_notInh: self.weight_str = 'exc'
        else: self.weight_str = 'inh'
        ## Take numbers and weights of all connections
        ## from/to central mits (premits) VS from/to lateral mits (postmits)
        self.CentralWts_notLateral = CentralWts_notLateral
        if CentralWts_notLateral: self.centlat_str = 'central'
        else: self.centlat_str = 'lateral'

        print "READ CAREFULLY: Using",self.centlat_str,"mitrals'",self.weight_str,"weights and numbers."


    def readNetworkMLFromFile(self,filename,params={}):
        print "reading file ... ", filename
        tree = ET.parse(filename)
        root_element = tree.getroot()
        #print "Tweaking model ... "
        tweak_model(root_element, params)
        #print "Parsing model for mitrals connectivity via joint granules ... "
        return self.readNetworkML(root_element,root_element.attrib['lengthUnits'])

    def readNetworkML(self,network,lengthUnits="micrometer"):
        if lengthUnits in ['micrometer','micron']:
            self.length_factor = 1e-6
        else:
            self.length_factor = 1.0
        self.network = network
        #print "loading mitral positions ... "
        self.loadMitrals() # create cells
        #print "figuring joints connections ... "
        self.figureJoints() # create list of connections

    def loadMitrals(self):
        for population in self.network.findall(".//{"+nml_ns+"}population"):
            cellname = population.attrib["cell_type"]
            populationname = population.attrib["name"]
            if populationname != 'mitrals': continue
            self.mitralsDict = {}
            for instance in population.findall(".//{"+nml_ns+"}instance"):
                instanceid = int(instance.attrib['id'])
                location = instance.find('./{'+nml_ns+'}location')
                x = float(location.attrib['x'])*self.length_factor
                y = float(location.attrib['y'])*self.length_factor
                z = float(location.attrib['z'])*self.length_factor
                self.mitralsDict[instanceid]=(x,y,z)
        #print self.mitralsDict

    def figureJoints(self):
        ## make a list of singles and joints connected to each mitral
        self.singlesDict={}
        self.jointsDict={}
        projections = self.network.find(".//{"+nml_ns+"}projections")
        for projection in self.network.findall(".//{"+nml_ns+"}projection"):
            projectionname = projection.attrib["name"]
            ## connections listed in joints and multis are used to get
            ## the mitral cells connected to each other
            ## Whether exc or inh part of the reciprocal synapses should be considered
            if self.Exc_notInh:
                proj_joints = ( projectionname == 'mitral_granule_main_exc_joints' )
                proj_multis = (INCLUDEMULTIS and projectionname == 'mitral_granule_main_exc_multis')
                proj_singles = ( projectionname == 'mitral_granule_main_exc_singles' )
                mitpositionstr = 'pre'
                granpositionstr = 'post'
            else:
                proj_joints = ( projectionname == 'granule_mitral_inh_joints' )
                proj_multis = (INCLUDEMULTIS and projectionname == 'granule_mitral_inh_multis')
                proj_singles = ( projectionname == 'granule_mitral_inh_singles' )
                mitpositionstr = 'post'
                granpositionstr = 'pre'
            if proj_joints or proj_multis:
                for connection in projection.findall(".//{"+nml_ns+"}connection"):
                    mit_cell_id = int(connection.attrib[ mitpositionstr+'_cell_id' ])
                    gran_cell_id = int(connection.attrib[ granpositionstr+'_cell_id' ])
                    ## joint-s and multi-s granule id-s both start with 0.
                    ## make multi-id-s go negative, starting from -1, to maintain uniqueness.
                    if proj_multis: gran_cell_id = -gran_cell_id-1
                    mit_seg_id = int(connection.attrib[ mitpositionstr+'_segment_id' ])
                    gran_seg_id = int(connection.attrib[ granpositionstr+'_segment_id' ])
                    props = connection.find(".//{"+nml_ns+"}properties")
                    weight = float(props.attrib['weight'])
                    ## pre cells are granules, post are mitrals
                    ## joints will have two <connection> tags in the projection
                    ## multis will have multiple <connection> tags in the projection
                    if mit_cell_id in self.jointsDict:
                        self.jointsDict[mit_cell_id].append((gran_cell_id,mit_seg_id,gran_seg_id,weight))
                    else:
                        self.jointsDict[mit_cell_id] = [(gran_cell_id,mit_seg_id,gran_seg_id,weight)]
            elif proj_singles:
                for connection in projection.findall(".//{"+nml_ns+"}connection"):
                    mit_cell_id = int(connection.attrib[ mitpositionstr+'_cell_id' ])
                    gran_cell_id = int(connection.attrib[ granpositionstr+'_cell_id' ])
                    mit_seg_id = int(connection.attrib[ mitpositionstr+'_segment_id' ])
                    gran_seg_id = int(connection.attrib[ granpositionstr+'_segment_id' ])
                    props = connection.find(".//{"+nml_ns+"}properties")
                    weight = float(props.attrib['weight'])
                    if mit_cell_id in self.singlesDict:
                        self.singlesDict[mit_cell_id].append((gran_cell_id,mit_seg_id,gran_seg_id,weight))
                    else:
                        self.singlesDict[mit_cell_id] = [(gran_cell_id,mit_seg_id,gran_seg_id,weight)]
        
    def calc_connections(self):
        ## collate those joints that are connected to mitrals 0 or 1
        ## Also calculate distance between mitrals
        ## also figure out how many singles and joints
        ## are connected to primary vs secondary dendrites
        self.connectivityDict = {0:{},self.mitB:{}}
        self.connectivityDictNums = {0:{},self.mitB:{}}
        mit_prim_joints = {0:0,self.mitB:0}
        mit_somadend_joints = {0:0,self.mitB:0}
        mit_sec_joints = {0:0,self.mitB:0}
        mit_prim_singles = {0:0,self.mitB:0}
        mit_sec_singles = {0:0,self.mitB:0}
        if TWOGLOMS: postmitrange = [0,self.mitB]
        else: postmitrange = range(self.numgloms*MIT_SISTERS)
        for postmitid in postmitrange:
            self.connectivityDict[0][postmitid] = [0,0] # distance, inh weights
            self.connectivityDict[self.mitB][postmitid] = [0,0] # distance, inh weights
            self.connectivityDictNums[0][postmitid] = [0,0] # distance, number of grans
            self.connectivityDictNums[self.mitB][postmitid] = [0,0] # distance, number of grans
            ### print list of segments at which each mitral is connected to joints/multis.
            #gran2mainmitids,presegids,postsegids,_ = zip(*self.jointsDict[postmitid])
            #print "Mitral",postmitid,"is connected to joints on its segments :",presegids,\
            #    "total =",len(presegids)

        ## print list of segments at which mit2/3 has joints/multis. 3 is directed, 2 in not.
        ## Be very careful in finding segid of directedness if frac_directed is small like 0.3%.
        grans2premit,presegids,postsegids,_ = zip(*self.jointsDict[3])
        presegids = list(presegids)
        presegids.sort()
        print "Connections to joints/multis of mitral",3,"are at segments",presegids,\
            "total =",len(presegids)
        print "Be very careful in finding segid of directedness if frac_directed < 1%."
        
        for premitid in [0,self.mitB]:
            #### collate info about joints to mitrals 0 and self.mitB
            if premitid in self.jointsDict:
                ## separate out the granule_id-s to this premit, pre_seg_id-s and post_seg_id-s into separate lists
                grans2premit,mitsegids,gransegids,premit_weights = zip(*self.jointsDict[premitid])
                #print "Connections to joints of mitral",premitid,"are at segments",presegids,\
                #    "total =",len(presegids)
                #### count total number of joints connected to the two main mits
                #### segregate into those connected to primary dendrite vs soma+neardend vs sec dend
                #### if TWOGLOMS, count only those that connect to the other main mit.
                if premitid==0: postmitid=self.mitB
                else: postmitid=0
                for i,mitsegid in enumerate(mitsegids):
                    ## if TWOGLOMS, do not count joints to other mitrals, only to the two main ones.
                    if TWOGLOMS:
                        gran_id = grans2premit[i]
                        grans2postmit = zip(*self.jointsDict[postmitid])[0]
                        if gran_id not in grans2postmit:
                            continue
                    ## far prim dend
                    far_prim_dend_segids = [16,17,18,19,20]
                    ## soma & nearest-prim and near-sec dends
                    ## 0 is soma, 15 is the nearest prim dend, rest are near sec dends
                    close_dend_segids = list( set(proximal_mitral_segids) - set(far_prim_dend_segids) )
                    if mitsegid in far_prim_dend_segids: mit_prim_joints[premitid]+=premit_weights[i]
                    elif mitsegid in close_dend_segids: mit_somadend_joints[premitid]+=premit_weights[i]
                    ## rest of the sec dend
                    else: mit_sec_joints[premitid]+=premit_weights[i]
                #### find number of joints from each of the main mits
                #### to each of the remaining mitrals noting distance
                premitposition = self.mitralsDict[premitid]
                for postmitid in postmitrange:
                    if postmitid == premitid: continue # find granules only between non-self!
                    (postx,posty,postz) = self.mitralsDict[postmitid]
                    distance = ((premitposition[0]-postx)**2 + (premitposition[1]-posty)**2 +\
                        (premitposition[2]-postz)**2)**0.5
                    self.connectivityDict[premitid][postmitid][0] = distance * 1e6 # microns!
                    self.connectivityDictNums[premitid][postmitid][0] = distance * 1e6 # microns!
                    if postmitid not in self.jointsDict: continue
                    grans2postmit,mitsegids,gransegids,postmit_weights = zip(*self.jointsDict[postmitid])
                    ## Take numbers and weights of connections from/to central mits (premits)
                    ## VS from/to lateral mits (postmits)
                    if self.CentralWts_notLateral:
                        grans2mitA = grans2premit
                        grans2mitB = grans2postmit
                        mit_weights = premit_weights
                    else:
                        grans2mitB = grans2premit
                        grans2mitA = grans2postmit
                        mit_weights = postmit_weights
                    for gran_id in set(grans2mitB): # can be many of the same gran_id-s, take unique set
                        if gran_id in grans2mitA:
                            granidxs = where(array(grans2mitA)==gran_id)[0] # multiple syns from granid to premit
                            for granidx in granidxs:
                                ## Be careful to take the premit weight! weight decays exponentially along dendrite!
                                self.connectivityDict[premitid][postmitid][1] += mit_weights[granidx]
                                self.connectivityDictNums[premitid][postmitid][1] += 1
            #### collate info about singles to mitrals 0 and self.mitB
            ## For singles, CentralWts_notLateral has no meaning, always central mit (premit) weights
            if premitid in self.singlesDict:
                ## separate out the gran_cell_id-s, mit_seg_id-s, gran_seg_id-s, weights into separate lists
                grans2premit,mitsegids,gransegids,premit_weights = zip(*self.singlesDict[premitid])
                ## count inh wt of the singles on to the primary dendrite vs sec dend
                for i,mitsegid in enumerate(mitsegids):
                    if mitsegid in [15,16,17,18,19,20]: mit_prim_singles[premitid]+=premit_weights[i]
                    else: mit_sec_singles[premitid]+=premit_weights[i]
        self.mit_posn_grannos = \
            (mit_prim_joints,mit_somadend_joints,mit_sec_joints,mit_prim_singles,mit_sec_singles)

        return (self.mitralsDict,self.connectivityDict,self.connectivityDictNums,self.mit_posn_grannos)

    def print_info(self):
        (mit_prim_joints,mit_somadend_joints,mit_sec_joints,mit_prim_singles,mit_sec_singles) = \
            self.mit_posn_grannos
        print self.weight_str,"central-mit\'s weight of joints at far primary dendrites of mitA and mitB is ",mit_prim_joints
        print self.weight_str,"central-mit\'s weight of joints at soma and near prim & sec dendrites of mitA and mitB is ",mit_somadend_joints
        print self.weight_str,"central-mit\'s weight of joints at secondary dendrites of mitA and mitB is ",mit_sec_joints
        print self.weight_str,"central-mit\'s weight of singles at primary dendrites of mitA and mitB is ",mit_prim_singles
        print self.weight_str,"central-mit\'s weight of singles at secondary dendrites of mitA and mitB is ",mit_sec_singles

    def calc_synaptic_decay(self,indistance,mitid=0):
        """ indistance is max distance of segment in m,
        indends is a list of lat/peripheral dendrites e.g. ['d1','p1',...]"""
        ## if lateral mitral, get weights for only the directed dendrite
        if mitid != 0:
            ## list of segments at which mit2/3 has joints/multis.
            count_dict = {}
            for _,mit_seg_id,_,_ in self.jointsDict[mitid]:
                if mit_seg_id not in count_dict: count_dict[mit_seg_id] = 0
                else: count_dict[mit_seg_id] += 1
            ## get the segment that has the most joints/multis, and is a lat dend
            while True:
                mit_seg_id_max = max(count_dict.iterkeys(), key=(lambda key: count_dict[key]))
                segname = cellSegmentDict['mitral'][str(mit_seg_id_max)][0]
                ## latdend of lateral mit that directs to mit0 soma
                ## it should one of d1-d4 or p1-p4
                ## typically a few p [proximal?] lat dend compts go on to d [distal?] compts
                dendnum = segname[14:15]
                if dendnum in ['1','2','3','4']:
                    indends = ['p'+dendnum,'d'+dendnum] # p# and d# belong to dend num #.
                    print 'Found compartment',segname,'with most syns.'
                    break
                else:
                    del count_dict[mit_seg_id_max] ## remove for next iteration
                    print 'Not using compartment',segname,'checking next most...'

        ## get the weights on each compartment
        prim_dist_weight_list = []
        sec_dist_weight_list = []
        for segment in cellSegmentDict['mitral'].values():
            ## segment = [ segname,(proximalx,proximaly,proximalz),\
            ##    (distalx,distaly,distalz),diameter,length,[potential_syn1, ... ] ]
            ## segment[0] is segment name, which has segid at the end after an underscore
            segname = segment[0]
            segid = int(string.split(segname,'_')[-1])
            ## use this seg only if granule_mitral is one of the potential synapses here
            if "granule_mitral" in segment[5]:
                ## average segment position from the proximal and distal points
                segx1 = segment[1][0]
                segx2 = segment[2][0]
                segx = ( segx1 + segx2 ) / 2.0
                segy1 = segment[1][1]
                segy2 = segment[2][1]
                segy = ( segy1 + segy2 ) / 2.0
                segz1 = segment[1][2]
                segz2 = segment[2][2]
                segz = ( segz1 + segz2 ) / 2.0
                ## segid is for the prototype cell whose soma x0,y0,z0 is at origin.
                ## The start of segment is taken in dend_decay()
                ## in generate_neuroml.py for exp decay,
                ## but here, we want the mid-segment for distance of synapses.
                distance = sqrt(segx**2+segy**2+segz**2)
                ## only consider distances up to indistance
                if distance>indistance: continue
                wt_list = []
                ## if lateral mitral, keep only lat directed dend or primary dendrite
                if mitid != 0:
                    dendname = segname[13:15] # take 2 chars from string secname
                    if dendname not in indends and segid not in [0,15,16,17,18,19,20]:
                        ## It is important to have empty wt_list for dendrites that are excluded,
                        ## since avg_synaptic_decay() assumes the same order of segments
                        ## across different network seeds.
                        sec_dist_weight_list.append((distance,segname,wt_list))
                        continue
                ## singly connected granules
                for (gran_cell_id,mit_seg_id,gran_seg_id,weight) in self.singlesDict[mitid]:
                    if mit_seg_id==segid:
                        if not self.Exc_notInh:
                            ## weight of each inh synapse on a singly-connected granule is actually
                            ## GRANS_CLUB_SINGLES/float(SYNS_PER_CLUBBED_SINGLE) times unit wt,
                            ## Each syn represents (GRANS_CLUB_SINGLES/SYNS_PER_CLUBBED_SINGLE) # of syns
                            ## wt_list is a list of individual synapse weights, undo clubbing below ...
                            wt_list.extend( \
                                [weight/(GRANS_CLUB_SINGLES/float(SYNS_PER_CLUBBED_SINGLE))] \
                                * (GRANS_CLUB_SINGLES/SYNS_PER_CLUBBED_SINGLE) )
                        else:
                            ## Exc weight from a given mitral to granules connected to it.
                            wt_list.append(weight)
                ## jointly and multiply connected granules
                for (gran_cell_id,mit_seg_id,gran_seg_id,weight) in self.jointsDict[mitid]:
                    if mit_seg_id==segid: wt_list.append(weight)
                if segid in [0,15,16,17,18,19,20]: # soma & prim dend segments
                    prim_dist_weight_list.append((distance,segname,wt_list))
                else: # sec dend segments
                    sec_dist_weight_list.append((distance,segname,wt_list))
        return prim_dist_weight_list,sec_dist_weight_list

latdistlim = 800 ## plot only lateral 800 microns

def avg_synaptic_decay(seednums,frac_dir_str,\
        Exc_notInh=False,CentralWts_notLateral=True,mitid=0):
    numgloms = 3
    allprimwtlists = None
    allsecwtlists = None
    for seednum in seednums:
        ### The _4syndecay named files below, have minimal synvariation=0.01 uniform,
        ### and allow_multi_granule_conn = False
        #filename = \
        #    '../netfiles/syn_conn_array_10000_singlesclubbed100_jointsclubbed1_numgloms'+\
        #    str(numgloms)+'_seed'+str(seednum)+'_directed'+frac_dir_str+'_proximal_4syndecay.xml'
        ## But Upi suggested I use my usual netfiles instead of _4syndecay above,
        ## and plot an average curve over an example.
        filename = \
            '../netfiles/syn_conn_array_10000_singlesclubbed100_jointsclubbed1_numgloms'+\
            str(numgloms)+'_seed'+str(seednum)+'_directed'+frac_dir_str+'_proximal.xml'
        netml = NetworkML(numgloms,Exc_notInh,CentralWts_notLateral)
        tweaks = {}
        netml.readNetworkMLFromFile(filename,params=tweaks)
        prim_dist_weight_list,sec_dist_weight_list = \
            netml.calc_synaptic_decay(latdistlim*1e-6,mitid) # microns to m
        ## zip(*  ) is the inverse of zip()
        ## * is the unpacking operator
        sortedprimdists,sortednames,sortedwtlists = zip(*sorted(prim_dist_weight_list))
        if allprimwtlists is None: allprimwtlists = sortedwtlists
        else:
            for i,wtlist in enumerate(sortedwtlists): allprimwtlists[i].extend(wtlist)
        ## sec_dist_weight_list = [(distance,segname,wt_list),...]
        ## sort by sec dend num and then distance, same sec dend stays together
        ## first 15 chars of segname contain dendrite num d1,d2,d3,d4. eg:Seg0_sec_dendd3_4_225
        sortedsecdists,sortednames,sortedwtlists = \
            zip(*sorted( sec_dist_weight_list,key=(lambda x: x[1][:15]+'%01.6f'%x[0]) ))
        if allsecwtlists is None: allsecwtlists = sortedwtlists
        else:
            for i,wtlist in enumerate(sortedwtlists): allsecwtlists[i].extend(wtlist)

    sortedprimwts = [ mean(wtlist) for wtlist in allprimwtlists ]
    sortedsecwts = [ mean(wtlist) for wtlist in allsecwtlists ]
    return sortedprimdists,sortedprimwts,sortedsecdists,sortedsecwts,sortednames

def bin_by_distance(dist_list,wt_list):
    ## add one bin at the end which only acts to set limits, but contains nothing
    distance_bins = array( \
        [15,50,80,120,200,270,350,420,485,575]+range(675,latdistlim+101,100) )*1e-6
    dist_numbins = len(distance_bins)
    binned_wts = [0.0]*len(distance_bins)
    binned_nums = [0.0]*len(distance_bins)
    for distnum,dist in enumerate(dist_list):
        wt = wt_list[distnum]
        if isnan(wt): continue
        ## small list, not bothering to do binary search, just sequential search
        bin_lower = 0
        i = 0
        while True:
            bin_upper = (distance_bins[i]+distance_bins[i+1])/2.0
            if dist>=bin_lower and dist<bin_upper:
                binned_wts[i] += wt
                binned_nums[i] += 1
                break
            else:
                bin_lower = bin_upper
                ## if dist is outside values in distance_bins, go to next dist
                if i>=dist_numbins: break
                i += 1
    ## do not return the last set-limit bin
    return distance_bins[0:-1],\
        array(binned_wts[0:-1])/array(binned_nums[0:-1]) # element-wise division

def plot_synaptic_decay_paperfigure():
    """The single seed scatter plot calculations
     and the directed only mitid=3 plot calculations
     may give some floating point warnings,
     as some of the compartments have no synapses, or wt_list is set to [],
     and so their weight values become nan-s on taking mean. """
    seednums = arange(750.0,760.0,1.0)
    eg_seednum_idx = 4
    latdisttick = 600 ## tick at lateral 600 microns, rest show with clip off

    ## inh weights on central mit as a function of distance on soma-prim and sec dend
    primdists,primwts,secdists,secwts,secnames = avg_synaptic_decay(seednums,'0.01')
    sec_distance_bins,binned_secwts = bin_by_distance(secdists,secwts)
    _,primwts_undir,_,secwts_undir,_ = avg_synaptic_decay(seednums,'0.0')
    _,binned_secwts_undir = bin_by_distance(secdists,secwts_undir)
    primdists_um = array(primdists)*1e6
    secdists_um = array(secdists)*1e6
    sec_distance_bins_um = array(sec_distance_bins)*1e6

    #### soma-primary dendrite inh gran-|mit plot
    fig = figure(figsize=(columnwidth*7/8.0,linfig_height/3.0),\
        dpi=fig_dpi,facecolor='none') # none = transparent
    ax = fig.add_subplot(111)
    ax.plot(primdists_um,array(primwts_undir)*granule_mitral_GABA_Gbar*1e9,\
        color='b',marker=',',ms=marker_size,linewidth=linewidth)
    ax.plot(primdists_um,array(primwts)*granule_mitral_GABA_Gbar*1e9,\
        color='r',marker=',',ms=marker_size,linewidth=linewidth)
    ## inh weights, directed, one network seed only
    _,primwts_undir,_,_,_ = avg_synaptic_decay([seednums[eg_seednum_idx]],'0.0')
    _,primwts,_,_,_ = avg_synaptic_decay([seednums[eg_seednum_idx]],'0.01')
    ax.scatter(primdists_um,array(primwts_undir)*granule_mitral_GABA_Gbar*1e9,\
        color='b',marker='s',s=marker_size)
    ax.scatter(primdists_um,array(primwts)*granule_mitral_GABA_Gbar*1e9,\
        color='r',marker='x',s=marker_size)
    beautify_plot(ax,xticksposn='bottom',yticksposn='left')
    #axes_labels(ax,'$\delta$ ($\mu m$)','$\sigma$ (nS)',fontsize=label_fontsize)
    add_scalebar(ax,matchx=False,matchy=False,hidex=False,hidey=False,\
        sizex=100,labelx='100 $\mu m$',sizey=1,labely='  1 nS',\
        bbox_to_anchor=[0.8,0.3],bbox_transform=ax.transAxes)
    fig_clip_off(fig)
    fig.tight_layout()
    fig.savefig('../figures/connectivity/soma-primdend-inh-vs-distance.svg',\
        dpi=fig.dpi,transparent=True)
    
    #### sec dend inh gran-|mit plot
    fig = figure(figsize=(columnwidth/2.0,linfig_height/3.0),\
        dpi=fig_dpi,facecolor='none') # none = transparent
    ax = fig.add_subplot(111)
    ax.plot(sec_distance_bins_um,array(binned_secwts_undir)*granule_mitral_GABA_Gbar*1e9,\
        color='b',marker=',',ms=marker_size,linewidth=linewidth)
    ax.plot(sec_distance_bins_um,array(binned_secwts)*granule_mitral_GABA_Gbar*1e9,\
        color='r',marker=',',ms=marker_size,linewidth=linewidth)
    _,_,_,ymax = beautify_plot(ax,xticksposn='bottom',yticksposn='left')
    ## inh weights, directed, one network seed only, all dendrites
    _,_,_,secwts_undir,_ = avg_synaptic_decay([seednums[eg_seednum_idx]],'0.0')
    _,_,_,secwts,_ = avg_synaptic_decay([seednums[eg_seednum_idx]],'0.01')
    ax.scatter(secdists_um,array(secwts_undir)*granule_mitral_GABA_Gbar*1e9,\
        color='b',marker='s',s=marker_size)
    ax.scatter(secdists_um,array(secwts)*granule_mitral_GABA_Gbar*1e9,\
        color='r',marker='x',s=marker_size)
    #axes_labels(ax,'$\delta$ ($\mu m$)','$\sigma$ (nS)',fontsize=label_fontsize)
    add_scalebar(ax,matchx=False,matchy=False,hidex=False,hidey=False,\
        sizex=250,labelx='250 $\mu m$',sizey=3,labely='  3 nS',\
        bbox_to_anchor=[0.9,0.3],bbox_transform=ax.transAxes)
    ax.set_xlim(0,latdisttick)
    ax.set_ylim(0,8)
    ax.set_yticks([0,8])
    fig.tight_layout()
    fig_clip_off(fig)
    fig.savefig('../figures/connectivity/secdend-inh-vs-distance.svg',\
        dpi=fig.dpi,transparent=True)

    #### sec dend exc mit->gran plot for central and lateral mitral
    fig = figure(figsize=(columnwidth/2.0,linfig_height*1.15),\
        dpi=fig_dpi,facecolor='none') # none = transparent
    ax = fig.add_subplot(3,1,1)
    ## exc weights of central mit as a function of distance on sec dend
    ## Not extra-directed, avg over all networks, all dendrites
    _,primwts,_,secwts,_ = avg_synaptic_decay(seednums,'0.0',Exc_notInh=True)
    _,binned_wts = bin_by_distance(secdists,secwts)
    ax.plot(sec_distance_bins_um,array(binned_wts)*mitral_granule_AMPA_Gbar*1e9,\
        color='b',linestyle='solid',marker=',',linewidth=linewidth,ms=marker_size)
    ## extra-directed, avg over all networks, all dendrites
    _,primwts,_,secwts,_ = avg_synaptic_decay(seednums,'0.01',Exc_notInh=True)
    distance_bins,binned_wts = bin_by_distance(secdists,secwts)
    ax.plot(sec_distance_bins_um,array(binned_wts)*mitral_granule_AMPA_Gbar*1e9,\
        color='r',linestyle='solid',marker=',',linewidth=linewidth,ms=marker_size)
    ## extra-directed, one network, all dendrites (lateral and peripheral)
    _,primwts,_,secwts_undir,_ = \
        avg_synaptic_decay([seednums[eg_seednum_idx]],'0.0',Exc_notInh=True)
    _,primwts,_,secwts,_ = \
        avg_synaptic_decay([seednums[eg_seednum_idx]],'0.01',Exc_notInh=True)
    ax.scatter(secdists_um,array(secwts_undir)*mitral_granule_AMPA_Gbar*1e9,color='b',\
        marker='s',s=marker_size)
    ax.scatter(secdists_um,array(secwts)*mitral_granule_AMPA_Gbar*1e9,color='r',\
        marker='x',s=marker_size)
    _,_,_,ymax = beautify_plot(ax,xticksposn='bottom',yticksposn='left')
    #axes_labels(ax,'$\delta$ ($\mu m$)','$\sigma$ (nS)',fontsize=label_fontsize)
    add_scalebar(ax,matchx=False,matchy=False,hidex=False,hidey=False,\
        sizex=250,labelx='250 $\mu m$',sizey=0.2,labely='0.2 nS',\
        bbox_to_anchor=[1.1,-1.0],bbox_transform=ax.transAxes)
    ax.set_xlim(0,latdisttick)
    ax.set_ylim(0,0.6)
    ax.set_xticks([0])
    ax.set_yticks([0,0.6])

    ax = fig.add_subplot(3,1,3)
    ## exc weights of lateral mit as a function of distance on its sec dend
    ## mitid=2 does the trick, CentralWts_notLateral is for later processing not needed here.
    ## get secdists again as now only directed lat dend is averaged over.
    _,primwts,secdists,secwts,_ = avg_synaptic_decay(seednums,'0.0',Exc_notInh=True,\
        CentralWts_notLateral=False,mitid=2)
    secdists_um = array(secdists)*1e6
    sec_distance_bins,binned_wts = bin_by_distance(secdists,secwts)
    sec_distance_bins_um = array(sec_distance_bins)*1e6
    ax.plot(-sec_distance_bins_um,array(binned_wts)*mitral_granule_AMPA_Gbar*1e9,\
        color='b',marker=',',linewidth=linewidth)
    _,primwts,_,secwts,_ = avg_synaptic_decay(seednums,'0.01',Exc_notInh=True,\
        CentralWts_notLateral=False,mitid=2)
    _,binned_wts = bin_by_distance(secdists,secwts)
    ax.plot(-sec_distance_bins_um,array(binned_wts)*mitral_granule_AMPA_Gbar*1e9,\
        color='g',marker=',',linewidth=linewidth)
    _,_,_,ymax = beautify_plot(ax,xticksposn='bottom',yticksposn='right',drawyaxis=False)
    ## extra-directed, one network, only one directed dendrite
    _,primwts,secdists,secwts_undir,_ = avg_synaptic_decay([seednums[eg_seednum_idx]],'0.0',\
        Exc_notInh=True,CentralWts_notLateral=False,mitid=2)
    secdists_um = array(secdists)*1e6
    ax.scatter(-secdists_um,array(secwts_undir)*mitral_granule_AMPA_Gbar*1e9,color='b',\
        marker='s',s=marker_size)
    # secdists for single seeds can be different, since compartments with no syns are ignored
    _,primwts,secdists,secwts,_ = avg_synaptic_decay([seednums[eg_seednum_idx]],'0.01',\
        Exc_notInh=True,CentralWts_notLateral=False,mitid=2)
    secdists_um = array(secdists)*1e6
    ax.scatter(-secdists_um,array(secwts)*mitral_granule_AMPA_Gbar*1e9,color='g',\
        marker='x',s=marker_size)
    #axes_labels(ax,'$\delta$ ($\mu m$)','$\sigma$ (nS)',fontsize=label_fontsize)
    ax.spines['right'].set_color('k') # draw right axis
    ## xlim and ylim must be same as for above axes, since scale bar is shared 
    ax.set_xlim(-latdisttick,0)
    ax.set_ylim(0,0.6)
    ax.set_xticks([0])
    ax.set_yticks([0,0.6])
    fig.tight_layout()
    fig_clip_off(fig)
    fig.savefig('../figures/connectivity/secdend-exc-vs-distance.svg',\
        dpi=fig.dpi,transparent=True)

    show()

def plot_synaptic_decay_paperfigure_v2():
    """The single seed scatter plot calculations
     and the directed only mitid=3 plot calculations
     may give some floating point warnings,
     as some of the compartments have no synapses, or wt_list is set to [],
     and so their weight values become nan-s on taking mean. """
    seednums = arange(750.0,760.0,1.0)
    eg_seednum_idx = 4
    latdisttick = 800 ## tick at lateral 600 microns, rest show with clip off

    ## inh weights on central mit as a function of distance on soma-prim and sec dend
    primdists,primwts,secdists,secwts,secnames = avg_synaptic_decay(seednums,'0.01')
    sec_distance_bins,binned_secwts = bin_by_distance(secdists,secwts)
    _,primwts_undir,_,secwts_undir,_ = avg_synaptic_decay(seednums,'0.0')
    _,binned_secwts_undir = bin_by_distance(secdists,secwts_undir)
    primdists_um = array(primdists)*1e6
    secdists_um = array(secdists)*1e6
    sec_distance_bins_um = array(sec_distance_bins)*1e6

    #### soma-primary dendrite inh gran-|mit plot
    fig = figure(figsize=(columnwidth/2.0,linfig_height/3.0*1.2),\
        dpi=fig_dpi,facecolor='w') # none = transparent
    ax = fig.add_subplot(111)
    ax.plot(primdists_um,array(primwts_undir)*granule_mitral_GABA_Gbar*1e9,\
        color='b',marker=',',ms=marker_size,linewidth=linewidth)
    ax.plot(primdists_um,array(primwts)*granule_mitral_GABA_Gbar*1e9,\
        color='r',marker=',',ms=marker_size,linewidth=linewidth)
    ## inh weights, directed, one network seed only
    _,primwts_undir,_,_,_ = avg_synaptic_decay([seednums[eg_seednum_idx]],'0.0')
    _,primwts,_,_,_ = avg_synaptic_decay([seednums[eg_seednum_idx]],'0.01')
    ax.scatter(primdists_um,array(primwts_undir)*granule_mitral_GABA_Gbar*1e9,\
        color='b',marker='s',s=marker_size)
    ax.scatter(primdists_um,array(primwts)*granule_mitral_GABA_Gbar*1e9,\
        color='r',marker='x',s=marker_size)
    beautify_plot(ax,xticksposn='bottom',yticksposn='left')
    ax.set_xlim(0,550)
    ax.set_xticks([0,550])
    axes_labels(ax,'distance ($\mu m$)',u'G─┤M (nS)',fontsize=label_fontsize,xpad=-4,ypad=2)
    #add_scalebar(ax,matchx=False,matchy=False,hidex=False,hidey=False,\
    #    sizex=100,labelx='100 $\mu m$',sizey=1,labely='  1 nS',\
    #    bbox_to_anchor=[0.8,0.3],bbox_transform=ax.transAxes)
    fig_clip_off(fig)
    fig.tight_layout()
    fig.savefig('../figures/connectivity/soma-primdend-inh-vs-distance_v2.svg',\
        dpi=fig.dpi,transparent=True)

    fig = figure(figsize=(columnwidth/2.0,linfig_height),\
        dpi=fig_dpi,facecolor='w') # none = transparent
    #### sec dend inh gran-|mit plot
    ax = fig.add_subplot(3,1,1)
    ax.plot(sec_distance_bins_um,array(binned_secwts_undir)*granule_mitral_GABA_Gbar*1e9,\
        color='b',marker=',',ms=marker_size,linewidth=linewidth)
    ax.plot(sec_distance_bins_um,array(binned_secwts)*granule_mitral_GABA_Gbar*1e9,\
        color='r',marker=',',ms=marker_size,linewidth=linewidth)
    _,_,_,ymax = beautify_plot(ax,xticksposn='bottom',yticksposn='left')
    ## inh weights, directed, one network seed only, all dendrites
    _,_,_,secwts_undir,_ = avg_synaptic_decay([seednums[eg_seednum_idx]],'0.0')
    _,_,_,secwts,_ = avg_synaptic_decay([seednums[eg_seednum_idx]],'0.01')
    ax.scatter(secdists_um,array(secwts_undir)*granule_mitral_GABA_Gbar*1e9,\
        color='b',marker='s',s=marker_size)
    ax.scatter(secdists_um,array(secwts)*granule_mitral_GABA_Gbar*1e9,\
        color='r',marker='x',s=marker_size)
    axes_labels(ax,'',u'G─┤M (nS)',fontsize=label_fontsize,ypad=9)
    #add_scalebar(ax,matchx=False,matchy=False,hidex=False,hidey=False,\
    #    sizex=250,labelx='250 $\mu m$',sizey=3,labely='  3 nS',\
    #    bbox_to_anchor=[0.9,0.3],bbox_transform=ax.transAxes)
    ax.set_xlim(0,latdisttick)
    ax.set_xticklabels(['',''])
    ax.set_ylim(0,8)
    ax.set_yticks([0,8])

    #### sec dend exc mit->gran plot for central and lateral mitral
    ax = fig.add_subplot(3,1,2)
    ## exc weights of central mit as a function of distance on sec dend
    ## Not extra-directed, avg over all networks, all dendrites
    _,primwts,_,secwts,_ = avg_synaptic_decay(seednums,'0.0',Exc_notInh=True)
    _,binned_wts = bin_by_distance(secdists,secwts)
    ax.plot(sec_distance_bins_um,array(binned_wts)*mitral_granule_AMPA_Gbar*1e9,\
        color='b',linestyle='solid',marker=',',linewidth=linewidth,ms=marker_size)
    ## extra-directed, avg over all networks, all dendrites
    _,primwts,_,secwts,_ = avg_synaptic_decay(seednums,'0.01',Exc_notInh=True)
    distance_bins,binned_wts = bin_by_distance(secdists,secwts)
    ax.plot(sec_distance_bins_um,array(binned_wts)*mitral_granule_AMPA_Gbar*1e9,\
        color='r',linestyle='solid',marker=',',linewidth=linewidth,ms=marker_size)
    ## extra-directed, one network, all dendrites (lateral and peripheral)
    _,primwts,_,secwts_undir,_ = \
        avg_synaptic_decay([seednums[eg_seednum_idx]],'0.0',Exc_notInh=True)
    _,primwts,_,secwts,_ = \
        avg_synaptic_decay([seednums[eg_seednum_idx]],'0.01',Exc_notInh=True)
    ax.scatter(secdists_um,array(secwts_undir)*mitral_granule_AMPA_Gbar*1e9,color='b',\
        marker='s',s=marker_size)
    ax.scatter(secdists_um,array(secwts)*mitral_granule_AMPA_Gbar*1e9,color='r',\
        marker='x',s=marker_size)
    _,_,_,ymax = beautify_plot(ax,xticksposn='bottom',yticksposn='left')
    axes_labels(ax,'',u'M→G (nS)',fontsize=label_fontsize,ypad=2)
    #add_scalebar(ax,matchx=False,matchy=False,hidex=False,hidey=False,\
    #    sizex=250,labelx='250 $\mu m$',sizey=0.2,labely='0.2 nS',\
    #    bbox_to_anchor=[1.1,-1.0],bbox_transform=ax.transAxes)
    ax.set_xlim(0,latdisttick)
    ax.set_xticks([0,latdisttick])
    ax.set_xticklabels(['',''])
    ax.set_ylim(0,0.6)
    ax.set_yticks([0,0.6])
    ax.set_yticklabels(['0','0.6'])

    ax = fig.add_subplot(3,1,3)
    ## exc weights of lateral mit as a function of distance on its sec dend
    ## mitid=2 does the trick, CentralWts_notLateral is for later processing not needed here.
    ## get secdists again as now only directed lat dend is averaged over.
    _,primwts,secdists,secwts,_ = avg_synaptic_decay(seednums,'0.0',Exc_notInh=True,\
        CentralWts_notLateral=False,mitid=2)
    secdists_um = array(secdists)*1e6
    sec_distance_bins,binned_wts = bin_by_distance(secdists,secwts)
    sec_distance_bins_um = array(sec_distance_bins)*1e6
    ax.plot(sec_distance_bins_um,array(binned_wts)*mitral_granule_AMPA_Gbar*1e9,\
        color='b',marker=',',linewidth=linewidth)
    _,primwts,_,secwts,_ = avg_synaptic_decay(seednums,'0.01',Exc_notInh=True,\
        CentralWts_notLateral=False,mitid=2)
    _,binned_wts = bin_by_distance(secdists,secwts)
    ax.plot(sec_distance_bins_um,array(binned_wts)*mitral_granule_AMPA_Gbar*1e9,\
        color='g',marker=',',linewidth=linewidth)
    _,_,_,ymax = beautify_plot(ax,xticksposn='bottom',yticksposn='left')
    ## extra-directed, one network, only one directed dendrite
    _,primwts,secdists,secwts_undir,_ = avg_synaptic_decay([seednums[eg_seednum_idx]],'0.0',\
        Exc_notInh=True,CentralWts_notLateral=False,mitid=2)
    secdists_um = array(secdists)*1e6
    ax.scatter(secdists_um,array(secwts_undir)*mitral_granule_AMPA_Gbar*1e9,color='b',\
        marker='s',s=marker_size)
    # secdists for single seeds can be different, since compartments with no syns are ignored
    _,primwts,secdists,secwts,_ = avg_synaptic_decay([seednums[eg_seednum_idx]],'0.01',\
        Exc_notInh=True,CentralWts_notLateral=False,mitid=2)
    secdists_um = array(secdists)*1e6
    ax.scatter(secdists_um,array(secwts)*mitral_granule_AMPA_Gbar*1e9,color='g',\
        marker='x',s=marker_size)
    axes_labels(ax,'distance ($\mu m$)',u'M→G (nS)',fontsize=label_fontsize,xpad=-4,ypad=2)
    ## xlim and ylim must be same as for above axes, since scale bar is shared 
    ax.set_xlim(0,latdisttick)
    ax.set_ylim(0,0.6)
    ax.set_xticks([0,latdistlim])
    ax.set_yticks([0,0.6])
    ax.set_yticklabels(['0','0.6'])

    fig.tight_layout()
    fig_clip_off(fig)
    fig.savefig('../figures/connectivity/secdend-vs-distance_v2.svg',\
        dpi=fig.dpi,transparent=True)

    show()

if __name__ == "__main__":
    err_string = "Please give me a NeuroML filename or SYN_DECAY as argument."
    if len(sys.argv) < 2:
        print err_string
        sys.exit(1)
    filename = sys.argv[1]

    if 'SYN_DECAY' in sys.argv:
        #plot_synaptic_decay_paperfigure()
        plot_synaptic_decay_paperfigure_v2()

    else:
        if TWOGLOMS: numgloms=2
        else:
            postnumgloms_str = filename.split('numgloms')[1] # take filename after 'numgloms'
            numgloms = int(postnumgloms_str.split('_')[0]) # take the integer between 'numgloms' and '_'
        netml = NetworkML(numgloms)
        #if TWOGLOMS:
        #    includeProjections = ['granule_baseline']
        #    tweaks = build_tweaks( mitralsclub=True, nospineinh=False, nosingles=False,
        #        nojoints=False, nomultis=False, nopgs=False, onlytwomits=True,
        #        includeProjections=includeProjections, twomitrals=(0,netml.mitB) )
        #else: tweaks = {}
        tweaks = {}
        netml.readNetworkMLFromFile(filename,params=tweaks)
        (mitralsDict,connectivityDict,connectivityDictNums,mit_posn_grannos) = netml.calc_connections()
        netml.print_info()
        ## dict.iteritems() returns a list of (key,value) pairs.
        ## sorted sorts above list after calling itemgetter(1)
        ## on each item which returns the second element i.e. value
        ## finally, a list of (key, value) pairs is returned sorted by value
        mitAconnectivityList = sorted(connectivityDict[0].iteritems(), key=operator.itemgetter(1))
        mitBconnectivityList = sorted(connectivityDict[netml.mitB].iteritems(), key=operator.itemgetter(1))
        mitAconnectivityNumsList = sorted(connectivityDictNums[0].iteritems(), key=operator.itemgetter(1))
        mitBconnectivityNumsList = sorted(connectivityDictNums[netml.mitB].iteritems(), key=operator.itemgetter(1))

        fig = figure()
        ax = fig.add_subplot(111)
        ## zip(*  ) is the inverse of zip()
        ## * is the unpacking operator
        ## key is postmitid; value is (distance,weights,#joints)
        sortedkeys,sortedvalues = zip(*mitAconnectivityList)
        for key in sortedkeys:
            ax.annotate(str(key), xy=connectivityDict[0][key])
        x,y = zip(*sortedvalues)
        ax.plot(x,y,'r-o',label='mitral 0')
        ## zip(*  ) is the opposite of zip() 
        ## key is postmitid; value is (distance,#joints)
        sortedkeys,sortedvalues = zip(*mitBconnectivityList)
        for key in sortedkeys:
            ax.annotate(str(key), xy=connectivityDict[netml.mitB][key])
        x,y = zip(*sortedvalues)
        ax.plot(x,y,'g-x',label='mitral '+str(netml.mitB))
        legend(loc='upper left')
        xlabel('distance (microns)')
        ylabel(netml.centlat_str+'-mit\'s '+netml.weight_str+' weight at joint+multi granules')
        title(netml.centlat_str+'-mit\'s '+netml.weight_str+' weight at all joints/multis between mits 0/1 <--> mit x')

        fig = figure()
        ax = fig.add_subplot(111)
        sortedkeys,sortedvalues = zip(*mitAconnectivityNumsList)
        for key in sortedkeys:
            ax.annotate(str(key), xy=connectivityDictNums[0][key])
        x,y = zip(*sortedvalues)
        ax.plot(x,y,'r-o',label='mitral 0')
        sortedkeys,sortedvalues = zip(*mitBconnectivityNumsList)
        for key in sortedkeys:
            ax.annotate(str(key), xy=connectivityDictNums[netml.mitB][key])
        x,y = zip(*sortedvalues)
        ax.plot(x,y,'g-x',label='mitral '+str(netml.mitB))
        legend(loc='upper left')
        xlabel('distance (microns)')
        ylabel('# of '+netml.centlat_str+'-mit\'s '+netml.weight_str+' conns on joint+multi granules')
        title('# of '+netml.centlat_str+'-mit\'s '+netml.weight_str+' conns on joints/multis bet. mits 0/1 & mit x')
    
        print "READ CAREFULLY: Using",netml.centlat_str,"mitrals'",netml.weight_str,"weights and numbers."
        show()

