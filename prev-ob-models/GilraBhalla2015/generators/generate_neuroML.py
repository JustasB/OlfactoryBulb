#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import math
import pickle

sys.path.extend(["..","../channels","../synapses",\
    "../cells","../networks","../neuroml","../simulations"])

from NeuroML_writer import *
from networkConstants import *
from stimuliConstants import *

from pylab import * # part of matplotlib that depends on numpy but not scipy

from load_channels import *
from load_synapses import *
from load_cells import *
from data_utils import * # has mpi import and variables

import collections # for defaultdict() to count number of occurences of items in a list.

#### USAGE
# (from node000) mpiexec -machinefile ~/hostfile -n <numprocs> ~/Python-2.6.4/bin/python2.6 generate_neuroMl.py [args]
# OR
# python2.6 generate_neuroMl.py [2MITS|2GLOMS] [INVITRO]

## Use the 2MITS option if you want to study connectivity / inhibition between sisters:
## use a small MIT_DISTANCE for 2MITS, else PGs would not be present above mit1;
## set mitralidx=0 and mitralsidekickidx=1 (in networkConstants.py) for 2MITS.

## Use 2GLOMS option for non-sisters which retains mit0 and mit2 of glom0 and glom1:
## But no choice here of distal connection between two gloms --
## mimic in sim by reversing mit0 and mit2;
## set mitralidx=0 and mitralsidekickidx=2 (in networkConstants.py) for 2GLOMS.

##### The two values below set randomization for synapses and mitral somas.
##### These are actually simulation parameters.
##### but I want all randomization only at generator level, not during simulation.
if 'MAXSYNS' in sys.argv: ## deprecated
    # While testing EPSPs, I just use the largest synaptic weights MAX_SYNS=True,
    # to set the max EPSP height below the threshold.
    # Of course KA is still varied to get async granule firing - 
    # just syn weights variation also gives async granule firing.
    # However varied KA causes no discernible difference to the cells soma EPSP!?.
    MAX_SYNS = True
else:
    MAX_SYNS = False
## Na gmax in soma is randomized ## could not be implemented easily using NeuroML, so ditched!!!
RANDOMIZE_MITRALS = False #True ## deprecated
# synapses can be sprinkled isotropically in the mitral dendritic area,
# and connected to the nearest mitral lateral segment.
# OR have 10^4 synapses equiprobable on the lateral dendrites
# and connect each mitral segment to its closest granules.
if 'ISOTROPIC' in sys.argv: ## deprecated
    ISOTROPIC = True
else:
    ISOTROPIC = False
if '2GLOMS' in sys.argv: ## two non-sister mitrals - see USAGE above
    TWOGLOMS = True
    TWOMITS = False
    NUM_GLOMS = 2 # num of modeled gloms
    central_glom = 0
    #DIRECTED_CONNS = [None,None,0,None,None,None]
    DIRECTED_CONNS = [None,None,None,0,0,None]
    #DIRECTED_CONNS = [None,None,None,None] ## Not directed for activity-dep inh
    ## club GRANS_CLUB_SINGLES number of single granules together
    GRANS_CLUB_SINGLES = GRANS_CLUB_SINGLES_2GLOMS
    GRANS_CLUB_JOINTS = GRANS_CLUB_JOINTS_2GLOMS
    mitdistancestr = '_mitdist'+str(mit_distance)
    print "Distance between mit0 and mit2 is",mit_distance,"microns."
    MIT_DISTANCE = mit_distance*1e-6 # meters # distance between the two non-sister mitrals
elif '2MITS' in sys.argv: ## 2 sister mitrals ## deprecated - use 2GLOMS above
    TWOMITS = True
    TWOGLOMS = False
    NUM_GLOMS = 1 # num of modeled gloms - see USAGE above
    central_glom = 0
    DIRECTED_CONNS = [None,0]
    #DIRECTED_CONNS = [None,None] ## Not directed for activity-dep inh
    GRANS_CLUB_SINGLES = 10 # club GRANS_CLUB_SINGLES number of single granules together
    GRANS_CLUB_JOINTS = 1 # we need enough joints to average over - do not club them
    ## don't increase beyond 50 (max 100) microns, as the PGs will not be above mit1.
    mitdistancestr = '_mitdist'+str(mit_distance)
    print "Distance between mit0 and mit1 is",mit_distance,"microns."
    MIT_DISTANCE = mit_distance*1e-6 # meters # distance between the two sister mitrals
else:
    TWOMITS = False
    TWOGLOMS = False
    mitdistancestr = ''
if 'INVITRO' in sys.argv:
    ## half the granules are lost in slicing :
    ## look at my old pre-CNS simulations with various r-dependent connectivity.
    ## even for columns, I assume half of the granules on the prim dend will die.
    ## Instead of half the number of syns, we prune grans more than healthy_slice_width away in y.
    #mit_syns = mit_syns/2 # half of mitral synapses are lost in slice
    INVITRO = True
else: INVITRO = False

synconnfilenamebase = '../netfiles/syn_conn_array_'+str(mit_syns)+\
    '_singlesclubbed'+str(GRANS_CLUB_SINGLES)+\
    '_jointsclubbed'+str(GRANS_CLUB_JOINTS)+\
    '_numgloms'+str(NUM_GLOMS)+'_seed'+str(stim_net_seed)+mitdistancestr
if DIRECTED:
    synconnfilenamebase += '_directed'+str(FRAC_DIRECTED)
    if PROXIMAL_CONNECTION:
        synconnfilenamebase += '_proximal'
    else:
        synconnfilenamebase += '_distal'
if mpisize>1:
    manyprocs = True
    synconnfilenamebase += '_proc'+str(mpirank)
else: manyprocs = False

for arg in sys.argv[1:]:
    synconnfilenamebase += '_'+arg
print "Generating", synconnfilenamebase, "with MAX_SYNS =",MAX_SYNS,\
    ", ISOTROPIC =", ISOTROPIC,\
    ", NUM_GLOMS =",NUM_GLOMS,", DIRECTED =",DIRECTED

class ConnGenerator():

    def __init__(self):
        self.AMPA_factor = AMPA_factor
        self.NMDA_factor = NMDA_factor
        self.GABA_factor = GABA_factor

        self.num_grans = 0
        self.make_gran_conn_arrays()

        # use "micrometer" just for length, but use SI elsewhere.
        self.xmlwriter = NeuroML_writer(length_units="micrometer",units="SI Units")
        self.xmlwriter.note(self.xmlwriter.neuroml,"Fully instantiated network model of Olfactory Bulb"+\
            " by Aditya Gilra, NCBS, Bangalore, India, 2010.")
        load_channels()
        # synchan_activation_correction depends on SIMDT hence typically passed all the way down
        # to granule_mitral_GABA synapse from the top level
        # But here these synapses are loaded just to load the cells,
        # which are loaded to get potential synaptic locations.
        load_synapses(synchan_activation_correction=1.0)
        self.cellSegmentDict = load_cells()
        ### parse the cellSegmentDict returned above to fill in segment lists
        ### that specify potential granule_mitral and PG_mitral synapses
        self.read_mitral_syn_segments()
        print "Total dendritic length for mit|-->gran syns is",self.lateral_length*1e6,"microns."
        self.read_PG_syn_segments()

        self.mitral_z = MITRAL_Z # all mitral somas are in a layer at this z
        self.granule_z = GRANULE_Z # all granule somas are in a layer at this z
        # all PG somas are in a layer at this z
        # I placed them slightly higher than mitrals for easy visualization
        self.PG_z = PG_Z 
        self.PGDict = {}
        self.make_and_xmlwrite_mitrals_and_virtualconnect_to_granules_PGs()

        self.granuleSinglesTable = [] # 1D array of loaded single granules, x,y,z in meters
        self.granuleJointsTable = [] # 1D array of loaded joint granules, x,y,z in meters
        self.granuleMultisTable = [] # 1D array of loaded multi granules, x,y,z in meters
        
        print "Total number of initial granules = ",self.num_grans
        self.club_mitrals_prune_granules()
        print "Clubbed mitrals."
        self.club_granules()
        print "Clubbed granules."
        
        self.set_xml_granules()
        self.set_xml_PGs()
        self.set_xml_granule_connections()
        self.set_xml_PG_connections()
        self.set_xml_ORN_connections()
        self.set_xml_SA_connections()
        self.xmlwriter.writeToFile(synconnfilenamebase+".xml")
        print "Wrote file ", synconnfilenamebase+".xml"

    def make_gran_conn_arrays(self):
        """
        makes a virtual 2D array (virtual=not all will be modeled) of granules
        makes arrays which will hold the connections between the mitral and granules
        see comments next to each array in the function definition
        direct connection to a virtual granule is from to be modeled mitrals i.e. first two in each glom
        indirect connection to a virtual granule is from NOT to be modeled mitrals
        i.e. sister mitrals of the first two in each glom
        """
        ## we require enough granules for gloms up to mit_dend_r radial spread about the central glom,
        ## and then mit_dend_r for the glom at the edge which can be upto glom_r away
        ## and then enough for granule dendrites at the boundaries.

        ## connection array is a 3D array [granx, grany, glomnum] which contains 
        ## whether there is a synapse between indexed granule and indexed mitral:
        ## stores mitral indices as a list
        ## below: perturbation of gran_dend_r is added to synaptic distribution:
        ## have enough slots in array
        granindexsize = int((4*mit_dend_r+2*glom_r+2*gran_dend_r)/gran_spacing)+1
        ## numpy 3D array initialized to None to hold connection lists
        self.connArray = empty( (granindexsize, granindexsize, NUM_GLOMS), dtype=object )
        ## The entry is one if a granule has been created, else None.
        self.granArray = empty( (granindexsize, granindexsize), dtype=object )
        ## numpy 3D array initialized to None-s to hold direct mitral connections
        self.directConnArray = empty( (granindexsize, granindexsize, NUM_GLOMS*MIT_SISTERS), dtype=object )
        ## numpy 3D array initialized to 0-s to hold indirect mitral connections
        self.indirectConnArray = zeros( (granindexsize, granindexsize, NUM_GLOMS), dtype=float32 )
        self.gran_xmaxindex = granindexsize
        self.gran_ymaxindex = granindexsize
        self.gran_xmidindex = granindexsize/2 # integer division
        self.gran_ymidindex = granindexsize/2 # integer division
        print "Made the conn array 3D size: ", granindexsize, granindexsize, NUM_GLOMS

    def make_and_xmlwrite_mitrals_and_virtualconnect_to_granules_PGs(self):
        """
        Makes a virtual number (not all are modeled) of mitrals per glom
        and connects all of them to the virtual granules (not all modeled)
        XML writes only two mitrals per glom which will be
        finally modeled into the the neuroml model.
        """
        self.xmlwriter.set_population(name="mitrals",cell_type="mitral")
        self.glom_positions = [(0.0,0.0)]*NUM_GLOMS
        self.mitral_positions = [(0.0,0.0,0.0)]*(2*NUM_GLOMS)
        ## set the central glom's mitral cells 0 and 1
        ## Need to set these first because 
        ## rest of the gloms must overlap these mitrals, if distal connections.
        for mitnum in [0,1]:
            if not TWOMITS and not TWOGLOMS:
                mitralx = uniform(-glom_r,glom_r)
                mitraly = uniform(-glom_r,glom_r)
                mitralzrotation = uniform(0.0,2*math.pi)
            else:
                ## mit 0's position is at origin if TWOMITS or TWOGLOMS
                if mitnum == 0:
                    mitralx = 0.0
                    mitraly = 0.0
                    mitralzrotation = uniform(0.0,2*math.pi)
                ## mit 1's position separated by MIT_DISTANCE if TWOMITS
                elif TWOMITS:
                    mitralx = 0.0 + MIT_DISTANCE
                    mitraly = 0.0
                    # if DIRECTED in the PROXIMAL_CONNECTION manner,
                    # rotate the central mit1 so that its dendrite overlaps with mit0 soma,
                    # else random rotation.
                    while True:
                        mitralzrotation = uniform(0.0,2*math.pi)
                        if DIRECTED and PROXIMAL_CONNECTION:
                            # important to not have 'else: break' clause for below!
                            if self.dend_thru_central_mit(DIRECTED_CONNS[mitnum],\
                                mitralx, mitraly, mitralzrotation):
                                break
                        else: break
                ## mit1's position is random if TWOGLOMS
                else:
                    mitralx = uniform(-glom_r,glom_r)
                    mitraly = uniform(-glom_r,glom_r)
                    mitralzrotation = uniform(0.0,2*math.pi)
            self.mitral_positions[central_glom*2+mitnum] = (mitralx,mitraly,mitralzrotation)

        ## generate positions of centers of other gloms, mitral positions in another loop after this.
        for glomnum in range(NUM_GLOMS):
            ## central_glom by default (above initialization) has position (0.0,0.0)
            ## set positions for other gloms
            if glomnum != central_glom:
                if TWOGLOMS:
                    if DIRECTED and not PROXIMAL_CONNECTION:
                        ## for distal connections and directed, place the other glom
                        ## on central mit dend, ensuring distance = MIT_DISTANCE.
                        while True:
                            glomx = uniform(-mit_dend_r/2.0, mit_dend_r/2.0)
                            ## set glomy to get lat glom MIT_DISTANCE from central glom
                            if abs(glomx)>MIT_DISTANCE: continue
                            else:
                                glomy = sqrt(MIT_DISTANCE**2 - glomx**2)
                                if self.glom_on_centralmit_dendrites(glomnum, glomx, glomy):
                                    break
                    else:
                        glomx = 0.0 + MIT_DISTANCE
                        glomy = 0.0
                else:
                    ## ensure that other gloms do not step on the central one.
                    ## if distal connections, ensure that other gloms
                    ## fall on central_glom's mit0 or mit1 's lateral dendrites
                    while True:
                        ## repeat until this glom is far enough away from central glom
                        ## force positions of first two gloms to be near central glom
                        ## essentially to check differential inhibition odor responses.
                        if NEAR_2GLOMS and glomnum in [1,2]:
                            glomx = uniform(-glom_r*2.0, glom_r*2.0)
                            glomy = uniform(-glom_r*2.0, glom_r*2.0)
                        else:
                            glomx = uniform(-mit_dend_r/2.0, mit_dend_r/2.0)
                            glomy = uniform(-mit_dend_r/2.0, mit_dend_r/2.0)
                        ## repeat until this glom is far enough away from central glom
                        if abs(glomx) > glom_r and abs(glomy) > glom_r:
                            if DIRECTED and not PROXIMAL_CONNECTION:
                                ## Distal directed connections
                                ##    ensure that other gloms lie on central glom's lateral dendrite
                                ## important to not have 'else: break' clause for below!
                                if self.glom_on_centralmit_dendrites(glomnum, glomx, glomy):
                                    break
                            else: break
                self.glom_positions[glomnum] = (glomx,glomy)

        ### xmlwrite the mitral cells and make the PG and granules connections.
        ### Each glom has its constellation of clubbed PG cells. For PGs I club first and then connect
        ### as I do not assume differential excitation of PG cells within a glomerulus.
        self.PGConnArray = empty( (NUM_GLOMS, PGindexsize, PGindexsize), dtype=object )
        for glomnum in range(NUM_GLOMS):
            glomx,glomy = self.glom_positions[glomnum]
            for mitralnum in range(mits_per_glom):
                print "generating for glomerulus ", glomnum, " mitral ", mitralnum
                # Do not generate mitral positions for the central glom's mit0 and mit1.
                if glomnum == central_glom and mitralnum in [0,1]:
                    mitralx,mitraly,mitralzrotation = self.mitral_positions[2*glomnum + mitralnum]
                else:
                    if TWOGLOMS and PROXIMAL_CONNECTION and mitralnum==0:
                        ## if proximal, set mit2 exactly MIT_DISTANCE away
                        ## actually, even for this, I should set mit3,
                        ## else jitter in tuftinh with distance.
                        mitralx = glomx
                        mitraly = glomy
                    elif TWOGLOMS and not PROXIMAL_CONNECTION and mitralnum==1:
                        ## if distal, set mit3 exactly MIT_DISTANCE away
                        mitralx = glomx
                        mitraly = glomy
                    else:
                        mitralx = glomx + uniform(-glom_r,glom_r)
                        mitraly = glomy + uniform(-glom_r,glom_r)
                    ## generate mitral rotation about z-axis
                    dend_on_centralmit_iter = 0
                    good_factor = 3.0   ## first aim for at least 3*gran_dend_r
                                        ## within rectangle of side 2*gran_dend_r
                    while True:
                        mitralzrotation = uniform(0.0,2*math.pi)
                        ## if DIRECTED, for every mit0, mit1 of 
                        ## glomerulus other than central glom,
                        ## ensure it intersects with central mit's soma.
                        if DIRECTED and PROXIMAL_CONNECTION and mitralnum in [0,1]:
                            centralmitnum = DIRECTED_CONNS[glomnum*2+mitralnum]
                            if centralmitnum is None: break
                            ## important to not have 'else: break' clause for below!
                            if self.dend_thru_central_mit(centralmitnum,\
                                    mitralx, mitraly, mitralzrotation, good_factor):
                                break
                        else: break
                        dend_on_centralmit_iter += 1
                        ## reduce the fitting factor by 0.1 if not finding a good enough overlap.
                        if (dend_on_centralmit_iter%1000)==0:
                            good_factor -= 0.1
                            if good_factor<2.5: print "Not such good overlap, factor =",good_factor
                if mitralnum in [0,1]:
                    self.mitral_positions[glomnum*2+mitralnum] = (mitralx,mitraly,mitralzrotation)
                    ## xml file has lengths in microns, hence 1e6 factor,
                    ## mitrals' somas are 300 microns above granules' somas
                    self.xmlwriter.set_instance(elementid=str(glomnum*2+mitralnum),\
                        x=str(mitralx*1e6),y=str(mitraly*1e6),z=str(self.mitral_z*1e6),\
                        zrotation=mitralzrotation)
                    ### Presently apart from translation, there is no way
                    ### to tweak some parts of a cell in neuroml
                    #if RANDOMIZE_MITRALS:
                    #    for NaChanID in get_matching_children(mitralsoma, ['Na']):
                    #        NaChan = moose.HHChannel(NaChanID)
                    #        NaChan.Gbar = NaChan.Gbar*uniform(1.0-SYN_VARIATION,1.0+SYN_VARIATION)
                    #    #mitral._mitralDendNa.Gbar = mitral._mitralDendNa.Gbar*\
                    #    #      uniform(1.0-SYN_VARIATION,1.0+SYN_VARIATION)
                    ##### Connect randomly uniformly to PGs in the glom.
                    self.PGConnect(glomx, glomy, glomnum, mitralnum)
                self.granuleConnect(mitralx, mitraly, glomnum, mitralnum, mitralzrotation, 1.0)

        ## after generating the positions of mitrals,
        ## and connecting mitrals randomly,
        ## connect some already created granules (will be mostly singles)
        ## directed-ly between centralmit and other mitrals (either proximally or distally).
        #### NOTE: granuleConnect allows multiple same-mit |--> same-gran connections
        #### But granuleConnectDirected does not
        ## connect the directed joints for modelled mitrals only
        for mitnum in range(NUM_GLOMS*MIT_SISTERS):
            if DIRECTED and (DIRECTED_CONNS[mitnum] is not None):
                centralmitnum = DIRECTED_CONNS[mitnum]
                print "connecting directed joints between mitral", mitnum,\
                    "and central mitral",centralmitnum
                self.granuleConnectDirected(centralmitnum,mitnum,FRAC_DIRECTED)

    def good_dendritic_overlap(self, sqdistances, sqdistcutoff,\
        seglengths, min_length, num_segs = None):
        """
        The total length of all the segments within sqdistcutoff in sqdistances
        should be > min_length (if num_segs is None) to qualify as good overlap.
        If num_segs is not None, num_segs number of segments must lie within sqdistcutoff.
        """
        goodsegidxs = where(array(sqdistances)<sqdistcutoff)[0]
        if len( goodsegidxs ) > 0:
            ## if num_segs is None, check if total length is sufficient
            if num_segs is None:
                goodseglengths = array(seglengths)[goodsegidxs]
                if sum(goodseglengths)>min_length:
                    return True
                else: return False
            ## else just check if number of segments is sufficient
            elif len( goodsegidxs ) >= num_segs:
                return True
            else: return False
        else: return False

    def glom_on_centralmit_dendrites(self, glomnum, glomx, glomy):
        """
        centralmitnum is DIRECTED_CONNS[mitnum], where mitnum is mit0 of glomnum.
        if (glomx,glomy) lie on central mitral's lateral dendrite, return True else False
        at least two segments of centralmitnum's
        lateral dendrites must be within glom_r of the (glomx,glomy)
        See also odor_responses.py:
         Lower half of gloms receive only odor A, upper half receive only odor B.
        """
        centralmitnum = DIRECTED_CONNS[2*glomnum]
        ## for 2GLOMS, mit3 (not mit2) is dir onto mit0, see DIRECTED_CONNS at top of file
        if centralmitnum is None: centralmitnum = DIRECTED_CONNS[2*glomnum+1]
        mitx, mity, mitzrot = self.mitral_positions[centralmitnum]
        distances = []
        seglengths = []
        for (segid,(segx1,segy1,segz1),(segx2,segy2,segz2),seglength,segdia) \
                in self.granule_mitral_syn_positions:
            segx,segy = self.rotate((segx1+segx2)/2.0,(segy1+segy2)/2.0,mitzrot)
            ## do not include z in the distance calculation since:
            ## 1) it was not used while assigning the granule connections,
            ## 2) if z is used, it creates inward pointing granule-mitral connections (towards mitral center).
            dist = sqrt((glomx-segx-mitx)**2+(glomy-segy-mity)**2)
            distances.append( dist )
            seglengths.append( seglength )
        ## the total lengths of segments within glom_dia must be > glom dia.
        return self.good_dendritic_overlap(distances, 2*glom_r, seglengths, 2*glom_r)

    def dend_thru_central_mit(self, centralmitnum, mitx, mity, mitzrot, goodfactor=2.6):
        """
        if enough dendritic length of a mitral at (mitx,mity) with rotation mitzrot lies on 
        centralmitnum's soma, return True else False
        See also odor_responses.py: Lower half of gloms receive only odor A, upper half receive only odor B.
        """
        centralmitx, centralmity, centralmitzrot = self.mitral_positions[centralmitnum]
        xmin = centralmitx-gran_dend_r
        ymin = centralmity-gran_dend_r
        xmax = centralmitx+gran_dend_r
        ymax = centralmity+gran_dend_r
        distances = []
        seglengths = []
        seglengthwithin = 0.0
        for (segid,(segx1,segy1,segz1),(segx2,segy2,segz2),seglength,segdia) \
                in self.granule_mitral_syn_positions:
            segx1,segy1 = self.rotate(segx1,segy1,mitzrot)
            segx2,segy2 = self.rotate(segx2,segy2,mitzrot)
            ### do not include z in the distance calculation since:
            ### 1) it was not used while assigning the granule connections,
            ### 2) if z is used, it creates inward pointing granule-mitral connections (towards mitral center).
            ### calculate the minimum distance (data_utils.py) from centralmit's soma to this dendritic segment
            #dist = minimum_distance((mitx+segx1,mity+segy1),(mitx+segx2,mity+segy2),(centralmitx,centralmity))
            #distances.append( dist )
            #seglengths.append( seglength )
            
            ## only take the length of segments within glom_dend_r of centralmit soma (rectangular region).
            accept,x1,y1,x2,y2 = clip_line_to_rectangle(\
                mitx+segx1,mity+segy1,mitx+segx2,mity+segy2,xmin,ymin,xmax,ymax)
            if accept: seglengthwithin += norm( array((x1,y1))-array((x2,y2)) )
        ### a mitral segment may only connect to a granule whose soma is gran_dend_r apart.
        ### thus, two mitral segments will barely connect via granules if they are 2*gran_dend_r apart.
        ### 500 micron apart mitrals are not able to satify gran_dend_r/2.0 distance bound
        ### the total lengths of segments within 1.5*gran_dend_r must be > 3*gran_dend_r
        ##return self.good_dendritic_overlap(distances, (1.5*gran_dend_r)**2, seglengths, 3*glom_r)
        ### just one segment within gran_dend_r/2.0 is required
        #return self.good_dendritic_overlap(distances, (gran_dend_r/2.0)**2, seglengths, 0.0, num_segs=1)
        
        ## The dendritic length across the rectangular patch of side 2*gran_dend_r
        ## must be at least 2.6*gran_dend_r = good_factor*gran_dend_r.
        ## A fork in the segments or a near diagonal traversal sqrt(2)=1.4 will achieve this.
        ## sqrt(2)*2*glom_dend_r = 2.828*glom_dend_r ~= goodfactor*gran_dend_r
        if seglengthwithin > goodfactor*gran_dend_r: return True
        else: return False

    def club_mitrals_prune_granules(self):
        """
        This function decimates away the extra mitrals in the glomerulus, keeping only two: 0 and 1.
        It populates two arrays: directConnArray and indirectConnArray.
        directConnArray[x][y][glomnum*2+0|1] is set to mit_seg_id if granule at x,y
            is connected to mitral 0|1 of glomerulus glomnum.
        indirectConnArray[x][y][glomnum*2+0|1] is set to the number of mitrals of glomerulus glomnum
            that connect to granule at x,y if it is already connected to mitral 0|1.
        granArray is pruned to remove non-useful granules that do not connect to mitrals 0|1 of any glomerulus.
        If in vitro, granArray is also pruned of those > healthy_slice_width away in y.
        #commented# granArray is also pruned of joints that are not between pairs given in DIRECTED_CONNS (networkConstants.py).
        """
        healthy_slice_y = int(healthy_slice_width/gran_spacing)
        for x in range(self.gran_xmaxindex):
            for y in range(self.gran_ymaxindex):
                ## If no granule here, go to next point in grid.
                if self.granArray[x][y] is None: 
                    continue

                ## if in vitro, prune granules that are more than healthy_slice_width away in y
                ## mits are separated parallel to slice cut in x
                if INVITRO:
                    if abs(y-self.gran_ymidindex)>healthy_slice_y:
                        self.granArray[x][y]=None
                        continue

                ### If columnar, prune away granules outside a glomerular column (singles, joints and multis)
                #if COLUMNAR:
                #    xcont=(x-self.gran_xmidindex)*gran_spacing
                #    ycont=(y-self.gran_ymidindex)*gran_spacing
                #    # make a list of True/False whether (x,y) falls inside a modeled glomerulus.
                #    within_column_list = \
                #        [ abs(self.glom_positions[glomnum][0]-xcont)<=2*glom_r \
                #        and abs(self.glom_positions[glomnum][1]-ycont)<=2*glom_r \
                #        for glomnum in range(NUM_GLOMS) ]
                #    # if granule at (x,y) is not inside a modelled glom,
                #    # discard it and carry on to next granule
                #    if True not in within_column_list:
                #        self.granArray[x][y] == None
                #        continue

                useful = False
                gran_conns = {}
                for glomnum,glom_mits_segs in enumerate(self.connArray[x][y]):
                    if glom_mits_segs is not None:
                        glom_mits, mit_segs = zip(*glom_mits_segs)
                        glom_mits = list(glom_mits)
                        mit_segs = list(mit_segs)
                        granxy_to_mit0 = False
                        granxy_to_mit1 = False
                        while True:
                            ## create reciprocal synapses for all the connections between mit0 and this gran
                            if (0 in glom_mits):
                                useful = True
                                ## fill in the seg for the mit0
                                ## note that glom_mits and mit_segs are parallel lists: see zip(* ) above
                                ## store the segid of mitral's connection to this gran
                                mit0_idx = glom_mits.index(0)
                                if not granxy_to_mit0: # first connection
                                    self.directConnArray[x][y][glomnum*2] = [ mit_segs[mit0_idx] ]
                                else: # later connections are appended
                                    self.directConnArray[x][y][glomnum*2].append(mit_segs[mit0_idx])
                                granxy_to_mit0 = True
                                ## remove the the mitnum and segid from the parallel lists for next round
                                del glom_mits[mit0_idx]
                                del mit_segs[mit0_idx]
                            else: break
                        while True:
                            ## create reciprocal synapses for all the connections between mit1 and this gran
                            if (1 in glom_mits):
                                useful = True
                                ## fill in the seg for the mit0
                                ## note that glom_mits and mit_segs are parallel lists: see zip(* ) above
                                ## store the segid of mitral's connection to this gran
                                mit1_idx = glom_mits.index(1)
                                if not granxy_to_mit1: # first connection
                                    self.directConnArray[x][y][glomnum*2+1] = [ mit_segs[mit1_idx] ]
                                else: # later connections are appended
                                    self.directConnArray[x][y][glomnum*2+1].append(mit_segs[mit1_idx])
                                granxy_to_mit1 = True
                                ## remove the the mitnum and segid from the parallel lists for next round
                                del glom_mits[mit1_idx]
                                del mit_segs[mit1_idx]
                            else: break
                        ## for the connections from unmodelled mitrals of this glom,
                        ## fill up the indirectConnArray.
                        for mitnum in glom_mits:
                            ## extra excitation from unmodelled mits of this glom is accounted for:
                            ## but reduced by factor mitspread_extraexc_redux to compensate
                            ## that I have not spatially spread the mitrals enough
                            ## (only 100 microns vs 200-300 microns spread - [Sosulski et al 2011])
                            ## while generating the mitral positions above.
                            self.indirectConnArray[x][y][glomnum] += mitspread_extraexc_redux
                                
                if not useful:
                    self.granArray[x][y] = None # if not connected to a useful mitral, just remove the granule

                ### prune joints that do not join required mitrals.
                #if DIRECTED and NUM_GLOMS>1:
                #    ## where() doesn't work with '!=None' or 'is not None' 
                #    #mitindices = where(self.directConnArray[x][y]!=None)[0]
                #    mitindices = []
                #    for i,segid in enumerate(self.directConnArray[x][y]):
                #        if segid is not None: mitindices.append(i)
                #    numconns = len(mitindices)
                #    if numconns == 2 : # jointly connected granules
                #        mit1index = mitindices[0] # index of the first mitral
                #        mit2index = mitindices[1] # index of the second mitral
                #        if DIRECTED_CONNS[mit1index]!=mit2index and DIRECTED_CONNS[mit2index]!=mit1index:
                #            # if not connected between required mitrals, just remove the granule
                #            self.granArray[x][y] = None

        ## after clubbing mitrals and sorting & pruning granules,
        ## there is no need to keep the detailed connArray.
        del self.connArray

    def get_seg_pos(self,mitnum,segid):
        """Returns the x,y of the midpoint of the segment segid of mitnum."""
        mitx,mity,mitzrot = self.mitral_positions[mitnum]
        mitz = self.mitral_z
        for seg in self.granule_mitral_syn_positions:
            ## seg = (segid,(segx1,segy1,segz1),(segx2,segy2,segz2),seglength,segdia)
            if segid == seg[0]:
                segx = (seg[1][0]+seg[2][0])/2.0
                segy = (seg[1][1]+seg[2][1])/2.0
                segx,segy = self.rotate(segx,segy,mitzrot)
        return mitx+segx,mity+segy

    def club_granules(self):
        """
        This function parses the directConnArray and indirectConnArray into granules 
            that only connect to a single useful mitral, or to two mitrals, or to more than two mitrals.
        granuleSinglesTable is a list of (x,y, mitindex, (mit_segment_id,num_conns), gran_segment_id, extra_excitation):
            each of which represents a clubbed granule
            i.e. GRANS_CLUB_SINGLES number of single granules connected to mitindex.
            x,y is the position of the last granule out of those clubbed.
            mit_segment_id is the first segment in the list of segments
            at which mitindex is connected to the above last granule.
            num_conns are the number of connections between mitindex and ALL clubbed granules.
            extra_excitation is the average coherent excitation
            that a clubbed granule receives from sister mitrals i.e. belonging to the same glomerulus.
        granuleJointsTable clubs granule connected to two mitrals. It is a list of:
            (x,y, mitindex1, (mit1_segment_id,num_conns1), gran_segment_id1, extraexc_from_mit1sisters,\
             mitindex2, (mit2_segment_id,num_conns2), gran_segment_id2, extraexc_from_mit2sisters)
        granuleMultisTable does not club granules connected to multiple i.e. > 2 mitrals
            as the combinatorics is too much to be useful for clubbing.
            It is a list of single granules each of the form:
            (x,y,[(mitindex,(mit_segment_id,num_conns),gran_segment_id,extra_exc_from_mitindex_sisters),...])
        """
        nummits = NUM_GLOMS*MIT_SISTERS
        gransnumarray = [0]*nummits # singles number
        joints = 0 # joints number
        multis = 0 # multiply connected number

        clubnumarray = [0]*nummits
        excnumarray = [0]*nummits
        clubnumarray_twos = zeros((nummits,nummits,2), dtype=float32)
        excnumarray_twos = zeros((nummits,nummits,2), dtype=float32)
        zcont=self.granule_z
        for x in range(self.gran_xmaxindex):
            for y in range(self.gran_ymaxindex):
                ## If no granule here, go to next point in grid.
                if self.granArray[x][y] is None:
                    continue
                ## where() doesn't work with '!=None' or 'is not None' 
                #mitindices = where(self.directConnArray[x][y]!=None)[0]
                mitindices = []
                for i,segids in enumerate(self.directConnArray[x][y]):
                    if segids is not None: mitindices.append((i,segids))
                numconns = len(mitindices)
                ## singly connected granules i.e. reciprocal syns to only one mit.
                ## extra excitation could be from many
                if numconns == 1 :
                    ## index of the single mitral connected
                    mitindex = mitindices[0][0]
                    ## len(segids) is the number of reciprocal syns
                    ## from this mitindex to this gran
                    clubnumarray[mitindex] += len(mitindices[0][1])
                    ## This granule is connected to mitindex mitral,
                    ## it may have coherent input from sister mitrals (same glomerulus)
                    ## Thus put in extra excitation from the same glomerulus,
                    ## the rest i.e. non-sister excitation goes into the background
                    ## extra one way excitatory synapses on the granule
                    excnumarray[mitindex] += self.indirectConnArray[x][y][mitindex/2]
                    ## make 1 granule for every GRANS_CLUB_SINGLES number of granules
                    ## could increment clubnumarray[mitindex] by larger than 1 per iteration, hence >= below
                    if clubnumarray[mitindex] >= GRANS_CLUB_SINGLES:
                        xcont=(x-self.gran_xmidindex)*gran_spacing
                        ycont=(y-self.gran_ymidindex)*gran_spacing
                        ## take only the first mit seg connected to this gran
                        mit_seg_id = self.directConnArray[x][y][mitindex][0]
                        num_conns = clubnumarray[mitindex]
                        extra_weight = float(excnumarray[mitindex])/num_conns
                        self.granuleSinglesTable.append(\
                            (xcont, ycont, mitindex, (mit_seg_id,num_conns),'1', extra_weight) )
                        gransnumarray[mitindex] += 1
                        if extra_weight>1.0:
                            print "Finished setting up single granule number",gransnumarray[mitindex],\
                                "to mitral",mitindex,"extra excitation weight =",extra_weight
                        clubnumarray[mitindex] = 0
                        excnumarray[mitindex] = 0
                ## Create the joint granules that are connected to two mitrals
                elif numconns == 2:
                    mit1index = mitindices[0][0] # index of the first mitral
                    mit2index = mitindices[1][0] # index of the second mitral
                    clubnumarray_twos[mit1index][mit2index][0] += len(mitindices[0][1]) # no of segids
                    clubnumarray_twos[mit1index][mit2index][1] += len(mitindices[1][1]) # no of segids
                    #### always mit1index<mit2index,
                    #### hence the above matrix will have zeroes on diagonal and one side.
                    ## This granule is near mitral1 and mitral2,
                    ## put in the extra excitation from other sister mitrals
                    ## extra one way excitatory synapses on the granule from mit1index's glom
                    excnumarray_twos[mit1index][mit2index][0] += self.indirectConnArray[x][y][mit1index/2]
                    ## extra one way excitatory synapses on the granule from mit2index's glom
                    excnumarray_twos[mit1index][mit2index][1] += self.indirectConnArray[x][y][mit2index/2]
                    ## Aggregate/club GRANS_CLUB_JOINTS number of joint granules together
                    ## when number of connections from mit1 or mit2 reaches GRANS_CLUB_JOINTS
                    conn_mit1 = clubnumarray_twos[mit1index][mit2index][0]
                    conn_mit2 = clubnumarray_twos[mit1index][mit2index][1]
                    if max(conn_mit1,conn_mit2) >= GRANS_CLUB_JOINTS:
                        xcont=(x-self.gran_xmidindex)*gran_spacing
                        ycont=(y-self.gran_ymidindex)*gran_spacing
                        ## take only the first mit1 and mit2 seg-s connected to this gran
                        mit1_seg_id = self.directConnArray[x][y][mit1index][0]
                        mit2_seg_id = self.directConnArray[x][y][mit2index][0]
                        extra_weight1 = \
                            float(excnumarray_twos[mit1index][mit2index][0])/conn_mit1
                        extra_weight2 = \
                            float(excnumarray_twos[mit1index][mit2index][1])/conn_mit2
                        ## if mit1 and mit2 are from same glom, 
                        ## extra exc is being counted twice, so / 2.0
                        if mit1index/2 == mit2index/2: # integer div
                            extra_weight1 /= 2.0
                            extra_weight2 /= 2.0
                        self.granuleJointsTable.append((xcont,ycont, mit1index, (mit1_seg_id,conn_mit1), '1',\
                            extra_weight1, mit2index, (mit2_seg_id,conn_mit2), '1', extra_weight2))
                        joints += 1
                        if extra_weight1>1.0 or extra_weight2>1.0:
                            print "Finished setting up joint granule number",joints,"between",\
                                mit1index,"and",mit2index,"with extraweights =",\
                                extra_weight1,"and",extra_weight2
                        clubnumarray_twos[mit1index][mit2index][0] = 0
                        clubnumarray_twos[mit1index][mit2index][1] = 0
                        excnumarray_twos[mit1index][mit2index][0] = 0
                        excnumarray_twos[mit1index][mit2index][1] = 0
                ## Create the multiply-connected granules
                ## that are connected to more than two mitrals
                elif numconns > 2:
                    synconns = []
                    xcont=(x-self.gran_xmidindex)*gran_spacing
                    ycont=(y-self.gran_ymidindex)*gran_spacing
                    for mitnum,seg_ids in mitindices:
                        extra_weight = self.indirectConnArray[x][y][mitnum/2]
                        ## Take only the first mit seg connected to this gran
                        synconns.append((mitnum, (seg_ids[0],len(seg_ids)), '1', extra_weight))
                    self.granuleMultisTable.append((xcont,ycont,synconns))
                    multis += 1
                    #print "Finished setting up multi granule number",multis,"between",synconns
        for i in range(nummits):
            print "The number of single granules to mitral", i,"after clubbing",\
                GRANS_CLUB_SINGLES,"together =",gransnumarray[i]
        print "The number of joint granules after clubbing",GRANS_CLUB_JOINTS,"together =",joints
        print "The number of multi granules is",multis

    def set_xml_granules(self):
        z=self.granule_z*1e6
        ## singly-connected granules
        self.xmlwriter.set_population(name="granules_singles",cell_type="granule")
        gran_num = 0
        for gran in self.granuleSinglesTable:
            # xml file has lengths in microns, hence 1e6 factor
            self.xmlwriter.set_instance(elementid=str(gran_num),\
                x=str(gran[0]*1e6),y=str(gran[1]*1e6),z=str(z))
            gran_num += 1
        ## jointly connected granules
        self.xmlwriter.set_population(name="granules_joints",cell_type="granule")
        gran_num = 0
        for gran in self.granuleJointsTable:
            # xml file has lengths in microns, hence 1e6 factor
            self.xmlwriter.set_instance(elementid=str(gran_num),\
                x=str(gran[0]*1e6),y=str(gran[1]*1e6),z=str(z))
            gran_num += 1
        ## multiply-connected granules
        ## For 2 mitrals there will not be any multi-s.
        if len(self.granuleMultisTable) != 0:
            self.xmlwriter.set_population(name="granules_multis",cell_type="granule")
            gran_num = 0
            for gran in self.granuleMultisTable:
                # xml file has lengths in microns, hence 1e6 factor
                self.xmlwriter.set_instance(elementid=str(gran_num),\
                    x=str(gran[0]*1e6),y=str(gran[1]*1e6),z=str(z))
                gran_num += 1

    def set_xml_PGs(self):
        z=self.PG_z*1e6
        self.xmlwriter.set_population(name="PGs",cell_type="PG")
        self.PG_num = 0
        self.PGPopulation = []
        for glomnum,PGfull in enumerate(self.PGConnArray):
            for xidx,PGx in enumerate(PGfull):
                for yidx,PG in enumerate(PGx):
                    if PG is not None:
                        # length in microns in the xml file, hence the 1e6 factor
                        self.xmlwriter.set_instance( elementid=str(self.PG_num),\
                            x=str((self.glom_positions[glomnum][0]-PG_halfextent+xidx*PG_spacing)*1e6),\
                            y=str((self.glom_positions[glomnum][1]-PG_halfextent+yidx*PG_spacing)*1e6),\
                            z=str(z) )
                        self.PG_num += 1
        print "Total number of PG cells (clubbed "+str(PG_CLUB)+") =", self.PG_num

    def get_delay(self,delay_distrib,delay_spread):
        if delay_distrib == 0: return delay_spread ## no distribution, assume actual delay passed in as delay_spread
        elif delay_distrib == 1: return uniform(0.0,delay_spread) ## uniform distribution
        elif delay_distrib == 2: return exponential(scale=delay_spread) ## exponential distribution
        elif delay_distrib == 3: return gamma(shape=2.0,scale=1.0)*delay_spread ## Gamma distribution

    def make_synproplist(self,synapselist,strongsynfactor=1.0):
        if not MAX_SYNS:
            ## if MAX_SYNS is false, distribute synaptic weights lognormally or uniformly
            if lognormal_weights:
                ## Song et al PLOS Biology 2005 find lognormal distribution of synaptic weights
                ## http://en.wikipedia.org/wiki/Log-normal_distribution 
                ## mean and sigma are for the underlying normal distribution of log(wt).
                ## if m and d be the mean and s.d. of wt, then:
                ## sigma = sqrt(log(1+d^2/m^2)) and mean = log(m) - sigma^2/2
                ## m = 1.0 and d = SYN_VARIATION*m
                sigma = sqrt(log(1+SYN_VARIATION**2))
                synproplist = [ (syn[0],\
                    lognormal(mean=log(1.0)-sigma**2/2.0,sigma=sigma)*syn[1]*strongsynfactor, syn[2])\
                    for syn in synapselist ] # synapse_type (name of synapse) and weight and delay
            else:
                ## uniform distribution between 1.0-SYN_VARIATION and 1.0+SYN_VARIATION
                synproplist = [ (syn[0],syn[1]*uniform(1.0-SYN_VARIATION,1.0+SYN_VARIATION)*strongsynfactor,syn[2])\
                    for syn in synapselist ] # synapse_type (name of synapse) and weight and delay
        else:
            # if MAX_SYNS is True, set synaptic weights to 1.0+SYN_VARIATION times normal
            # I still make a call to the random number generator,
            # so that synaptic connections remain the same!
            synproplist = []
            for syn in synapselist:
                if lognormal_weights:
                    sigma = sqrt(log(1+SYN_VARIATION**2))
                    lognormal(log(1.0)-sigma**2/2.0,sigma=sigma)
                else:
                    uniform(1.0-SYN_VARIATION,1.0+SYN_VARIATION)
                synproplist.append( (syn[0],syn[1]*(1.0+SYN_VARIATION)*strongsynfactor,syn[2]) )
                # synapse_type (name of synapse) and weight and delay
        return synproplist

    def set_xml_PG_connections(self):
        m_pg_AMPA_weight = self.AMPA_factor
        #m_pg_NMDA_weight = self.NMDA_factor
        pg_m_GABA_weight = self.GABA_factor
        ##################### excitatory mitral->PG
        self.xmlwriter.set_projection(name='mitral_PG',source='mitrals',target='PGs')
        self.xmlwriter.set_synapse_props(synapse_type='mitral_PG',\
            weight=str(m_pg_AMPA_weight),threshold=THRESHOLD,delay=exc_synaptic_delay)
        PG_num = 0
        conn_num = 0
        for glomnum,PGfull in enumerate(self.PGConnArray):
            for xidx,PGx in enumerate(PGfull):
                for yidx,PG in enumerate(PGx):
                    if PG is not None:
                        if INVITRO:
                            num_pg_stim = 1
                        else:
                            ## distribute the unconnected mit->PG spines
                            ## to get inputs from the connected mit->PG spines
                            ## as a proxy for the unmodelled mitrals
                            num_pg_stim = (NUM_PG_to_M_T_ET_SPINES-len(PG))/len(PG) # integer division
                            if num_pg_stim<1: num_pg_stim=1 # have at least one mit->PG synapse
                        for conn in PG:
                            ## 'num_pg_stim' number of excitatory mitral->PG synapses
                            ## are created from the modelled mitral to clubbed PG
                            ## to also take into account unmodelled mitrals' and ET's input.
                            ## these are randomly delayed so as not to be synchronous.
                            for multisyn in range(num_pg_stim):
                                delay = self.get_delay(exc_delay_distrib,mitral_PG_delay_spread)
                                synapselist = [ ('mitral_PG',m_pg_AMPA_weight,delay) ]
                                synproplist = self.make_synproplist(synapselist)
                                self.xmlwriter.set_connection(elementid=str(conn_num),\
                                    pre_cell_id=str(conn[0]),post_cell_id=str(PG_num),\
                                    pre_segment_id=conn[1],post_segment_id=conn[2], synproplist=synproplist)
                                conn_num += 1
                        PG_num += 1
        ################ inhibitory PG-|mitral
        self.xmlwriter.set_projection(name='PG_mitral',source='PGs',target='mitrals')
        ## set default values, but these are overridden later.
        self.xmlwriter.set_synapse_props(synapse_type='PG_mitral',\
            weight=str(pg_m_GABA_weight),threshold=THRESHOLD,delay=inh_synaptic_delay)
        PG_num = 0
        conn_num = 0
        for glomnum,PGfull in enumerate(self.PGConnArray):
            for xidx,PGx in enumerate(PGfull):
                for yidx,PG in enumerate(PGx):
                    if PG is not None:
                        for conn in PG:
                            ## extra inhibitory PG-|mitral synapses are created
                            ## from the modelled PG to take into account unmodelled PGs' input.
                            ## these are randomly delayed so as not to be synchronous.
                            ###### IMPORTANT: For the granule cells, the first spike has minimal delay
                            ## whereas here, I've set delay for the first spike too.... 
                            ## Along with the spread of the excitatory synapses, this too can be thought
                            ## of as taking care of spread in inputs due to ET cells and various mitral cells.
                            ## Though, ideally, the first spike should not be delayed, but it doesn't change
                            ## the tuft inhibition results at least, probably only more linear!!
                            for multisyn in range(PG_CLUB):
                                ## uniform or exponential delay over PG_mitral_delay_spread
                                delay = self.get_delay(inh_delay_distrib,PG_mitral_delay_spread)
                                synapselist = [ ('PG_mitral',pg_m_GABA_weight,delay) ]
                                synproplist = self.make_synproplist(synapselist)
                                self.xmlwriter.set_connection(elementid=str(PG_num),\
                                    pre_cell_id=str(PG_num),post_cell_id=str(conn[0]),\
                                    pre_segment_id=conn[2],post_segment_id=conn[1],synproplist=synproplist)
                                conn_num += 1
                        PG_num += 1
        ################ inhibitory file-|mitral connections at tuft
        if random_inh:
            self.xmlwriter.set_projection(name='tuft_mitral',source='file',target='mitrals')
            self.xmlwriter.set_synapse_props(synapse_type='PG_mitral',\
                weight=str(pg_m_GABA_weight),threshold=THRESHOLD,delay=inh_synaptic_delay)
            PG_num = 0
            conn_num = 0
            for glomnum,PGfull in enumerate(self.PGConnArray):
                for xidx,PGx in enumerate(PGfull):
                    for yidx,PG in enumerate(PGx):
                        if PG is not None:
                            for conn in PG:
                                for multisyn in range(PG_CLUB):
                                    ## uniform delay over PG_mitral_delay_spread
                                    delay = self.get_delay(inh_delay_distrib,PG_mitral_delay_spread)
                                    synapselist = [ ('PG_mitral',pg_m_GABA_weight,delay) ]
                                    synproplist = self.make_synproplist(synapselist)
                                    self.xmlwriter.set_connection(elementid=str(PG_num),\
                                        pre_cell_id='file',post_cell_id=str(conn[0]),\
                                        pre_segment_id=str(int(uniform(num_gran_baseline_files))),\
                                        post_segment_id=conn[1],synproplist=synproplist)
                                    conn_num += 1
                            PG_num += 1

    def set_xml_granule_connections(self):
        m_g_AMPA_weight = self.AMPA_factor
        m_g_NMDA_weight = self.NMDA_factor
        g_m_GABA_weight = self.GABA_factor
        ##### Separate connection lists are made for single and joint granules
        ################## Excitatory
        for multiplicity in ['singles','joints','multis']:
            ## extra excitatory connections AMPA and NMDA from pruned sister mitral cells
            ## proxied by extra excitation from modeled mitral cells 0 and 1.
            ## By setting CLUB_MITRALS = True/False in simset_*.py, you can switch at runtime
            ## between above behaviour / provide extra Poisson background as proxy.
            self.set_xml_connections_by_multiplicity('mitral_granule_extra_exc_'+\
                multiplicity,2,multiplicity,\
                [("mitral_granule_AMPA",m_g_AMPA_weight,THRESHOLD,exc_synaptic_delay),\
                ("mitral_granule_NMDA",m_g_NMDA_weight,THRESHOLD,exc_synaptic_delay)],\
                exc_delay_distrib, mitral_granule_delay_spread)
            ## main excitation of mitrals on to granules.
            if SATURATING_SINGLE_INHIBITION:
                ## excitatory saturating AMPA and NMDA from mitral to granule
                self.set_xml_connections_by_multiplicity('mitral_granule_main_exc_'+\
                    multiplicity,1,multiplicity,\
                    [("mitral_granule_saturatingAMPA",m_g_AMPA_weight,THRESHOLD,exc_synaptic_delay),\
                    ("mitral_granule_saturatingNMDA",m_g_NMDA_weight,THRESHOLD,exc_synaptic_delay)],\
                    0, exc_synaptic_delay)
            else:
                ## excitatory NON-saturating i.e. usual AMPA and NMDA from mitral to granule
                self.set_xml_connections_by_multiplicity('mitral_granule_main_exc_'+\
                    multiplicity,1,multiplicity,\
                    [("mitral_granule_AMPA",m_g_AMPA_weight,THRESHOLD,exc_synaptic_delay),\
                    ("mitral_granule_NMDA",m_g_NMDA_weight,THRESHOLD,exc_synaptic_delay)],\
                    0, exc_synaptic_delay)

        ################### Inhibitory
        ## singles granule-mitral inhibition is asynchronously delayed
        ## GRANS_CLUB_SINGLES/SYNS_PER_CLUBBED_SINGLE is the weight of single GABA synapse
        ## SYNS_PER_CLUBBED_SINGLE number of randomly delayed synapses are made
        self.set_xml_connections_by_multiplicity('granule_mitral_inh_singles',-1,'singles',\
                [("granule_mitral_GABA",g_m_GABA_weight,THRESHOLD,inh_synaptic_delay)],\
                inh_delay_distrib, granule_mitral_delay_spread)
        ## singles granule-spine/self-mitral inhibition is asynchronously delayed
        ## GRANS_CLUB_SINGLES/SYNS_PER_CLUBBED_SINGLE is the weight of single GABA synapse
        ## SYNS_PER_CLUBBED_SINGLE number of randomly delayed synapses are made
        self.set_xml_connections_by_multiplicity('granule_mitral_inh_spinesingles',0,'singles',\
                [("mitral_self_GABA",1.0,THRESHOLD,inh_synaptic_delay)],\
                inh_delay_distrib, granule_mitral_delay_spread)
        ## joints granule-mitral inhibition is asynchronously delayed
        ## GRANS_CLUB_JOINTS number of synapses are made
        self.set_xml_connections_by_multiplicity('granule_mitral_inh_joints',-1,'joints',\
                [("granule_mitral_GABA",g_m_GABA_weight,THRESHOLD,inh_synaptic_delay)],\
                inh_delay_distrib, granule_mitral_delay_spread)
        ## joints granule-spine/self-mitral inhibition is asynchronously delayed
        ## GRANS_CLUB_JOINTS number of synapses are made
        self.set_xml_connections_by_multiplicity('granule_mitral_inh_spinejoints',0,'joints',\
                [("mitral_self_GABA",1.0,THRESHOLD,inh_synaptic_delay)],\
                inh_delay_distrib, granule_mitral_delay_spread)
        ## For 2 mitrals there will not be any multi-s.
        if len(self.granuleMultisTable) != 0:
            ## multis granule-mitral inhibition is asynchronously delayed
            self.set_xml_connections_by_multiplicity('granule_mitral_inh_multis',-1,'multis',\
                    [("granule_mitral_GABA",g_m_GABA_weight,THRESHOLD,inh_synaptic_delay)],\
                    inh_delay_distrib, granule_mitral_delay_spread)
            ## multis granule-spine/self-mitral inhibition is asynchronously delayed
            self.set_xml_connections_by_multiplicity('granule_mitral_inh_spinemultis',0,'multis',\
                    [("mitral_self_GABA",1.0,THRESHOLD,inh_synaptic_delay)],\
                    inh_delay_distrib, granule_mitral_delay_spread)

        ################### Baseline input to each family of granules
        self.set_xml_granule_baseline('singles','files','granules_singles')
        self.set_xml_granule_baseline('joints','files','granules_joints')
        ## For 2 mitrals there will not be any multi-s.
        if len(self.granuleMultisTable) != 0:
            self.set_xml_granule_baseline('multis','files','granules_multis')
        
        ## For noise in mitral responses
        ################## Baseline inhibition to mitral cells
        if random_inh: # default False
            self.set_xml_mitral_baseline('singles','files','mitrals')
            self.set_xml_mitral_baseline('joints','files','mitrals')
            ## For 2 mitrals there will not be any multi-s.
            if len(self.granuleMultisTable) != 0:
                self.set_xml_mitral_baseline('multis','files','mitrals')

    def set_xml_connections_by_multiplicity(self,\
        projectionname, direction, multiplicity, synapselist,\
        delay_distrib, delay_spread):

        if direction == 1 or direction == 2: # forward excitation
            source = 'mitrals'
            target='granules_'+multiplicity
        elif direction == -1: # reverse inhibition
            source = 'granules_'+multiplicity
            target='mitrals'
        elif direction == 0: # spine-based granule-mitral inhibition hack-modeled as self-inhibition
            source = 'mitrals'
            target='mitrals'
        gran_num = 0
        conn_num = 0
        self.xmlwriter.set_projection(name=projectionname,source=source,target=target)
        for synapse in synapselist:
            self.xmlwriter.set_synapse_props(synapse_type=synapse[0],\
                weight=synapse[1],threshold=synapse[2],delay=synapse[3])
        ## This holds a histogram of the # of simultaneous exc syn-s from premitid to any granule
        occurences = collections.defaultdict(int)
        premitid = 0
        strong_singles = [0]*(NUM_GLOMS*2-2)
        if multiplicity == 'singles':
            ## gran = (xcont, ycont, mitindex, (mit_seg_ids,num_conns), '1', extra_weight)
            for gran in self.granuleSinglesTable:
                weight_fac = 1.0
                ## num_conns i.e. gran[5] is due to clubbing different singly-connected granules and
                ## simultaneous connections to same granule up to GRANS_CLUB_SINGLES number.
                ## But due to wide spread of dendrites,
                ## I expect number of simultaneous connections to same single to average ~1.
                ## direction: 2 = extra exc; 1 = main exc; 0 = self-inh; -1 = main inh
                if direction == 2: num_multisyns = int(gran[5]+0.5) # extra excitation rounded
                elif direction == 1: num_multisyns = 1
                else:
                    ## Since SYNS_PER_CLUBBED_SINGLE <> GRANS_CLUB_SINGLES, we do some massaging.
                    num_multisyns = SYNS_PER_CLUBBED_SINGLE
                    weight_fac = gran[3][1]/float(SYNS_PER_CLUBBED_SINGLE)
                if direction in [1,2]: firstcell_delay = exc_synaptic_delay
                else: firstcell_delay = inh_synaptic_delay
                pre_cell_id = int(gran[2]) # mit id
                pre_seg_id = int(gran[3][0]) # gran[3][0] is always mitral segment, of type str
                ## if mitral's segid in soma & prim and nearest sec dends (column), increase inh gran-|mit weight
                ## columns i.e. strong gran-|mit near soma of lat mits is not made for distal connections
                ## exc mit->gran weight is not increased within column by a check in set_connection_by_direction()
                if STRONG_SYNAPSES and PROXIMAL_CONNECTION and pre_cell_id not in [0,1] \
                    and pre_seg_id in proximal_mitral_segids \
                    and strong_singles[pre_cell_id-2]*GRANS_CLUB_SINGLES < frac_directed*mit_syns:
                    strongsynon = True
                    strong_singles[pre_cell_id-2] += 1
                else: strongsynon = False
                for multisyn in range(num_multisyns):
                    ## multisyns represent spikes due to many cells represented by this clubbed cell
                    ## don't spread the delay of the first cell's spike
                    if multisyn == 0: delay = firstcell_delay
                    else: delay = self.get_delay(delay_distrib,delay_spread)
                    synapselist_delayed = [ (syn[0],syn[1]*weight_fac,delay) for syn in synapselist ]
                    self.set_connection_by_direction(direction=direction,elementid=str(conn_num),\
                        pre_cell_id=str(gran[2]),post_cell_id=str(gran_num),\
                        pre_segment_id=gran[3][0],post_segment_id=gran[4],\
                        synapselist=synapselist_delayed,strongsynon=strongsynon)
                    conn_num += 1
                gran_num += 1
        elif multiplicity == 'joints':
            ## gran = (xcont,ycont, mit1index, (mit1_seg_id,conn_mit1), '1', extra_weight1,\
            ##      mit2index, (mit2_seg_id,conn_mit2), '1', extra_weight2)
            for gran in self.granuleJointsTable:
                ## mit1index=gran[2] < mit2index=gran[6] by construction: see club_granules() above
                ## if the mit2index is connected directedly to mit1index as per DIRECTED_CONNS,
                ## then make stronger all synapses (exc and inh) from both mitrals to/from this granule.
                if STRONG_SYNAPSES and DIRECTED_CONNS[gran[6]]==gran[2]: strongsynon = True
                else: strongsynon = False
                ## conn_mit1 (conn_mit2) is due to (1) clubbing different jointly connected granules and
                ## (2) simultaneous (multiple same mit) connections to same granule up to GRANS_CLUB_JOINTS number.
                ## So exc_weight *= conn_mit1/float(GRANS_CLUB_JOINTS), number of exc synapses should remain 1
                ## ALso inh_weight *= conn_mit1/float(GRANS_CLUB_JOINTS),
                ##   number of inh synapses should be GRANS_CLUB_JOINTS .
                ## Above scheme works well only if GRANS_CLUB_JOINTS ~ 1,
                ## if large like GRANS_CLUB_SINGLES, can't separate the two contributions here.
                ## direction: 2 = extra exc; 1 = main exc; 0 = self-inh; -1 = main inh
                conn_mit1 = gran[3][1]
                if direction == 2: num_multisyns = int(gran[5]+0.5) # extra excitation rounded
                elif direction == 1: num_multisyns = 1
                else: num_multisyns = GRANS_CLUB_JOINTS
                if direction != 2 : simultaneous_wt_fac = conn_mit1/float(GRANS_CLUB_JOINTS)
                else: simultaneous_wt_fac = 1
                ## histogram of number of connections of a given premitid to any granule
                ## gran[2] always < gran[6], hence works only for premitid=0
                if gran[2]==premitid: occurences[simultaneous_wt_fac] += 1
                if direction in [1,2]: firstcell_delay = exc_synaptic_delay
                else: firstcell_delay = inh_synaptic_delay
                for multisyn in range(num_multisyns):
                    ## multisyns represent spikes due to many cells represented by this clubbed cell
                    ## don't spread the delay of the first cell's spike
                    if multisyn == 0: delay = firstcell_delay
                    else: delay = self.get_delay(delay_distrib,delay_spread)
                    ## syn = (synapse_type,weight,delay)
                    synapselist_delayed = [ (syn[0],syn[1]*simultaneous_wt_fac,delay) for syn in synapselist ]
                    self.set_connection_by_direction(direction=direction,elementid=str(conn_num),\
                        pre_cell_id=str(gran[2]),post_cell_id=str(gran_num),\
                        pre_segment_id=gran[3][0],post_segment_id=gran[4],\
                        synapselist=synapselist_delayed, strongsynon=strongsynon)
                    conn_num += 1
                ## direction: 2 = extra exc; 1 = main exc; 0 = self-inh; -1 = main inh
                conn_mit2 = gran[7][1]
                if direction == 2: num_multisyns = int(gran[9]+0.5) # extra excitation rounded
                elif direction == 1: num_multisyns = 1
                else: num_multisyns = GRANS_CLUB_JOINTS
                if direction != 2 : simultaneous_wt_fac = conn_mit2/float(GRANS_CLUB_JOINTS)
                else: simultaneous_wt_fac = 1
                if direction in [1,2]: firstcell_delay = exc_synaptic_delay
                else: firstcell_delay = inh_synaptic_delay
                for multisyn in range(num_multisyns):
                    ## multisyns represent spikes due to many cells represented by this clubbed cell
                    ## don't spread the delay of the first cell's spike
                    if multisyn == 0: delay = firstcell_delay
                    else: delay = self.get_delay(delay_distrib,delay_spread)
                    synapselist_delayed = [ (syn[0],syn[1]*simultaneous_wt_fac,delay) for syn in synapselist ]
                    self.set_connection_by_direction(direction=direction,elementid=str(conn_num),\
                        pre_cell_id=str(gran[6]),post_cell_id=str(gran_num),\
                        pre_segment_id=gran[7][0],post_segment_id=gran[8],\
                        synapselist=synapselist_delayed, strongsynon=strongsynon)
                    conn_num += 1
                gran_num += 1
        elif multiplicity == 'multis':
            ## granuleMultisTable is an array of
            ## (x,y, granconn)
            for gran in self.granuleMultisTable:
                ## list of granconn-s in gran[2]
                ## granconn = (mitindex,(mit_segment_id,numconns),gran_segment_id,extra_exc_sisters)
                ## zip(* ) makes parallel arrays of mitindices, etc. & [0] selects mitindices array
                mitindices = zip(*gran[2])[0]
                for granconn in gran[2]:
                    ## direction: 2 = extra exc; 1 = main exc; 0 = self-inh; -1 = main inh
                    ## strengthen the exc and inh synapses if this multi mediates connections
                    ## between directedly connected cells
                    ## by default synapse is standard, unless strong-ified by conditions below
                    strongsynon = False
                    if STRONG_SYNAPSES:
                        current_mitindex = granconn[0]
                        ## if current mitral is one of the central sisters,
                        ## and this granule is also connected to
                        ## a mitral directedly connected to the sister, strong-ify the synapse.
                        if current_mitindex in [central_glom,central_glom+1]:
                            for mitindex in mitindices:
                                if DIRECTED_CONNS[mitindex]==current_mitindex:
                                    strongsynon = True
                        ## if current mitral is not a central sister,
                        ## and it is directedly connected to a sister,
                        ## which is connected to this granule, strong-ify the synapse.
                        elif DIRECTED_CONNS[current_mitindex] in [central_glom,central_glom+1] \
                             and DIRECTED_CONNS[current_mitindex] in mitindices:
                                strongsynon = True
                    ## Since multi-s are not clubbed, numconns are simultaneous.
                    ## Hence increase synaptic weight for exc & inh (syn[1]*numconns below),
                    ## but keep the number of synapses = 1 for the main exc and inhibition.
                    numconns = granconn[1][1]
                    if direction == 2: num_multisyns = int(granconn[3]+0.5) # extra excitation rounded
                    elif direction == 1: num_multisyns = 1
                    else: num_multisyns = 1
                    if direction != 2 : simultaneous_wt_fac = numconns
                    else: simultaneous_wt_fac = 1
                    if direction in [1,2]: firstcell_delay = exc_synaptic_delay
                    else: firstcell_delay = inh_synaptic_delay
                    ## histogram of number of connections of a given premitid to any granule
                    if granconn[0]==premitid: occurences[simultaneous_wt_fac] += 1
                    for multisyn in range(num_multisyns):
                        ## multisyns represent spikes due to many cells represented by this clubbed cell
                        ## don't spread the delay of the first cell's spike
                        if multisyn == 0: delay = firstcell_delay
                        else: delay = self.get_delay(delay_distrib,delay_spread)
                        ## multiple connections
                        ## syn = (synapse_type,weight,delay)
                        synapselist_delayed = [ (syn[0],syn[1]*simultaneous_wt_fac,delay) for syn in synapselist ]
                        self.set_connection_by_direction(direction=direction,elementid=str(conn_num),\
                            pre_cell_id=str(granconn[0]),post_cell_id=str(gran_num),\
                            pre_segment_id=granconn[1][0],post_segment_id=granconn[2],\
                            synapselist=synapselist_delayed, strongsynon=strongsynon)
                        conn_num += 1
                gran_num += 1
        print projectionname, 'direction =', direction, multiplicity
        print 'To mit',premitid,'histogram of simultaneous connections',occurences

    def set_connection_by_direction(self, direction, elementid, pre_cell_id,\
        post_cell_id, pre_segment_id, post_segment_id, synapselist, strongsynon=False):
        ## modify below to change directed inhibition mit->gran vis a vis directed excitation gran-|mit
        ## pre_cell_id is of type str!!!
        if strongsynon:
            if direction in [0,-1]: # self and lateral inhibition
                if STRONG_IHNSYNS_CENTRAL_ONLY:
                    ## only strengthen the central mit synapses...
                    if pre_cell_id in ['0','1']:
                        strongsynfactor = strongsynfactorinh
                    else: strongsynfactor = 1.0
                else:
                    ## strengthen central and lat synapses
                    strongsynfactor = strongsynfactorinh
            elif direction in [1,2]: # main and extra excitation
                #strongsynfactor = 1.0
                ### only strengthen the non-central mit synapses...
                #if pre_cell_id not in ['0','1']:
                ### only strengthen the mit->gran synapses that are distal from mit soma
                if int(pre_segment_id) not in proximal_mitral_segids: # pre_segment_id is of type str
                    strongsynfactor = strongsynfactorexclateral
                else: strongsynfactor = strongsynfactorexccentral
        else: strongsynfactor = 1.0
        ## This sets the randomized weights (lognormal or uniform or MAXSYNS)
        synproplist = self.make_synproplist(synapselist,strongsynfactor)
        if direction == 1 or direction == 2: # forward - mitral to granule
            self.xmlwriter.set_connection(elementid=elementid,\
                pre_cell_id=pre_cell_id,post_cell_id=post_cell_id,\
                pre_segment_id=pre_segment_id,post_segment_id=post_segment_id,\
                synproplist = synproplist)
        elif direction == -1: # reverse - granule to mitral
            ## reverse (reciprocal) synapse, here pre is post and vice versa compared to above forward
            synproplistdecayed = [ (syn[0],syn[1]*self.dend_decay(pre_segment_id),syn[2]) \
                for syn in synproplist]
            self.xmlwriter.set_connection(elementid=elementid,\
                pre_cell_id=post_cell_id,post_cell_id=pre_cell_id,\
                pre_segment_id=post_segment_id,post_segment_id=pre_segment_id,\
                synproplist = synproplistdecayed)
        if direction == 0: # self/spine inhibition - mitral to self
            synproplistdecayed = [ (syn[0],syn[1]*self.dend_decay(pre_segment_id),syn[2]) \
                for syn in synproplist]
            self.xmlwriter.set_connection(elementid=elementid,\
                pre_cell_id=pre_cell_id,post_cell_id=pre_cell_id,\
                pre_segment_id=pre_segment_id,post_segment_id=pre_segment_id,\
                synproplist = synproplistdecayed)

    def dend_decay(self,segid):
        """ Lowe 2002 found that mit--|gran inh synapse conductance density (per area)
        drops exponentially, differently for prim and sec dends. Hence an exponential decay,
        and a scaling with dendritic diameter. This returns a decay factor for a given segment. """
        for seg in self.granule_mitral_syn_positions:
            ## seg = (segid,(segx1,segy1,segz1),(segx2,segy2,segz2),seglength,segdia)
            ## exponential decay taking start of segment, not center of segment, else too huge decay.
            if segid == seg[0]:
                x=seg[1][0]
                y=seg[1][1]
                z=seg[1][2]
                segdia = seg[4]
                break
        ## all dia-s are normalized to dia of dendrite's base.
        ## soma is same as for primary dendrite base, and no exp decay.
        if segid == '0': # if soma
            return 1.0
        elif segid in ['15','16','17','18','19','20']: # if primary dendrite
            return segdia/self.mit_prim_decay_dia_base * math.exp(-z/mit_prim_decay_l)
        else: # secondary dendrite
            if USE_SECDEND_DECAY:
                ## Choose one of the below lines for receptor areal dendity decaying exp or constant.
                return segdia/self.mit_sec_decay_dia_base * math.exp(-sqrt(x**2+y**2+z**2)/mit_sec_decay_l)
                #return segdia/self.mit_sec_decay_dia_base
            else: return 1.0

    def set_xml_granule_baseline(self, multiplicity, source, target):
        if not MAX_SYNS: # if MAX_SYNS is True, set synaptic weights to 1.0+SYN_VARIATION times normal, else 1.0
            gmax_factor = 1.0 # This is an average synapse, so do not randomize it.
        else:
            gmax_factor = 1.0+SYN_VARIATION

        gran_num = 0
        self.xmlwriter.set_projection(name='granule_baseline_'+multiplicity,source=source,target=target)
        self.xmlwriter.set_synapse_props(synapse_type='mitral_granule_AMPA',\
            weight=self.AMPA_factor*gmax_factor*exc_avg_factor,threshold=0,delay=exc_synaptic_delay)
        self.xmlwriter.set_synapse_props(synapse_type='mitral_granule_NMDA',\
            weight=self.NMDA_factor*gmax_factor*exc_avg_factor,threshold=0,delay=exc_synaptic_delay)
        if multiplicity == 'singles':
            for gran in self.granuleSinglesTable:
                self.xmlwriter.set_connection(elementid=str(gran_num),\
                    pre_cell_id='file',post_cell_id=str(gran_num),\
                    pre_segment_id=str(int(uniform(num_gran_baseline_files))),post_segment_id='1')
                gran_num += 1
        elif multiplicity == 'joints':
            for gran in self.granuleJointsTable:
                self.xmlwriter.set_connection(elementid=str(gran_num),\
                    pre_cell_id='file',post_cell_id=str(gran_num),\
                    pre_segment_id=str(int(uniform(num_gran_baseline_files))),post_segment_id='1')
                gran_num += 1
        elif multiplicity == 'multis':
            for gran in self.granuleMultisTable:
                self.xmlwriter.set_connection(elementid=str(gran_num),\
                    pre_cell_id='file',post_cell_id=str(gran_num),\
                    pre_segment_id=str(int(uniform(num_gran_baseline_files))),post_segment_id='1')
                gran_num += 1

    def set_xml_mitral_baseline(self, multiplicity, source, target):
        """ For setting a recurrent inhibition to increase noise in the mitral cell. """
        syn_num = 0
        self.xmlwriter.set_projection(name='mitral_baseline_mit'+multiplicity[0],\
            source=source,target=target)
        self.xmlwriter.set_synapse_props(synapse_type='mitral_self_GABA',\
            weight=1.0,threshold=0,delay=exc_synaptic_delay)
        if multiplicity == 'singles':
            for gran in self.granuleSinglesTable:
                for synnum in range(GRANS_CLUB_SINGLES):
                    synproplist = self.make_synproplist([('mitral_self_GABA',1.0,exc_synaptic_delay)])
                    self.xmlwriter.set_connection(elementid=str(syn_num),\
                        pre_cell_id='file',post_cell_id=str(gran[2]),\
                        pre_segment_id=str(int(uniform(num_gran_baseline_files))),\
                        post_segment_id=str(gran[3]), synproplist=synproplist)
                    syn_num += 1
        elif multiplicity == 'joints':
            for gran in self.granuleJointsTable:
                for synnum in range(GRANS_CLUB_JOINTS):
                    synproplist = self.make_synproplist([('mitral_self_GABA',1.0,exc_synaptic_delay)])
                    self.xmlwriter.set_connection(elementid=str(syn_num),\
                        pre_cell_id='file',post_cell_id=str(gran[2]),\
                        pre_segment_id=str(int(uniform(num_gran_baseline_files))),\
                        post_segment_id=str(gran[3]), synproplist=synproplist)
                    syn_num += 1
                    synproplist = self.make_synproplist([('mitral_self_GABA',1.0,exc_synaptic_delay)])
                    self.xmlwriter.set_connection(elementid=str(syn_num),\
                        pre_cell_id='file',post_cell_id=str(gran[6]),\
                        pre_segment_id=str(int(uniform(num_gran_baseline_files))),\
                        post_segment_id=str(gran[7]), synproplist=synproplist)
                    syn_num += 1
        elif multiplicity == 'multis':
            for gran in self.granuleMultisTable:
                for syn in gran[2]:
                    synproplist = self.make_synproplist([('mitral_self_GABA',1.0,exc_synaptic_delay)])
                    self.xmlwriter.set_connection(elementid=str(syn_num),\
                        pre_cell_id='file',post_cell_id=str(syn[0]),\
                        pre_segment_id=str(int(uniform(num_gran_baseline_files))),\
                        post_segment_id=str(syn[1]), synproplist=synproplist)
                    syn_num += 1

    def set_xml_ORN_connections(self):
        ## ORN_mitral synapses
        self.xmlwriter.set_projection(name='ORN_mitral',\
            source='file',target='mitrals')
        self.xmlwriter.set_synapse_props(synapse_type='ORN_mitral',\
            weight=1,threshold=0,delay=ORN_synaptic_delay)
        num_ORN_mitral_syn_positions = len(self.ORN_mitral_syn_positions)
        for mitnum in range(len(self.mitral_positions)):
            mit_seg_files = {}
            for synnum in range(NUM_ORN_MITRAL_SYNS):
                ## choose random one out of ORN_mitral_syn_positions and take the segid [0]
                post_segment_id = self.ORN_mitral_syn_positions\
                    [int(uniform(0,num_ORN_mitral_syn_positions))][0]
                pre_segment_id = str(int(uniform(NUM_ORN_FILES_PER_GLOM)))
                if post_segment_id not in mit_seg_files: # key does not exist
                    mit_seg_files[post_segment_id] = [pre_segment_id]
                else:
                    mit_seg_files[post_segment_id].append(pre_segment_id)
            # write one connection per post_segment that has the list of files
            for i,post_segment_id in enumerate(mit_seg_files.keys()):
                # id's not sequential, they jump every mitnum
                self.xmlwriter.set_connection(elementid=str(mitnum*NUM_ORN_MITRAL_SYNS+i),\
                    pre_cell_id='file',post_cell_id=str(mitnum),\
                    pre_segment_id='_'.join(mit_seg_files[post_segment_id]),\
                    post_segment_id=post_segment_id)
        # ORN_PG synapses
        self.xmlwriter.set_projection(name='ORN_PG',source='file',target='PGs')
        self.xmlwriter.set_synapse_props(synapse_type='ORN_PG',weight=1,threshold=0,delay=ORN_synaptic_delay)
        num_ORN_PG_syn_positions = len(self.ORN_PG_syn_positions)
        for PGnum in range(self.PG_num):
            PG_seg_files = {}
            for synnum in range(NUM_ORN_PG_SYNS):
                # choose random one out of ORN_PG_syn_positions and take the segid [0]
                post_segment_id = self.ORN_PG_syn_positions[int(uniform(0,num_ORN_PG_syn_positions))][0]
                pre_segment_id=str(int(uniform(NUM_ORN_FILES_PER_GLOM)))
                if post_segment_id not in PG_seg_files: # key does not exist
                    PG_seg_files[post_segment_id] = [pre_segment_id]
                else:
                    PG_seg_files[post_segment_id].append(pre_segment_id)
            # write one connection per post_segment that has the list of files
            for i,post_segment_id in enumerate(PG_seg_files.keys()):
                # id's not sequential, they jump every PGnum
                self.xmlwriter.set_connection(elementid=str(PGnum*NUM_ORN_PG_SYNS+i),\
                    pre_cell_id='file',post_cell_id=str(PGnum),\
                    pre_segment_id='_'.join(PG_seg_files[post_segment_id]),\
                    post_segment_id=post_segment_id)

    def set_xml_SA_connections(self):
        # SA->PG excitatory synapses which feed in global 
        self.xmlwriter.set_projection(name='SA_PG',source='file',target='PGs')
        self.xmlwriter.set_synapse_props(synapse_type='SA_PG',weight=1,threshold=0,delay=SA_delay)
        num_SA_PG_syn_positions = len(self.SA_PG_syn_positions)
        for PGnum in range(self.PG_num):
            PG_seg_files = {}
            for synnum in range(NUM_SA_SYNS_PER_PG):
                # choose random one out of SA_PG_syn_positions and take the segid [0]
                post_segment_id = self.SA_PG_syn_positions[int(uniform(0,num_SA_PG_syn_positions))][0]
                pre_segment_id=str(int(uniform(NUM_ORN_FILES_PER_GLOM)))
                if post_segment_id not in PG_seg_files: # key does not exist
                    PG_seg_files[post_segment_id] = [pre_segment_id]
                else:
                    PG_seg_files[post_segment_id].append(pre_segment_id)
            # write one connection per post_segment that has the list of files
            for i,post_segment_id in enumerate(PG_seg_files.keys()):
                # id's not sequential, they jump every PGnum
                self.xmlwriter.set_connection(elementid=str(PGnum*NUM_SA_SYNS_PER_PG+i),\
                    pre_cell_id='file',post_cell_id=str(PGnum),\
                    pre_segment_id='_'.join(PG_seg_files[post_segment_id]),\
                    post_segment_id=post_segment_id)

    def rotate(self,x,y,ztheta):
        return (x*cos(ztheta)-y*sin(ztheta), x*sin(ztheta)+y*cos(ztheta))

    def get_nearest_mitral_syn_segment(self,x,y,z,mitnum=0,type='granule'):
        mitx,mity,mitzrot = self.mitral_positions[mitnum]
        mitz = self.mitral_z
        sqdistances = []
        close_ones = 0
        if 'granule' in type:
            cutoff_distsq = glom_r**2
            for seg in self.granule_mitral_syn_positions:
                ## seg = (segid,(segx1,segy1,segz1),(segx2,segy2,segz2),seglength,segdia)
                segid = seg[0]
                segx = (seg[1][0]+seg[2][0])/2.0
                segy = (seg[1][1]+seg[2][1])/2.0
                segx,segy = self.rotate(segx,segy,mitzrot)
                ## do not include z in the distance calculation since:
                ## 1) it was not used while assigning the granule connections,
                ## 2) if z is used, it creates inward pointing granule-mitral connections (towards mitral center).
                #distsq = ((x-segx-mitx)**2+(y-segy-mity)**2+(z-segz-mitz)**2)
                distsq = ((x-segx-mitx)**2+(y-segy-mity)**2)
                sqdistances.append( (distsq,segid) )
                if distsq<cutoff_distsq: close_ones += 1
        elif type == 'PG':
            ## my PG cell is flat in the x-y plane while the mitral tuft compartments are quite spread out in z
            ## if I include z distance, it dominates and only two segments are closest always.
            ## thus do not include z distance below.
            cutoff_distsq = (2*glom_r)**2
            for (segid,segx,segy,segz) in self.PG_mitral_syn_positions:
                segx,segy = self.rotate(segx,segy,mitzrot)
                distsq = (x-segx-mitx)**2+(y-segy-mity)**2
                sqdistances.append( (distsq,segid) )
                if distsq<cutoff_distsq: close_ones += 1
        ## If there are two or more mitral segments within cutoff distance of the granule,
        ## choose one among them randomly, else choose the closest segment.
        ## This is needed because primary dendrite segments will all be close in x,y distance (z not considered)
        ## For these primary dendrite segments, one needs to choose any of these close ones.
        sqdistances.sort()
        if close_ones >= 2:
            random_idx = int( uniform(0,close_ones) )
            chosen = sqdistances[random_idx]
        else:
            chosen = sqdistances[0]
        ## If PRIMDEND_FORCE then if DIRECTED and joint granule,
        ## force connection to a primary dendrite for mit0
        if PRIMDEND_FORCE and DIRECTED and ('joint' in type) and mitnum==0:
            chosen = ( 0.0, str(int(uniform(15,20.9999))) )
        return chosen

    def read_mitral_syn_segments(self):
        """
        Reads which mitral segments can form synapses with ORNs/granule/PG cells
        from the info given by the MorphML_reader.
        Also calculate total lateral dendritic length
        This also includes primary dendrite length since reciprocal synapses are present on it too.
        """
        self.granule_mitral_syn_positions = []
        self.PG_mitral_syn_positions = []
        self.ORN_mitral_syn_positions = []
        self.lateral_length = 0.0 # total lateral dendritic length
        self.memb_area = 0.0 # total membrane area on dendrite
        for segment in self.cellSegmentDict['mitral'].values():
            # segment = [ segname,(proximalx,proximaly,proximalz),\
            #    (distalx,distaly,distalz),diameter,length,[potential_syn1, ... ] ]
            ## segment[0] is segment name, which has segid at the end after an underscore
            segid = string.split(segment[0],'_')[-1]
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
            if "granule_mitral" in segment[5]: # is granule_mitral one of the potential synapses in this segment?
                segdia = segment[3]
                seglength = segment[4]
                self.granule_mitral_syn_positions.append( \
                    (segid,(segx1,segy1,segz1),(segx2,segy2,segz2),seglength,segdia))
                self.lateral_length += seglength
                self.memb_area += pi*segdia*seglength
                ## base dia of starting dendrites is needed for decay of syn wt with dia.
                if segid == '15': self.mit_prim_decay_dia_base = segdia # for prim dend
                if segid == '115': self.mit_sec_decay_dia_base = segdia # for sec dend
            if "PG_mitral" in segment[5]: # is PG_mitral one of the potential synapses in this segment?
                self.PG_mitral_syn_positions.append((segid,segx,segy,segz))
            if "ORN_mitral" in segment[5]: # is ORN_mitral one of the potential synapses in this segment?
                self.ORN_mitral_syn_positions.append((segid,segx,segy,segz))

    def read_PG_syn_segments(self):
        """
        Reads which PG segments can form synapses with ORNs or mitrals
        from the info given by the MorphML_reader.
        """
        self.mitral_PG_syn_positions = []
        self.ORN_PG_syn_positions = []
        self.SA_PG_syn_positions = []
        for segment in self.cellSegmentDict['PG'].values():
            # segment = [ segname,(proximalx,proximaly,proximalz),(distalx,distaly,distalz),\
            #    diameter,length,[potential_syn1, ... ] ]
            # segment[0] is segment name, which has segid at the end  after an underscore
            segid = string.split(segment[0],'_')[-1]
            # average segment position from the proximal and distal points
            segx = ( segment[1][0] + segment[2][0] ) / 2.0
            segy = ( segment[1][1] + segment[2][1] ) / 2.0
            segz = ( segment[1][2] + segment[2][2] ) / 2.0
            if "mitral_PG" in segment[5]: # is mitral_PG one of the potential synapses in this segment?
                self.mitral_PG_syn_positions.append((segid,segx,segy,segz))
            if "ORN_PG" in segment[5]: # is ORN_PG one of the potential synapses in this segment?
                self.ORN_PG_syn_positions.append((segid,segx,segy,segz))
            if "SA_PG" in segment[5]: # is ORN_PG one of the potential synapses in this segment?
                self.SA_PG_syn_positions.append((segid,segx,segy,segz))

    # this returns a random r distributed according to 1/r between soma_r and mit_dend_r
    # unlike Egger and Urban 2006, there are no synapses at mitral soma.
    def oneoverr_dist(self): # use http://en.wikipedia.org/wiki/Inverse_transform_sampling
        # synapses between soma_r and mit_dend_r are distributed as 1/r
        # int_soma_r^mit_dend_r k/r dr = 1
        k = 1 / (math.log(mit_dend_r) - math.log(soma_r))
        lnrs = math.log(soma_r)
        x = uniform(0.0,1.0)
        return math.exp(x/k + lnrs)
        
    def lineardrop_dist(self):
        ##### Christie et al 2001: Fig 4B: At 800 microns (i.e. roughly mit_dend_r),
        ##### lateral charge goes to zero. 100% at 0 microns. Linear inbetween.
        ##### See also my onenote labnotes of 8th July 2010
        y = uniform(0.0,1.0)
        return mit_dend_r*(1-sqrt(1-y))

    def pick_granules(self, syns, xmitral, ymitral, glomnum, mitnum, mitzrot, mitseg=()):
        for synnum in range(syns):
            ## This is a do-while that executes until a random granule is found
            ## that is not already connected to this mitral (if allow_multi_gran_conn is false)
            ## (if allow_multi_gran_conn is true (default): see below).
            while True:
                if ISOTROPIC:
                    #r = self.lineardrop_dist()
                    ## see onenote labnotes of 27 May 2010 for the reason for choosing 1/r over uniform r.
                    #r = self.oneoverr_dist()
                    ## Uniform over r rather than 1/r seems more physical and inbetween 1/r and uniform over area.
                    ## uniform over r and theta does not make it uniform over the circular area!
                    ## It still remains clumped towards the center
                    r = uniform(soma_r, mit_dend_r)
                    theta = uniform(0.0,2*math.pi)
                    x = xmitral + r*math.cos(theta)
                    y = ymitral + r*math.sin(theta)
                    mit_seg_id = self.get_nearest_mitral_syn_segment( \
                        x,y,self.granule_z,mitnum=mitnum,type='granule')[1]
                else:
                    ## mitseg = (segid,(segx1,segy1,segz1),(segx2,segy2,segz2),seglength,segdia)
                    mit_seg_id = mitseg[0]
                    xseg1, yseg1 = self.rotate(mitseg[1][0], mitseg[1][1], mitzrot)
                    xseg2, yseg2 = self.rotate(mitseg[2][0], mitseg[2][1], mitzrot)
                    x1 = xmitral + xseg1
                    y1 = ymitral + yseg1
                    x2 = xmitral + xseg2
                    y2 = ymitral + yseg2
                    ## pick a random synaptic location on this mitral segment
                    ## picking an x in (x1,x2) and y in (y1,y2) actually picks a point in a rectangle
                    ## to get a point on the segment, pick a parameter t, insert in parametric line
                    t = uniform(0.0,1.0)
                    x = x1 + t*(x2-x1)
                    y = y1 + t*(y2-y1)
                ## picking a uniform r and theta of a granule's dendritic field has a drop off with r.
                ## unlike mitral's 1/r drop off, I assume granule has uniform density over area.
                ## so choose uniform gran_dend_r around above x and y and assign to the nearest granule.
                x += uniform(-gran_dend_r,gran_dend_r)
                y += uniform(-gran_dend_r,gran_dend_r)
                xindex = int(round(x/gran_spacing))+self.gran_xmidindex
                yindex = int(round(y/gran_spacing))+self.gran_ymidindex
                if allow_multi_gran_conn: # connect to this gran even if this mit is connected before
                    break
                else: # only connect if not connected to this mit before
                    ## if this granule is not connected to this mitral already,
                    ## only then break out of loop and connect these two.
                    ## self.connArray[xindex][yindex][glomnum] is a list of (mitnum,segid)-s
                    ## zip(* ) makes parallel arrays of mitnum-s and segid-s & [0] selects mitnum-s array
                    mits_segs = self.connArray[xindex][yindex][glomnum]
                    if mits_segs is None: break
                    elif mitnum not in zip(*mits_segs)[0]:
                        break
            ## Below alert is ok when doing activity dependent inhibition sims at very large distances
            ## Else some bug in program.
            if xindex<0 or xindex>=self.gran_xmaxindex or yindex<0 or yindex>=self.gran_ymaxindex:
                print "Alert: Some mitral segments are outside granule bed range!"
                print "This alert is ok for ADI at mitA-mitB separations > ~1700microns, else bug."
                continue
            if self.granArray[xindex][yindex] is None:
                self.granArray[xindex][yindex] = 1
                self.num_grans += 1
            ## append (mitnum,segid)
            ## self.connArray[xindex][yindex][glomnum] holds an object:
            ## here a python list, not numpy array; hence use list.append()
            if self.connArray[xindex][yindex][glomnum] is None:
                self.connArray[xindex][yindex][glomnum] = [(mitnum,mit_seg_id)]
            else:
                self.connArray[xindex][yindex][glomnum].append((mitnum,mit_seg_id))

    def granuleConnect(self, xmitral, ymitral, glomnum, mitnum, mitzrot, frac_syns):
        """
        Connect given mitral to a randomly chosen granule cell
         at uniformly distributed r and theta from center of mitral.
        xmitral and ymitral are the x and y coordinates of mitral,
        mitzrot is rotation in degrees about z-axis;
        glomnum, mitnum are the indices for the mitral to be connected.
        mit_syns*frac_syns are connected either isotropically or uniformly along dendrites.
        If allow_multi_gran_conn is false,
            pick_granule avoids connecting a granule to the same mitral twice.
        """
        if ISOTROPIC:
            ## if isotropic, the nearest mit segs are chosen in pick_granules.
            self.pick_granules(int(frac_syns*mit_syns), xmitral, ymitral, glomnum, mitnum, mitzrot)
        else:
            for mitseg in self.granule_mitral_syn_positions:
                ## mitseg = (segid,(segx1,segy1,segz1),(segx2,segy2,segz2),seglength,segdia)
                ## numsyns in this segment is proportional to its length
                ## I add uniform(0,1) before taking int,
                ## so that total number of synapses is on average 'mit_syns*frac_syns', not less.
                #### Use below if-else to set 80 syns on soma, else about 20 get set up.
                if mitseg[0]=='0': ## soma
                    numsyns = 80
                else:
                    numsyns = int( mitseg[3]/self.lateral_length*mit_syns*frac_syns + uniform(0,1) )
                #### about 20 get set up on soma with below line
                #numsyns = int( mitseg[3]/self.lateral_length*mit_syns*frac_syns + uniform(0,1) )
                ## comment above and uncomment below for membrane area distribution instead of length
                #numsyns = int( pi*mitseg[3]*mitseg[4]/self.memb_area*mit_syns*frac_syns + uniform(0,1) )
                self.pick_granules(numsyns, xmitral, ymitral, glomnum, mitnum, mitzrot, mitseg)

    def granuleConnectDirected(self, centralmitnum, mitnum, frac_joints):
        """
        Connect mitral mitnum to a granule cell which is already
            connected to a central mitral at uniformly distributed r and theta.
            within glom_r of central mitral's soma (if proximal).
        If distal, roles of centralmit and mit are reversed.
        If allow_multi_gran_conn is false,
            it avoids connecting a granule to the same mitral twice.
        """
        if not PROXIMAL_CONNECTION:
            ## if connected distally, swap the roles of centralmitnum and mitnum
            centralmitnum, mitnum = mitnum, centralmitnum
        mitralx,mitraly,mitralzrotation = self.mitral_positions[centralmitnum]
        latmitralx,latmitraly,latmitralzrotation = self.mitral_positions[mitnum]
        if PROXIMAL_CONNECTION and not self.dend_thru_central_mit\
                (centralmitnum,latmitralx,latmitraly,latmitralzrotation,goodfactor=0.0):
            ## If B's dend is not going through A's soma even with good_factor=0.0,
            ## it means no overlap.
            ## Hence, do not connected directed granules between these, and return.
            print "No overlap between",mitnum,"and",centralmitnum,\
                " hence not connecting extra directed granules."
            return
        ## I must append modelled mitralnums as 0,1 in self.connArray[xindex][yindex][glomnum]
        ## Hence obtain mitral labels are (glom,mitnum) using div-mod.
        cglomnum = centralmitnum / MIT_SISTERS # integer division
        centralmitnum_modulo = centralmitnum % MIT_SISTERS # mod
        glomnum = mitnum / MIT_SISTERS # integer division
        mitnum_modulo = mitnum % MIT_SISTERS # mod
        for synnum in range(int(frac_joints*mit_syns)):
            numiters_to_nonjointgran = 0
            while True:
                ## connect around soma of centralmitnum within area which is size of a glomerulus (column)
                ## uniform r and theta do not give a uniform distribution over the circular area: use square
                #r = uniform(0.0, glom_r)
                #theta = uniform(0.0,2*math.pi)
                #x = mitralx + r*math.cos(theta)
                #y = mitraly + r*math.sin(theta)
                x = mitralx + uniform(-glom_r,glom_r)
                y = mitraly + uniform(-glom_r,glom_r)
                xindex = int(round(x/gran_spacing))+self.gran_xmidindex
                yindex = int(round(y/gran_spacing))+self.gran_ymidindex
                ########## I am no longer considering only granules connected to central mitral cell.
                ########## was doing that to ensure exactly 10K synapses to premitral,
                ########## but postmitral developed more. So might as well allow >10K syns on pre-mit too.
                ########## So skip the rest of the checking for connection to central mitral.
                #break

                ## self.connArray[xindex][yindex][glomnum] is a list of (mitnum,segid)-s
                ## mitnums are 0 to 50 belonging to glomnum
                ## zip(* ) converts to parallel lists of mitnum-s and segid-s; [0] selects mitnum-s list
                mits_segs = self.connArray[xindex][yindex][glomnum]
                c_mits_segs = self.connArray[xindex][yindex][cglomnum]
                ## check if centralmitnum is connected to this granule or not
                if c_mits_segs is None: c_mitnum_conn = False
                elif centralmitnum_modulo in zip(*c_mits_segs)[0]: c_mitnum_conn = True
                else: c_mitnum_conn = False
                ## NOT allowing multiple connections for directed connectivity,
                ## else lots of one-sided multiple connections are created,
                ## since I connect only to latmit and not to centralmit (see 3rd Sep 2012 notes).

                ## connect only if connected to c_mitnum and not connected to mitnum,
                ## (but if no such is found after many iterations, connect anyways -- below)
                ## check if mitnum is connected to this granule or not
                if mits_segs is None: mitnum_notconn = True
                elif mitnum_modulo not in zip(*mits_segs)[0]: mitnum_notconn = True
                else: mitnum_notconn = False
                ## if this granule is connected to central mitral centralmitnum already,
                ## and not to mitnum; then break out of loop and connect it further to mitnum.
                ## of course if distally connected, the roles are reversed as done above
                if mitnum_notconn and c_mitnum_conn:
                    break
                ## if too many iterations pass without getting 'cmit but not latmit connected' gran,
                ## then just connect to a new granule, if allow_multi_gran_conn is True.
                elif numiters_to_nonjointgran>10000 and allow_multi_gran_conn:
                    ## connect central mitnum to this granule, choosing nearest segment of mitnum (not mitnum_modulo!!!)
                    c_mit_seg_id = self.get_nearest_mitral_syn_segment(x,y,self.granule_z,\
                            mitnum=centralmitnum,type='granule')[1]
                    ## self.connArray[xindex][yindex][glomnum] holds an object:
                    ## here a python list, not numpy array; hence use list.append()
                    if self.connArray[xindex][yindex][cglomnum] is None:
                        self.connArray[xindex][yindex][cglomnum] = [(centralmitnum_modulo,c_mit_seg_id)]
                    else: self.connArray[xindex][yindex][cglomnum].append((centralmitnum_modulo,c_mit_seg_id))
                    break
                numiters_to_nonjointgran += 1

            ## connect mitnum to this granule, choosing nearest segment of mitnum (not mitnum_modulo!!!)
            mit_seg_id = self.get_nearest_mitral_syn_segment(x,y,self.granule_z,\
                    mitnum=mitnum,type='granule')[1]
            ## self.connArray[xindex][yindex][glomnum] holds an object:
            ## here a python list, not numpy array; hence use list.append()
            if self.connArray[xindex][yindex][glomnum] is None:
                self.connArray[xindex][yindex][glomnum] = [(mitnum_modulo,mit_seg_id)]
            else: self.connArray[xindex][yindex][glomnum].append((mitnum_modulo,mit_seg_id))

    def PGConnect(self, glomx, glomy, glomnum, mitralnum):
        """
        Just connect the present mitral to any PG within this glomerulus
         - assumes that a mitral tuft covers the whole glomerulus.
        This doesn't matter since in my model,
         a glomerulus is a single entity receiving an average odor excitation.
        Unlike for granules, a mitral can connect to the same PG twice!
        There are no intra-glomerular differentials by design.
        """
        mitnum = glomnum*2+mitralnum
        for synnum in range(PGSYNS_PER_MITRAL):
            synx = uniform(-PG_halfextent,PG_halfextent)
            syny = uniform(-PG_halfextent,PG_halfextent)
            xidx = int( (synx + PG_halfextent)/PG_spacing )
            yidx = int( (syny + PG_halfextent)/PG_spacing )
            synx += glomx
            syny += glomy
            ## Don't choose the nearest mitral segment, a whole side of the mitral tuft gets left out
            #mitsegid = self.get_nearest_mitral_syn_segment(synx,syny,self.PG_z,mitnum,'PG')[1]
            ## Choose a mitral segment randomly from the tuft
            mitsegid = self.PG_mitral_syn_positions[int(uniform(0,len(self.PG_mitral_syn_positions)))][0]
            # randomly choose one out of the potential synaptic segments of the PG:
            PGsegid = self.mitral_PG_syn_positions[int(uniform(0,len(self.mitral_PG_syn_positions)))][0]
            ## self.PGconnArray[xindex][yindex][glomnum] holds an object:
            ## here a python list, not numpy array; hence use list.append()
            if self.PGConnArray[glomnum][xidx][yidx] is None:
                self.PGConnArray[glomnum][xidx][yidx] = [(mitnum,mitsegid,PGsegid)]
            else:
                self.PGConnArray[glomnum][xidx][yidx].append((mitnum,mitsegid,PGsegid))

if __name__ == "__main__":
    ### Seed only if called directly, else do not seed.
    ### Also seeding this way ensures seeding after importing other files that may set seeds.
    ### Thus this seed overrides other seeds.
    if manyprocs: seed([stim_net_seed*mpirank])
    else: seed([stim_net_seed]) ##### Seed numpy's random number generator. If no parameter is given, it uses current system time
    gen = ConnGenerator()
