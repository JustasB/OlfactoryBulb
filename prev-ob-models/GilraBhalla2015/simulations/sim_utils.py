import sys
sys.path.extend(["..","../networks","../generators","../simulations"])

from moose_utils import * # imports moose, and 'moose.utils import *'
from networkConstants import *
from stimuliConstants import *
from data_utils import *

def setupTables(network, nopgs, nosingles, nojoints, nomultis, args={}, spikes=False):
    pgtables = []
    singletables = []
    jointtables = []
    multitables = []
    ## if spikes == True i.e. only spiketimes are recorded,
    ## then record all cells, else less number of cells.
    allcells = spikes
    if not nopgs:
        ## Setup the tables to pull data from PGs
        ## get cells connected to mitrals specified in args. if not specified get arbitrary.
        PGlist = getCellsByMitralConnection(args, network, 'PG_mitral', 'PGs', allcells)
        for i,PG in enumerate(PGlist):
            # assumes at least one soma and takes the first!
            PG.soma = moose.Compartment(get_matching_children(PG, ['Soma','soma'])[0])
            PG._vmTableSoma = setupTable("vmTableSoma",PG.soma,'Vm')
            pgtables.append(PG._vmTableSoma)
            if spikes:
                PG._vmTableSoma.stepMode = TAB_SPIKE
                PG._vmTableSoma.stepSize = THRESHOLD
    if not nosingles:
        ## Setup the tables to pull data from single granules
        ## get cells connected to mitrals specified in args. if not specified get arbitrary.
        singlesList = getCellsByMitralConnection(args, network,
            'granule_mitral_inh_singles', 'granules_singles', allcells)
        for i,gran in enumerate(singlesList):
            ## assumes at least one soma and takes the first!
            gran.soma = moose.Compartment(get_matching_children(gran, ['Soma','soma'])[0])
            gran._vmTableSoma = setupTable("vmTableSoma",gran.soma,'Vm')
            singletables.append(gran._vmTableSoma)
            if spikes:
                gran._vmTableSoma.stepMode = TAB_SPIKE
                gran._vmTableSoma.stepSize = THRESHOLD
    if not nojoints:
        ## Setup the tables to pull data from joint granules
        ## get cells connected to mitrals specified in args.
        ## if not specified get arbitrary.
        jointsList = getCellsByMitralConnection(args, network,
            'granule_mitral_inh_joints', 'granules_joints', allcells)
        for i,gran in enumerate(jointsList):
            ## assumes at least one soma and takes the first!
            gran.soma = moose.Compartment(get_matching_children(gran, ['Soma','soma'])[0])
            gran._vmTableSoma = setupTable("vmTableSoma",gran.soma,'Vm')
            jointtables.append(gran._vmTableSoma)
            if spikes:
                gran._vmTableSoma.stepMode = TAB_SPIKE
                gran._vmTableSoma.stepSize = THRESHOLD
    if not nomultis:
        ## Setup the tables to pull data from multi granules
        ## get cells connected to mitrals specified in args.
        ## if not specified get arbitrary.
        multisList = getCellsByMitralConnection(args, network,
            'granule_mitral_inh_multis', 'granules_multis', allcells)
        for i,gran in enumerate(multisList):
            ## assumes at least one soma and takes the first!
            gran.soma = moose.Compartment(get_matching_children(gran, ['Soma','soma'])[0])
            gran._vmTableSoma = setupTable("vmTableSoma",gran.soma,'Vm')
            multitables.append(gran._vmTableSoma)
            if spikes:
                gran._vmTableSoma.stepMode = TAB_SPIKE
                gran._vmTableSoma.stepSize = THRESHOLD
    return (pgtables, singletables, jointtables, multitables)

def exportTable(network, neuronProj, neuronPop, colours, \
            args={}, spikes=True, allcells=True):
    """ Return tables of cells in neuronPop population name
        connected to mitrals specified in args,
        via neuronProj projection name
        if args is not specified get all.
        NOTE: All tables for all cells' somas should be present if allcells is True.
        If allcells is False, cells called here should have tables present.
    """
    exportDict = {'spikes':spikes,'data_tables':[]}
    if array(colours).shape == (3,): coloursList = False
    else: coloursList = True
    ## get cells connected to mitrals specified in args. if not specified get arbitrary.    
    celllist = getCellsByMitralConnection(args, network, neuronProj, neuronPop, allcells)
    for i,cell in enumerate(celllist):
        cellname = cell.path.split('/')[-1]
        ## assumes at least one soma and takes the first!
        cell.soma = moose.Compartment(get_matching_children(cell,['Soma','soma'])[0])
        cellTablePath = cell.soma.path+"/data/vmTableSoma"
        if moose.context.exists(cellTablePath):
            cell._vmTableSoma = moose.Table(cellTablePath)
        else:
            print "SimError: Did not find "+cellTablePath
            sys.exit(1)
        if coloursList: colour=colours[i]
        else: colour=colours
        exportDict['data_tables'].append((cellname,colour,array(cell._vmTableSoma)))
    return exportDict

def exportTables(network, nopgs, nosingles, nojoints, nomultis, mitspikes, args, colours):
    ## mitrals first each with a different colour
    allTables = [ exportTable(network, '', 'mitrals', colours, \
            args=args, spikes=mitspikes, allcells=True) ]
    ## Each colour-entry below is a tuple of (baseline/initial colour, spiking/peak colour, colourmap)
    ## Each colour is a tuple of (r,g,b,a)
    ## Set in Moogli's config file, whether to change color, and/or alpha, or use colourmap.
    if not nopgs:
        allTables.append( exportTable(network,\
            'PG_mitral', 'PGs', ((0.3,0,0,0.3),(1,0,0,1),'jet'), args) )
    if not nosingles:
        allTables.append( exportTable(network,\
            'granule_mitral_inh_singles', 'granules_singles', \
            ((0.3,0,0.3,0.3),(1,0,1,1),'jet'), args) )
    if not nojoints:
        allTables.append( exportTable(network,\
            'granule_mitral_inh_joints', 'granules_joints', \
            ((0.3,0.3,0,0.3),(1,1,0,1),'jet'), args) )
    if not nomultis:
        allTables.append( exportTable(network,\
            'granule_mitral_inh_multis', 'granules_multis', \
            ((0,0.3,0.3,0.3),(0,1,1,1),'jet'), args) )
    return allTables

def getCellsByMitralConnection(args, network, projection, population, allcells=False):
    """
    Returns a list of pre-synaptic MOOSE Cell-s in 'projection' that are connected to mitrals
    specificied in args={'mitrals':[mitid1, mitid2, ...]}.
    If args does not have 'mitrals' key, an arbitrary list of cells from 'population' are returned.
    if allcells is False, a few cells are returned, else all cells are returned.
    """
    cellList = []
    cellUniques = []
    if args.has_key('mitrals'):
        for mitid in args['mitrals']:
            mitpath = 'mitrals_'+str(mitid)
            cellnum = 0
            if projection in network.projectionDict:
                for conn in network.projectionDict[projection][2]:
                    if mitpath in conn[2]: # if mitrals_<mitid> is substring of post_seg_path of this connection
                        cellpath = string.split(conn[1],'/')[1] # take out cellname from '/cellname/segmentname
                        ## Take only those cells that have not been taken before.
                        if cellpath not in cellUniques:
                            cell = moose.Cell(cellpath)
                            cellList.append(cell)
                            cellUniques.append(cellpath)
                            cellnum += 1
                            if not allcells and cellnum == 30: break
    else:
        if allcells: cellList = network.populationDict[population][1].values()
        else: cellList = network.populationDict[population][1].values()[0:60]
    return cellList

def exportMitralTables(mitMOOSETables,mitcolours,mitspikes):
    mitTables = []
    for miti,(mitname,mitCellTables) in enumerate(mitMOOSETables):
        mitTables.append((mitname,mitspikes,mitcolours[miti],\
            [(compname,array(compTable)) for compname,compTable in mitCellTables]))
    return mitTables

def setupMitralTables(network,mitspikes):
    mitTables = []
    for mit in network.populationDict['mitrals'][1].values():
        mitTables.append((mit.name,setupCellTables(mit,mitspikes)))
    return mitTables

def setupCellTables(cellObj,spikes):
    cellTables = []
    for compartmentid in cellObj.getChildren(cellObj.id): # compartments
        comp = moose.Compartment(compartmentid)
        cellTable = setupTable('VmTableExtra',comp,'Vm')
        if spikes:
            cellTable.stepMode = TAB_SPIKE
            cellTable.stepSize = THRESHOLD
        cellTables.append((comp.name,cellTable))
    return cellTables

def numpy_convert_tables(tables):
    return [[array(table) for table in onetypetables] for onetypetables in tables]

def plot_extras(timevec, tables, nopgs, nosingles, nojoints, nomultis, plottitle=''):
    if not nopgs:
        figure()
        title("PGs "+plottitle)
        for pgtable in tables[0]:
            plot(timevec[:len(pgtable)], pgtable, ',')
    if not nosingles:
        figure()
        title("single granules "+plottitle)
        for singletable in tables[1]:
            plot(timevec[:len(singletable)], singletable, ',')
    if not nojoints:
        figure()
        title("joint granules "+plottitle)
        for jointtable in tables[2]:
            plot(timevec[:len(jointtable)], jointtable, ',')
    if not nomultis:
        figure()
        title("multi granules "+plottitle)
        for multitable in tables[3]:
            plot(timevec[:len(multitable)], multitable, ',')

def plot_extras_spikes(binvec, tables, nopgs, nosingles,\
    nojoints, nomultis, bins, runt, settlet, plottitle=''):
    if not nopgs:
        figure()
        title("PGs "+plottitle)
        for pgtable in tables[0]:
            plot(binvec, plotBins(pgtable, bins, runt, settlet), '.-')
        axes_labels(gca(),"time (s)","rate (Hz)",fontsize=label_fontsize+2)
    if not nosingles:
        figure()
        title("single granules "+plottitle)
        for singletable in tables[1]:
            plot(binvec, plotBins(singletable, bins, runt, settlet), '.-')
        axes_labels(gca(),"time (s)","rate (Hz)",fontsize=label_fontsize+2)
    if not nojoints:
        figure()
        title("joint granules "+plottitle)
        for jointtable in tables[2]:
            plot(binvec, plotBins(jointtable, bins, runt, settlet), '.-')
        axes_labels(gca(),"time (s)","rate (Hz)",fontsize=label_fontsize+2)
    if not nomultis:
        figure()
        title("multi granules "+plottitle)
        for multitable in tables[3]:
            plot(binvec, plotBins(multitable, bins, runt, settlet), '.-')
        axes_labels(gca(),"time (s)","rate (Hz)",fontsize=label_fontsize+2)

def print_extras_activity(tables, nopgs, nosingles, nojoints, nomultis, contextstr):
    spikestable = {}
    if not nopgs:
        numspikes = 0
        numcells = 0
        for totcells,pgtable in enumerate(tables[0]):
            pgtable = array(pgtable)
            ## MOOSE often inserts one or two spiketime = 0.0 entries
            ## when storing spikes, so discount those
            numspikescell = len(where(pgtable>0.0)[0])
            if numspikescell>0:
                numspikes += numspikescell
                numcells += 1
        print "The number of PG cells spiking for",contextstr,"is",\
            numcells,"of",totcells+1,"; number of spikes =",numspikes
        pgstable = (numcells,totcells,numspikes)
        spikestable['PGs'] = pgstable
    if not nosingles:
        numspikes = 0
        numcells = 0
        for totcells,singletable in enumerate(tables[1]):
            singletable = array(singletable)
            ## MOOSE often inserts one or two spiketime = 0.0 entries
            ## when storing spikes, so discount those
            numspikescell = len(where(singletable>0.0)[0])
            if numspikescell>0:
                numspikes += numspikescell
                numcells += 1
        print "The number of single granule cells spiking for",contextstr,"is",\
            numcells,"of",totcells+1,"; number of spikes =",numspikes
        singlestable = (numcells,totcells,numspikes)
        spikestable['singles'] = singlestable
    if not nojoints:
        numspikes = 0
        numcells = 0
        for totcells,jointtable in enumerate(tables[2]):
            jointtable = array(jointtable)
            ## MOOSE often inserts one or two spiketime = 0.0 entries
            ## when storing spikes, so discount those
            numspikescell = len(where(jointtable>0.0)[0])
            if numspikescell>0:
                numspikes += numspikescell
                numcells += 1
        print "The number of joint granule cells spiking for",contextstr,"is",\
            numcells,"of",totcells+1,"; number of spikes =",numspikes
        jointstable = (numcells,totcells,numspikes)
        spikestable['joints'] = jointstable
    if not nomultis:
        numspikes = 0
        numcells = 0
        for totcells,multitable in enumerate(tables[3]):
            multitable = array(multitable)
            ## MOOSE often inserts one or two spiketime = 0.0 entries
            ## when storing spikes, so discount those
            numspikescell = len(where(multitable>0.0)[0])
            if numspikescell>0:
                numspikes += numspikescell
                numcells += 1
        print "The number of multi granule cells spiking for",contextstr,"is",\
            numcells,"of",totcells+1,"; number of spikes =",numspikes
        multistable = (numcells,totcells,numspikes)
        spikestable['multis'] = multistable
    return spikestable

## non-overlapping bins
def rebin_pulses(mitral_responses_list, numbins, runtime, settletime):
    numtrials = len(mitral_responses_list)
    numtrains = len(mitral_responses_list[0])
    nummits = len(mitral_responses_list[0][0])
    return [[[ plotBins( mitral_responses_list[trialnum][trainnum][mitnum],\
        numbins, runtime, settletime) \
        for mitnum in range(nummits)] \
            for trainnum in range(numtrains)] \
               for trialnum in range(numtrials)]

### overlapping bins
#def rebin_pulses(mitral_responses_list, numbins, runtime, settletime, bin_width_time):
#    numtrials = len(mitral_responses_list)
#    numtrains = len(mitral_responses_list[0])
#    nummits = len(mitral_responses_list[0][0])
#    return [[[ plotOverlappingBins( mitral_responses_list[trialnum][trainnum][mitnum],\
#        numbins, runtime, settletime, bin_width_time ) \
#        for mitnum in range(nummits)] \
#            for trainnum in range(numtrains)] \
#                for trialnum in range(numtrials)]


def rebin(mitral_responses_list, numbins, bin_width_time, numresps=1):
    numtrials = len(mitral_responses_list)
    numodors = len(mitral_responses_list[0])
    nummits = len(mitral_responses_list[0][0])
    ## NUM_RESPS is number of RESPIRATION cycles simulated, we rebin and return the last 'numresps'
    return [[[  plotOverlappingBins( mitral_responses_list[trialnum][odornum][mitnum],\
        numbins, RESPIRATION*numresps, SETTLETIME+(NUM_RESPS-numresps)*RESPIRATION, bin_width_time )
        for mitnum in range(nummits)] \
            for odornum in range(numodors)] \
                for trialnum in range(numtrials)]

def build_tweaks(mitralsclub, nospineinh, nosingles,
    nojoints, nomultis, nopgs, onlytwomits, 
    includeProjections=[], twomitrals=(0,2), nolateral=False):
    """
    Excludes singles/joints/multis granules and their projections as appropriate
    Excludes extra-excitation from unmodeled sisters as appropriate
    Include only the two required mitrals (and its connected granules) if ONLY_TWO_MITS is True
    """
    excludePopulations = []
    excludeProjections = ['SA']
    ## In odor_pulses, odor_morphs, scaled_pulses, I have not specified to include 
    ## file-based inputs to 2nd order cells as below. If not specified, force include:
    if 'granule_baseline' not in includeProjections: includeProjections.append('granule_baseline')
    if 'ORN_PG' not in includeProjections: includeProjections.append('ORN_PG')
    if not mitralsclub:
        excludeProjections.append('mitral_granule_extra_exc')
    if nospineinh:
        excludeProjections.append('_spinesingles')
        excludeProjections.append('_spinejoints')
        excludeProjections.append('_spinemultis')
    if nosingles:
        excludePopulations.append('singles')
        excludeProjections.append('_singles') # _ to avoid deleting spinesingles
    if nojoints:
        excludePopulations.append('joints')
        excludeProjections.append('_joints') # _ to avoid deleting spinejoints
    if nomultis:
        excludePopulations.append('multis')
        excludeProjections.append('_multis') # _ to avoid deleting spinemultis
    if nopgs:
        excludePopulations.append('PGs')
        excludeProjections.append('PG')
    if onlytwomits:
        onlyInclude = {'includePopulation':('mitrals',[str(twomitrals[0]),str(twomitrals[1])]),
            'includeProjections':includeProjections}
        return {'excludePopulations':excludePopulations,
            'excludeProjections':excludeProjections,'onlyInclude':onlyInclude}
    else:
        if nolateral:
            ## remove other mitrals so that there is no lateral inhibition
            ## differs from nojoints, in keeping the joints self-inhibition
            print "EXCLUDING OTHER MITS, KEEPING ONLY mits 0 and 1"
            onlyInclude = {'includePopulation':('mitrals',['0','1']),
                'includeProjections':includeProjections}
            return {'excludePopulations':excludePopulations,
                'excludeProjections':excludeProjections,'onlyInclude':onlyInclude}
        else:
            return {'excludePopulations':excludePopulations,\
                'excludeProjections':excludeProjections}
