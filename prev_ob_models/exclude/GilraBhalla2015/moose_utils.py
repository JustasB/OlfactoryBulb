import string, os
import moose
from moose.utils import *

## mpi_uniqueid could be for example: 'pulses'+str(mpirank)
def attach_spikes(filebase, timetable, mpi_uniqueid):
    ## read the file that contains all the ORN firing times for this glom, odor and avgnum
    filehandle = open(filebase+'.txt','r')
    spiketimelists = filehandle.readlines()
    filehandle.close()

    filenums = string.split(timetable.getField('fileNumbers'),'_')
    ## Merge all the filenums into a temp file, load it and delete it.
    spiketimes = []
    for filenum in filenums: # loop through file numbers
        timesstr = spiketimelists[int(filenum)]
        if timesstr != '\n':
            timestrlist = string.split(timesstr,' ')
            ## convert to float for sorting else '10.0'<'6.0'
            spiketimes.extend([float(timestr) for timestr in timestrlist])
    spiketimes.sort()
    ## ensure that different processes do not write to the same file by using mpi_uniqueid
    ## mpi_uniqueid could be for example: 'pulses'+str(mpirank)
    fn = os.getenv('HOME')+'/tempspikes_'+str(mpi_uniqueid)+'.txt'
    filehandle = open(fn,'w')
    filehandle.write('\n'.join([str(spiketime) for spiketime in spiketimes]))
    filehandle.close()

    ############# OB model specific hack to give all ORN inputs to tuft-base compartment
    ### tt_path = postcomp.path+'/'+syn_name_full+glomstr+'_tt' ## glomstr is '' for us
    #tt_split = timetable.path.split('/')
    #if tt_split[-1]=='ORN_mitral_tt':
    #    tt_path = string.join(tt_split[:-2],'/')+'/Seg0_glom_1_22/ORN_mitral_tt_'+tt_split[-2] # unique timetable
    #    ## Choose one of the below two
    #    #syn_path = string.join(tt_split[:-2],'/')+'/Seg0_glom_1_22/ORN_mitral' # connect new tt to glom[1] synapse
    #    syn_path = string.join(tt_split[:-1],'/')+'/ORN_mitral' # connect new tt to original synapse
    #    syn = moose.SynChan(syn_path) # wrapping created synapse
    #    tt = moose.TimeTable(tt_path) # new timetable
    #    # Be careful to connect the timetable only once while creating it as below:
    #    tt.connect("event", syn, "synapse")
    #    print "Connecting",timetable.path,"to",syn_path,"via",tt_path
    #    timetable = tt

    timetable.filename = fn
    os.remove(fn)
