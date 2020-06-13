from neuron import h, gui

from net_mitral_centric import mkgranule
import granules

gcs = []
for i in xrange(1): #All gs are same
    print('building ',i)
    gc = mkgranule(granules.gid2pos.keys()[i])
    gcs.append(gc)

from hoc2swc import neuron2swc, swc_type_from_section_name

def new_map(section_name):
    '''
    Returns an integer string of an SWC point type in response to a string name of a NEURON section.

    Override this method to map custom section names to parts of SWC cell.

    See column 2 of http://www.neuronland.org/NLMorphologyConverter/MorphologyFormats/SWC/Spec.html

    :param section_name: name string of a NEURON section
    :return: integer string e.g. "1" or "3" that corresponds to a SWC point type
    '''

    if "tuft" in section_name or "pri" in section_name:
        return "4"

    if "den" in section_name:
        if "apical" in section_name:
            return "4"

        return "3"

    if "axon" in section_name or "hillock" in section_name or "initial" in section_name:
        return "2"

    if "soma" in section_name:
        return "1"

    return "5"

swc_type_from_section_name.__code__ = new_map.__code__

neuron2swc("gc_morphology/gc.swc")