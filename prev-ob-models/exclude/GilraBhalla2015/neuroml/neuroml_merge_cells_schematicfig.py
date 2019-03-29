## not using minidom as ElementTree is faster by order of magnitude or more.
#import xml.dom.ext
#import xml.dom.minidom

from xml.etree import ElementTree as ET
import sys, string
neuromlfilename = sys.argv[1]

sys.path.extend(["..","../neuroml","../simulations","../synapses"])
from moose.neuroml.neuroml_utils import *
from sim_utils import *

TWEAK = True
NO_SPINE_INH = True
NO_SINGLES = False
NO_JOINTS = False
NO_MULTIS = False
NO_PGS = False

## If unmodeled mitrals connect to same granule as modeled mitrals,
## provide extra excitation to granule from modeled mitral
## as compensation for unmodelled mitrals!
CLUB_MITRALS = True
ONLY_TWO_MITS = False # Only two mitrals to test lateral inhibition.
## set to 1 for 2MITS / 2 for 2GLOMS option set in generate_neuroML.py
mitralsidekickidx = 2

## mits to include is ignored unless only_two_mits is true
if TWEAK:
    tweaks = build_tweaks(CLUB_MITRALS, NO_SPINE_INH,\
        NO_SINGLES, NO_JOINTS, NO_MULTIS, NO_PGS,\
        ONLY_TWO_MITS, [0,2,4], mitralsidekickidx)

## load neuromlfile (essentially networkml, but root tag is neuroml)
print "loading",neuromlfilename
## create an ElementTree instance from an XML file
neuromldoc = ET.parse(neuromlfilename)
neuromlroot = neuromldoc.getroot()

## ElementTree.parse() replaces original namespaces shorthands with its own ns0, ns1, etc.
## but in moose.neuroml.utils.py, I have defined _namespace_map
## further, the root <neuroml> tag needs to get the namespace defining attributes:

## Somehow, with the new cElementTree.parse(), the neuromlroot does not have namespace attribs,
## below function then adds those attribs, bit while write(), cElementTree again adds namespace attribs,
## causing duplicates and later xslt proc errors.
## So just don't add these namespaces.

#set_neuroml_namespaces_attribs(neuromlroot)

if TWEAK:
    ## tweak the model
    print "Tweaking model ... "
    print tweaks
    tweak_model(neuromlroot, tweaks)
######## minidom
## load neuromlfile (essentially networkml, but root tag is neuroml)
#neuromldoc = xml.dom.minidom.parse(neuromlfilename)
#neuromlroot = neuromldoc.getElementsByTagName("neuroml")[0]

cells = None

def import_cell_file(name, filename):
    """method to include external descriptions of cell"""
    global cells, neuromlroot, neuromldoc # needed to modify global variables
    if cells is None:
        cells = ET.Element('{'+neuroml_ns+'}'+"cells")
        neuromlroot.insert(0,cells) # don't append cells - xml documents have a specified order.
        ###### minidom
        #cells = neuromldoc.createElement("cells")
        #cells.setAttribute("xmlns",neuroml_ns)
        ## WRONG - don't append cells - xml documents have a specified order.
        #neuromlroot.appendChild(cells)
    celldoc = ET.parse(filename)
    for cell in celldoc.findall(".//{"+neuroml_ns+"}cell"):
        if cell.attrib['name'] == name:
            cells.append(cell)
    ###### minidom
    #celldoc = xml.dom.minidom.parse(filename)
    #for cell in celldoc.getElementsByTagName("cell"):
    #    if cell.getAttribute("name")==name:
    #        cells.appendChild(cell)

## add cells
print "adding cells"
mitralfilename = "../cells/mitral_bbmit1993davisonMS_neuroML_L1_L2_L3.xml"
granulefilename = "../cells/granule_granadityaMS2007_neuroML_L1_L2_L3.xml"
PGfilename = "../cells/PG_aditya2010_neuroML_L1_L2_L3.xml"
import_cell_file("mitral", mitralfilename)
import_cell_file("granule", granulefilename)
import_cell_file("PG", PGfilename)

## write neuroml file with cells
dotpieces = string.split(neuromlfilename,'.')
## <filebase>_withcells.<fileextension>
outfilename = '.'.join(dotpieces[0:-1]) + '_withcells' + '.' + dotpieces[-1] 
print "writing",outfilename
#write out XML from the ElementTree instance neuromldoc
# doc.write() doesn't pretty print, hence add indents
# from moose.neuro.utils.py; copied from from http://effbot.org/zone/element-lib.htm
indent(neuromlroot, level=0)
neuromldoc.write(outfilename)
####### minidom
#neuromlfile = open(filename,'w')
#xml.dom.ext.PrettyPrint(neuromldoc,neuromlfile)
#neuromlfile.close()
