#from elementtree.SimpleXMLWriter import XMLWriter
## I am not using XMLWriter as it is sequential
## we want to construct a document object model dom in non-sequential order and finally write to xml
## I should use 'from xml.etree import ElementTree as ET', it is far faster than 'minidom' and non-sequential.
import xml.dom.minidom

class NeuroML_writer():
    
    def __init__(self,length_units,units):
        """
        serial order of setting sub-elements is assumed. Thus there are lots of current_<element> variables.
        However, you can simultaneously build up say instances and connections (sub-sub-sub-elements of populations and projections),
        since populations and projections elements are disjoint siblings and single in number.
        The populations, projections, inputs, instances elements are only created when the first population, projection, input, instance is set.
        This is because neuroml does not allow empty elements.
        units can be "SI Units" or "Physiological Units" - Note capitalization.
        length_units can "meter" or "micrometer".
        """
        self.doc = xml.dom.minidom.Document()
        self.neuroml = self.doc.createElement("neuroml")
        self.neuroml.setAttribute("xmlns","http://morphml.org/neuroml/schema")
        self.neuroml.setAttribute("xmlns:xsi","http://www.w3.org/2001/XMLSchema-instance")
        self.neuroml.setAttribute("xmlns:meta","http://morphml.org/metadata/schema")
        self.neuroml.setAttribute("xmlns:mml","http://morphml.org/morphml/schema")
        self.neuroml.setAttribute("xmlns:nml","http://morphml.org/networkml/schema")
        self.neuroml.setAttribute("xmlns:bio","http://morphml.org/biophysics/schema")
        self.neuroml.setAttribute("xmlns:cml","http://morphml.org/channelml/schema")
        self.neuroml.setAttribute("xsi:schemaLocation","http://morphml.org/neuroml/schema")
        self.neuroml.setAttribute("lengthUnits",length_units)
        self.doc.appendChild(self.neuroml)
        self.units = units
        self.cells = None
        self.populations = None
        self.projections = None
        self.inputs = None
        
    def note(self, parent, text):
        note_el = self.doc.createElement("meta:notes")
        parent.appendChild(note_el)
        text_node = self.doc.createTextNode(text)
        note_el.appendChild(text_node)
                
    def writeToFile(self, filename):
        neuromlfile = open(filename,'w')
        #xml.dom.ext.PrettyPrint(self.doc,neuromlfile)
        prettyxmlstr = self.doc.toprettyxml()
        neuromlfile.write(prettyxmlstr)
        neuromlfile.close()

    # methods to specify instances of cells:
    
    def set_population(self, name, cell_type):
        if self.populations == None:
            self.populations = self.doc.createElement("populations")
            self.populations.setAttribute("xmlns","http://morphml.org/networkml/schema")
            self.neuroml.appendChild(self.populations)
        self.current_population = self.doc.createElement("population")
        self.current_population.setAttribute("name",name)
        self.current_population.setAttribute("cell_type",cell_type)
        self.populations.appendChild(self.current_population)
        self.current_instances = None

    def set_instance(self, elementid,x,y,z,zrotation=None):
        """
        size is a required attribute of instances. So everytime an instance is created, the size attribute is updated.
        x,y,z is a translation. zrotation is a rotation about z axis.
        No support for zrotation in neuroml, but I put it in as meta:notes.
        """
        if self.current_instances == None:
            self.current_instances = self.doc.createElement("instances")
            self.current_instances_size = 0
            self.current_instances.setAttribute("size","0")
            self.current_population.appendChild(self.current_instances)
        instance = self.doc.createElement("instance")
        instance.setAttribute("id",elementid)
        self.current_instances.appendChild(instance)
        location = self.doc.createElement("location")
        location.setAttribute("x",x)
        location.setAttribute("y",y)
        location.setAttribute("z",z)
        instance.appendChild(location)
        if zrotation is not None:
            self.note(instance,"zrotation="+str(zrotation))
        self.current_instances_size += 1
        self.current_instances.setAttribute("size",str(self.current_instances_size))

    # methods to specify connections:
        
    def set_projection(self, name, source, target):
        if self.projections == None:
            self.projections = self.doc.createElement("projections")
            self.projections.setAttribute("units",self.units)
            self.projections.setAttribute("xmlns","http://morphml.org/networkml/schema")
            self.neuroml.appendChild(self.projections)
        self.current_projection = self.doc.createElement("projection")
        self.current_projection.setAttribute("name",name)
        self.current_projection.setAttribute("source",source)
        self.current_projection.setAttribute("target",target)
        self.projections.appendChild(self.current_projection)
        self.current_connections = None
        
    def set_synapse_props(self, synapse_type, weight, threshold, delay):
        synapse_props = self.doc.createElement("synapse_props")
        synapse_props.setAttribute("synapse_type",synapse_type)
        synapse_props.setAttribute("weight",str(weight))
        synapse_props.setAttribute("threshold",str(threshold))
        synapse_props.setAttribute("prop_delay",str(delay))
        self.current_projection.appendChild(synapse_props)

    def set_connection(self, elementid,pre_cell_id,post_cell_id,pre_segment_id='0',post_segment_id='0',synproplist=[]):
        """
        size is not a required attribute of connections unlike for instances.
        But we do set the size anyway. Everytime a connection is created, the size attribute is updated.
        """
        if self.current_connections == None:
            self.current_connections = self.doc.createElement("connections")
            self.current_connections_size = 0
            self.current_connections.setAttribute("size","0")
            self.current_projection.appendChild(self.current_connections)
        connection = self.doc.createElement("connection")
        connection.setAttribute("id",elementid)
        connection.setAttribute("pre_cell_id",pre_cell_id)
        connection.setAttribute("pre_segment_id",pre_segment_id)
        connection.setAttribute("post_cell_id",post_cell_id)
        connection.setAttribute("post_segment_id",post_segment_id)
        self.current_connections.appendChild(connection)
        for synprops in synproplist:
            properties = self.doc.createElement("properties")
            properties.setAttribute("synapse_type",synprops[0])
            properties.setAttribute("weight",str(synprops[1]))
            if len(synprops)>2:
                properties.setAttribute("internal_delay",str(synprops[2]))                
            connection.appendChild(properties)
        self.current_connections_size += 1
        self.current_connections.setAttribute("size",str(self.current_connections_size))
