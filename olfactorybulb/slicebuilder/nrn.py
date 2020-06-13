from prev_ob_models.Birgiolas2020.isolated_cells import *
from olfactorybulb.database import CellModel, Odor, OdorGlom


class SliceBuilderNRN:
    """
    Establishes an interface for Blender to load cells into NEURON
    """

    def __init__(self):
        """
        Starts NEURON, creates a BlenderNEURON node, and adds methods
        to the node that Blender can call.
        """

        # Load NRN
        from neuron import h, gui

        # Connect to Blender addon
        from blenderneuron.neuronstart import BlenderNEURON as addon

        # Add methods that can be called from Blender
        addon.server.register_function(self.get_base_model_info)
        addon.server.register_function(self.create_cell)
        addon.server.register_function(self.get_odor_gloms)

        # Keep track of cells
        self.cells = {
            "MC": [],
            "TC": [],
            "GC": [],
        }

    def get_base_model_info(self):
        """
        Returns metadata information from the model database about each base cell model.

        A base cell model is a NEURON model of a cell with a representative morphology
        centered at the origin and rotated so the apical dendrites are aligned with the
        Z axis.

        The based models are positioned within layers, and oriented appropriately. MC and TC
        lateral dendrites are also modified so they follow the curvature of the external
        plexiform layer.

        :return: A tuple of mc, tc, and gc dictionaries of metadata
        """

        mc_base_models = self.get_model_info('MC')
        tc_base_models = self.get_model_info('TC')
        gc_base_models = self.get_model_info('GC')

        return mc_base_models, tc_base_models, gc_base_models

    def get_model_info(self, cell_type='MC'):
        """
        Returns the list of metadata of base cell models of a given cell type.

        See olfactorybulb.database.CellModel for metadata field names.

        :param cell_type: One of 'MC', 'TC', or 'GC'
        :return: A dictionary of base cell model class names and their metadata
        """

        return {
            cm["class_name"]: cm for cm in CellModel \
            .select() \
            .where((CellModel.source_id == 'Birgiolas (2020)') & (CellModel.cell_type == cell_type)) \
            .order_by(CellModel.class_name)
            .dicts() \
            }

    def create_cell(self, type, class_name):
        """
        Creates an instance of a base cell model in NEURON.

        :param type: One of 'MC', 'TC', or 'GC'
        :param class_name: One of the cell classes in prev_ob_models.Birgiolas2020.isolated_cells
        :return: The name that NEURON gives to the cell soma segment e.g. 'MC5[0].soma'
        """

        exec ("cell = " + class_name + "()")

        self.cells[type].append(cell)

        return str(cell.soma.name())

    def get_odor_gloms(self, odors):
        """
        Gets the list of glomerular ids that are activated by the specified list of odors

        :param odors: A list of odor names, or 'all'
        :return: List of glomerulus ids, or None when 'all' odors were specified
        """

        if odors == 'all':
            return None

        result = set()

        for odor in odors:
            odor_gloms = OdorGlom \
                .select(OdorGlom.glom_id) \
                .distinct() \
                .join(Odor) \
                .where(Odor.name == odor)

            for glom in odor_gloms:
                result.add(glom.glom_id)

        return list(result)
