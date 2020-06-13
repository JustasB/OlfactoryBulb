"""
This is the main class where NEURON cells are imported into Blender, where they are positioned and oriented within
bulbar layers. Transformed cell coordinates are saved for later instatiation by NEURON during simulations.
Synapse locations are identified and saved into synapse set files.
"""

import bpy, sys, os, time, re
import random
import os, sys
import bpy
import numpy as np
from collections import OrderedDict
from mathutils import Vector
import random
from math import pi, acos
import json
sys.path.append(os.getcwd())
from olfactorybulb import slices
from blenderneuron.blender.utils import fast_get, make_safe_filename
from blenderneuron.blender.views.vectorconfinerview import VectorConfinerView
from blenderneuron.blender.views.synapseformerview import SynapseFormerView
import blenderneuron.blender

'''
Sources:
Kikuta et. al. 2013 - TC soma distance from Glom center ~200 um
Witman and Greer 2007 - GC spine reach - from digitized figure 5.5 um
'''

def auto_start(scene):
    """
    A Blender startup script that starts the SliceBuilder on Blender startup
    """

    # Remove auto-execute command after starting
    bpy.app.handlers.scene_update_post.remove(auto_start)

    # Assuming starting at repo root
    sys.path.append(os.getcwd())

    # Create a slice builder class
    sbb = bpy.types.Object.SliceBuilder = SliceBuilderBlender()

    # from line_profiler import LineProfiler
    # lp = LineProfiler()
    # lp.add_function(sbb.add_mc)
    # lp.add_function(sbb.add_tc)
    # lp.add_function(sbb.add_gc)
    # lp.add_function(sbb.import_instance)
    # profiled_build = lp(sbb.build)
    # profiled_build()
    # lp.print_stats()

    sbb.build()



class SliceBuilderBlender:
    @property
    def node(self):
        """
        Returns the BlenderNEURON node that runs within Blender
        """

        return bpy.types.Object.BlenderNEURON_node

    @property
    def neuron(self):
        """
        Returns the client class that communicates with NEURON via methods defined in `olfactorybulb.slicebuilder.nrn`
        """
        return self.node.client

    @property
    def slice_dir(self):
        """
        Returns the path to the directory where virtual slice cell and synapse files will be saved
        """
        slice_dir = os.path.abspath(os.path.dirname(slices.__file__))
        return os.path.join(slice_dir, self.slice_name)

    def __init__(self,
                 odors=['Apple'],
                 slice_object_name='DorsalColumnSlice',
                 max_mcs=10, max_tcs=None, max_gcs=300,  # Uses mouse ratios if None
                 mc_particles_object_name='2 ML Particles',
                 tc_particles_object_name='1 OPL Particles',
                 gc_particles_object_name='4 GRL Particles',
                 glom_particles_object_name='0 GL Particles',
                 glom_layer_object_name='0 GL',
                 outer_opl_object_name='1 OPL-Outer',
                 inner_opl_object_name='1 OPL-Inner'):
        """
        Prepares the slice builder

        :param odors: A list of odors whose glomeruli to include (e.g. ['Apple', 'Mint']), use 'all' for all gloms.
        :param slice_object_name: The name of the Blender mesh that defines the shape of the virtual slice
        :param max_mcs: The maximum number of MCs to include in the model
        :param max_tcs: Maximum number of TCs. Use None to use mouse MC-TC ratio.
        :param max_gcs: Maximum number of GCs. Use None to use mouse MC-GC ratio.
        :param mc_particles_object_name: The name of the Blender object that defines MC soma locations
        :param tc_particles_object_name: The name of the Blender object that defines TC soma locations
        :param gc_particles_object_name: The name of the Blender object that defines GC soma locations
        :param glom_particles_object_name: The name of the Blender object that defines glomerulus locations
        :param glom_layer_object_name: The name of the Blender object that defines the geometry of glomerular layer
        :param outer_opl_object_name: The name of the Blender object that defines the outer boundary of the OPL layer
        :param inner_opl_object_name: The name of the Blender object that defines the inner boundary of the OPL layer
        """

        # In mouse, for each MC, there are:
        # See: model-data.sql > measurement table > mc/gc/tc_count entries for sources
        #  2.36 TCs
        if max_tcs is None:
            max_tcs = int(round(max_mcs * 2.36))

        #  16.97 GCs
        if max_gcs is None:
            max_gcs = int(round(max_mcs * 16.97))


        self.odors = odors
        self.glom_cells = {}

        self.slice_name = slice_object_name

        self.max_mcs = max_mcs
        self.max_tcs = max_tcs
        self.max_gcs = max_gcs

        self.glom_particles_name = glom_particles_object_name
        self.tc_particles_name = tc_particles_object_name
        self.mc_particles_name = mc_particles_object_name
        self.gc_particles_name = gc_particles_object_name

        self.glom_layer_object_name = glom_layer_object_name
        self.outer_opl_object_name = outer_opl_object_name
        self.inner_opl_object_name = inner_opl_object_name

        # Within slice
        self.get_cell_locations()

        self.get_cell_base_model_info()

        # Show as section objects
        self.create_groups()

        # Clear slice files
        self.clear_slice_files()

    def build(self, seed=0):
        """
        Positions MC/TC/GC models within the OB layers, identifies synapse locations, and saves the model
        for later simulation in NEURON

        :param seed: The random seed to use (to assist reproducibility)
        """

        random.seed(seed)

        self.max_alignment_angle = 35

        for i, loc in enumerate(self.mc_locs):
            print('Adding MC %s' % i)
            self.add_mc(loc)

        for i, loc in enumerate(self.tc_locs):
            print('Adding TC %s' % i)
            self.add_tc(loc)

        for i, loc in enumerate(self.gc_locs):
            print('Adding GC %s' % i)
            self.add_gc(loc)

        # Add synapse sets
        self.add_synapse_sets()

        # Select all cells in groups
        self.node.groups['MCs'].select_roots('All','mc*')
        self.node.groups['TCs'].select_roots('All','tc*')
        self.node.groups['GCs'].select_roots('All','gc*')

        gcs_with_syns = set()

        # Find and save syns in all synapse sets
        for syn_set in self.node.synapse_sets:
            file = os.path.join(self.slice_dir, make_safe_filename(syn_set.name)+'.json')
            print('Saving synapse set "'+syn_set.name+'" saved to: ' + file)

            pairs = syn_set.get_synapse_locations()

            # Note which GCs have synapses
            for pair in pairs:
                source_cell = pair.source.section_name[:pair.source.section_name.find(']')+1]
                gcs_with_syns.add(source_cell)

            syn_set.save_synapses(file)

        # Remove unconnected GCs
        # they won't contribute to simulation output
        # removing them makes the simulation smaller
        self.node.groups['GCs'].include_roots_by_name(
            [cell + '.soma' for cell in gcs_with_syns],
            exclude_others=True
        )

        # Save all cells
        for group in self.node.groups.values():
            file = os.path.join(self.slice_dir, make_safe_filename(group.name)+'.json')
            print('Saving cell group %s %s to: %s'%(len(group.roots.keys()), group.name, file))
            group.to_file(file)

        # Save glom-cell associations
        file = os.path.join(self.slice_dir, 'glom_cells.json')
        with open(file, 'w') as f:
            print('Saving glomerulus-cells links to: ' + file)
            json.dump(self.glom_cells, f)

        # Initially, reduce group display detail levels - All can be changed in Blender GUI
        self.node.groups['MCs'].interaction_granularity = 'Cell' # Clicking any cell dendrite will select the whole cell
        self.node.groups['TCs'].interaction_granularity = 'Cell'
        self.node.groups['GCs'].interaction_granularity = 'Cell'
        self.node.groups['MCs'].as_lines = True # Dendrites will be shown as 0-width lines
        self.node.groups['TCs'].as_lines = True
        self.node.groups['GCs'].as_lines = True

        # Show all group cells
        print('Creating blender scene...')
        bpy.ops.blenderneuron.display_groups()
        print('DONE')

    def add_synapse_sets(self):
        """
        Creates synapse sets between MCs-GCs and TCs-GCs
        """

        # Delete the default set
        self.node.synapse_sets.remove(0)

        self.create_synapse_set('GCs', 'MCs')
        self.create_synapse_set('GCs', 'TCs')

    def create_synapse_set(self, group_from, group_to):
        """
        Defines a synapse set between two BlenderNEURON groups of cells

        :param group_from: One of 'MCs', 'TCs', 'GCs'
        :param group_to: One of 'MCs', 'TCs', 'GCs'
        """

        new_set = self.node.add_synapse_set(group_from + '->' + group_to)
        new_set.group_source = group_from
        new_set.group_dest = group_to

        new_set.max_distance = 5
        new_set.use_radius = True
        new_set.max_syns_per_pt = 1
        new_set.section_pattern_source = "*apic*"
        new_set.section_pattern_dest = "*dend*"
        new_set.synapse_name_dest = 'GabaSyn'
        new_set.synapse_params_dest = str({
            'gmax': 0.005,  # uS,
            'tau1': 1,  # ms
            'tau2': 100  # ms
        })
        new_set.is_reciprocal = True
        new_set.synapse_name_source = 'AmpaNmdaSyn'
        new_set.synapse_params_source = str({'gmax': 0.1})
        new_set.create_spines = False
        new_set.spine_neck_diameter = 0.2
        new_set.spine_head_diameter = 1
        new_set.spine_name_prefix = 'Spine'
        new_set.conduction_velocity = 1
        new_set.initial_weight = 1
        new_set.threshold = 0

    def clear_slice_files(self):
        """
        Deletes previously saved virtual slice .json files from the slice directory
        (e.g. olfactorybulb/slices/DorsalColumnSlice/.json)
        """

        dir = self.slice_dir

        # Match e.g. 'MCs.json'
        pattern = re.compile('.+json')

        if os.path.exists(dir):
            for file in os.listdir(dir):
                if pattern.match(file) is not None:
                    os.remove(os.path.abspath(os.path.join(dir, file)))

    def get_cell_locations(self):
        """
        Identifies the locations of glomeruli, mcs, tcs, and gcs that are contained by the virtual slice.
        """

        self.globalize_slice()

        odor_glom_ids = self.neuron.get_odor_gloms(self.odors)
        self.glom_locs = self.get_locs_within_slice(self.glom_particles_name, self.slice_name, odor_glom_ids)

        self.inner_opl_locs = self.get_opl_locs(self.inner_opl_object_name, self.slice_name)
        self.outer_opl_locs = self.get_opl_locs(self.outer_opl_object_name, self.slice_name)

        self.mc_locs = self.get_locs_within_slice(self.mc_particles_name, self.slice_name, limit=self.max_mcs)
        self.tc_locs = self.get_locs_within_slice(self.tc_particles_name, self.slice_name, limit=self.max_tcs)
        self.gc_locs = self.get_locs_within_slice(self.gc_particles_name, self.slice_name, limit=self.max_gcs)

        print('Gloms:', len(self.glom_locs))
        print('TCs:', len(self.tc_locs))
        print('MCs:', len(self.mc_locs))
        print('GCs:', len(self.gc_locs))

    def get_opl_locs(self, opl_name, slice_obj_name):
        """
        Gets the coordinates of inner or outer boundaries of the OPL layer that are contained by the virtual slice
        """

        obj = bpy.data.objects[opl_name]
        wm = obj.matrix_world

        locs = fast_get(obj.data.vertices, 'co', 3)

        def globalize(loc):
            return np.array(wm * Vector(loc))

        # Globalize the coordinates
        locs = np.array(list(map(globalize, locs)))

        slice_obj = bpy.data.objects[slice_obj_name]

        return np.array([pt
                         for pt in locs
                         if self.is_inside(Vector(pt), slice_obj)])

    def get_cell_base_model_info(self):
        """
        Gets metadata info about each of the base (untransformed) MC, TC, and GC cell models
        """

        self.mc_base_models, self.tc_base_models, self.gc_base_models = \
            self.neuron.get_base_model_info()

        self.mc_base_models = OrderedDict(sorted(self.mc_base_models.items(), key=lambda i: i[0]))
        self.tc_base_models = OrderedDict(sorted(self.tc_base_models.items(), key=lambda i: i[0]))
        self.gc_base_models = OrderedDict(sorted(self.gc_base_models.items(), key=lambda i: i[0]))

        self.max_apic_mc_info = self.get_longest_apic_model(self.mc_base_models)
        self.mc_apic_lengths = self.get_apic_lengths(self.mc_base_models)

        self.max_apic_tc_info = self.get_longest_apic_model(self.tc_base_models)
        self.tc_apic_lengths = self.get_apic_lengths(self.tc_base_models)

        self.max_apic_gc_info = self.get_longest_apic_model(self.gc_base_models)
        self.gc_apic_lengths = self.get_apic_lengths(self.gc_base_models)

    @staticmethod
    def get_apic_lengths(base_models):
        """
        Gets the apical dendrite lengths of the specified base models

        :param base_models: List with metadata of base cell models
        :return: A numpy array of apical dendrite lengths
        """

        return np.array([c["apical_dendrite_reach"] for c in base_models.values()])

    def create_groups(self):
        """
        Creates empty BlenderNEURON cell groups for MCs, TCs, and GCs.
        """

        # Remove the default group
        self.node.groups['Group.000'].remove()

        # Create empty cell groups
        groups = [self.node.add_group(name, False) for name in ['MCs', 'TCs', 'GCs']]

        # show each section as blender objects - necessary for dend alignment
        for group in groups:
            group.interaction_granularity = 'Section'
            group.recording_granularity = 'Cell'
            group.record_activity = False

        # Add some color
        groups[0].default_color = [0.15, 0.71, 0.96]       # MCs - neon blue
        groups[1].default_color = [1,    0.73, 0.82]       # TCs - pink
        groups[2].default_color = [1,    0.80, 0.11]       # GCs - gold

    def globalize_slice(self):
        """
        Converts all points of the slice object to global coordinates (relative to scene origin)
        """

        # Apply all/any transformations to the slice
        slice = bpy.data.objects[self.slice_name]
        slice.select = True
        bpy.context.scene.objects.active = slice
        bpy.ops.object.transform_apply(location=True, scale=True, rotation=True)
        slice.select = False

    def add_mc(self, mc_pt):
        """
        Instantiates in NEURON, places, and orients a mitral cell within the mitral cell layer
        and confines its lateral dendrites to the curvature of the surrounding internal part
        of the outer plexiform layer.

        :param mc_pt: A dict whose 'loc' key contains a numpy array of xyz coordinates of the cell soma
        """

        # find the closest glom layer loc - cell will be pointed towards it
        closest_glom_loc, dist_to_gl = \
            self.closest_point_on_object(mc_pt['loc'], bpy.data.objects[self.glom_layer_object_name])

        longest_apic_reach = self.max_apic_mc_info["apical_dendrite_reach"]

        # Apics are too short use longest apic MC
        if dist_to_gl > longest_apic_reach:
            mc = self.max_apic_mc_info
        else:
            # get mcs with apics longer than dist to GL
            longer_idxs = np.where(self.mc_apic_lengths > dist_to_gl)[0]

            # pick a random mc from this list
            mc = self.get_random_model(self.mc_base_models, longer_idxs)

        # find a glom whose distance is as close to the length of the mc apic
        matching_glom_loc, matching_glom_id = self.find_matching_glom(mc_pt['loc'], mc)

        base_class = mc["class_name"]
        apic_glom_loc = matching_glom_loc

        # Create the selected MC in NRN
        instance_name = self.neuron.create_cell('MC', base_class)

        # Associate the cell with the glomerulus
        self.link_cell_to_glom(instance_name, matching_glom_id)

        # Import cell into Blender
        self.import_instance(instance_name, 'MCs')

        mc_soma, mc_apic_start, mc_apic_end = \
            self.get_key_mctc_section_objects(self.mc_base_models, base_class, instance_name)

        # Align apical towards the closest glom
        self.position_orient_align_mctc(mc_soma,
                                        mc_apic_start,
                                        mc_apic_end,
                                        mc_pt['loc'],
                                        closest_glom_loc,
                                        apic_glom_loc)

        # Retain the reoriented cell
        bpy.ops.blenderneuron.update_groups_with_view_data()

        self.confine_dends(
            'MCs',
            self.inner_opl_object_name,
            self.outer_opl_object_name,
            max_angle=self.max_alignment_angle,
            height_start=0,
            height_end=0.6
        )

        bpy.ops.blenderneuron.update_groups_with_view_data()

    def link_cell_to_glom(self, instance_name, matching_glom_id):
        """
        Adds a cell instance to a list of cells that belong to the specified glomerulus

        :param instance_name: The name of the cell as returned by NEURON
        :param matching_glom_id: The id of the glomerulus, with which to associate the cell
        """

        glom_cells = self.glom_cells.get(matching_glom_id, [])

        glom_cells.append(instance_name.replace('.soma',''))

        self.glom_cells[matching_glom_id] = glom_cells

    def import_instance(self, instance_name, group_name):
        """
        Imports and shows a cell instantiated in NEURON into Blender using a BlenderNEURON group

        :param instance_name: The name of the cell model to import (as named by NEURON)
        :param group_name: The name of the group to associate the cell with
        """

        # Get updated list of NRN cells in Blender
        bpy.ops.blenderneuron.get_cell_list_from_neuron()

        # Select the created instance
        group = self.node.groups[group_name]
        group.include_roots_by_name([instance_name], exclude_others=True)

        # Import group with the created cell and show it
        group.import_group()
        group.show()

    def add_tc(self, tc_pt):
        """
        Similar to `add_mc(pt)`. Instantiates, places, and orients a tufted cell model and
        confines its lateral dendrites to the outer portion of the outer plexiform layer.

        :param tc_pt: A dict with 'loc' that contains a numpy array of xyz coordinates of the TC soma
        """

        # find the closest glom layer loc - cell will be pointed towards it
        closest_glom_loc, dist_to_gl = \
            self.closest_point_on_object(tc_pt['loc'], bpy.data.objects[self.glom_layer_object_name])

        longest_apic_reach = self.max_apic_tc_info["apical_dendrite_reach"]

        # Apics are too short use longest apic TC
        if dist_to_gl > longest_apic_reach:
            tc = self.max_apic_tc_info
        else:
            # Apics are longer than distance
            # get tcs with apics longer than the closest glom,
            # but no further than ~200 um from glom (Source: Kikuta et. al. 2013)
            longer_idxs = np.where((self.tc_apic_lengths > dist_to_gl) &
                                   (self.tc_apic_lengths - dist_to_gl < 200))[0]

            # pick a random tc from this list
            tc = self.get_random_model(self.tc_base_models, longer_idxs)

        # find a glom whose distance is as close to the length of the tc apic
        matching_glom_loc, matching_glom_id = self.find_matching_glom(tc_pt['loc'], tc)

        base_class = tc["class_name"]
        apic_glom_loc = matching_glom_loc

        # Create the selected TC in NRN
        instance_name = self.neuron.create_cell('TC', base_class)

        # Associate the cell with the glomerulus
        self.link_cell_to_glom(instance_name, matching_glom_id)

        # Import it into Blender
        self.import_instance(instance_name, 'TCs')

        soma, apic_start, apic_end = \
            self.get_key_mctc_section_objects(self.tc_base_models, base_class, instance_name)

        # Align apical towards the closest glom
        self.position_orient_align_mctc(soma,
                                        apic_start,
                                        apic_end,
                                        tc_pt['loc'],
                                        closest_glom_loc,
                                        apic_glom_loc)

        # Retain the reoriented cell
        bpy.ops.blenderneuron.update_groups_with_view_data()

        self.confine_dends(
            'TCs',
            self.inner_opl_object_name,
            self.outer_opl_object_name,
            max_angle=self.max_alignment_angle,
            height_start=0.4,
            height_end=1.0
        )

        bpy.ops.blenderneuron.update_groups_with_view_data()

    def find_closest_glom(self, cell_loc):
        """
        Finds the closest glomerulus to the specified coordinates

        :param cell_loc: A numpy array of xyz coordinates
        :return: A tuple with the location of the closest glomerulus and distance to it
        """

        # Get distances to individual gloms
        glom_dists = self.dist_to_gloms(cell_loc)

        matching_glom_idx = np.argmin(glom_dists)
        matching_glom = self.glom_locs[matching_glom_idx]

        return matching_glom['loc'], glom_dists[matching_glom_idx]

    def find_matching_glom(self, cell_loc, cell_model_info):
        """
        Finds a glomerulus that is approximately the same distance from the soma
        as the length of the cells apical dendrite

        :param cell_loc: A numpy xyz array with coordinates
        :param cell_model_info: Cell model metadata dict with 'apical_dendrite_reach' key that contains the apical dendrite length
        :return: The xyz location and the id of the matching glomerulus
        """

        # Get distances to individual gloms
        glom_dists = self.dist_to_gloms(cell_loc)

        matching_glom_idx = np.argmin(np.abs(glom_dists - cell_model_info["apical_dendrite_reach"]))
        matching_glom = self.glom_locs[matching_glom_idx]

        return matching_glom['loc'], matching_glom['id']

    def get_opl_distance_info(self, cell_loc, pts):
        """
        Computes distances to and the closest point of an outer plexiform layer to a given point

        :param cell_loc: The xyz location of a point
        :param pts: The list of xyz coordinates of the OPL layer mesh
        :return: A tuple with closest location on the layer, the distance to that point, and a list of distances to all OPL points
        """

        dists = self.dist_to(pts, cell_loc)
        closest_idxs = np.argsort(dists)

        closest_loc = pts[closest_idxs][0]
        closest_dist = dists[closest_idxs][0]

        return closest_loc, closest_dist, dists

    @staticmethod
    def closest_point_on_object(global_pt, mesh_obj):
        """
        Gets the closest point and distance of a mesh object from a specified point

        :param global_pt: The point from which to measure distance
        :param mesh_obj: A mesh object
        :return: A tuple with an xyz numpy array of the closest point and distance to it
        """

        local_pt = mesh_obj.matrix_world.inverted() * Vector(global_pt)

        _, mesh_pt, _, _ = mesh_obj.closest_point_on_mesh(local_pt)

        mesh_pt_global = mesh_obj.matrix_world * mesh_pt

        dist = Vector(global_pt - mesh_pt_global).length

        return np.array(mesh_pt_global), dist

    def add_gc(self, gc_pt):
        """
        Instantiates, places, and orients a granule cell in the granule cell layer

        :param gc_pt: XYZ array of GC soma location
        """

        # find the closest glom layer loc - cell will be pointed towards it
        glom_loc, glom_dist = \
            self.closest_point_on_object(gc_pt['loc'], bpy.data.objects[self.glom_layer_object_name])
        # glom_loc, glom_dist = self.find_closest_glom(gc_pt['loc'])

        # find the closest inner opl loc - apics must go beyond this
        closest_iopl_loc, dist_to_iopl = \
            self.closest_point_on_object(gc_pt['loc'], bpy.data.objects[self.inner_opl_object_name])

        # find the closest outer opl loc - apics must stay under this
        closest_oopl_loc, dist_to_oopl = \
            self.closest_point_on_object(gc_pt['loc'], bpy.data.objects[self.outer_opl_object_name])

        # Find such GC models whose apics are confined to the OPL
        # Specifically:
        # Get gcs with apics longer than the closest opl
        # AND
        # the apic does not exceed outer opl (an external margin is allowed)
        min_length = dist_to_iopl
        max_length = dist_to_oopl + 30

        matching_idxs = np.where((self.gc_apic_lengths > min_length) &
                                 (self.gc_apic_lengths < max_length))[0]

        # If no cell can be confined to OPL, leave that location blank
        if len(matching_idxs) == 0:
            return

        # pick a random gc from this list
        gc = self.get_random_model(self.gc_base_models, matching_idxs)

        base_class = gc["class_name"]
        apic_target_loc = glom_loc

        # Create the selected GC in NRN
        instance_name = self.neuron.create_cell('GC', base_class)

        # Import it into Blender
        self.import_instance(instance_name, 'GCs')

        soma, apic_start, apic_end = \
            self.get_key_mctc_section_objects(self.gc_base_models, base_class, instance_name)

        self.position_orient_cell(soma, apic_end, gc_pt['loc'], apic_target_loc)

        # Retain the reoriented cell
        bpy.ops.blenderneuron.update_groups_with_view_data()

    def get_random_model(self, base_models, longer_idxs):
        """
        Given a list of indices stored in longer_idxs, selects a random element of base_models

        :param base_models: A list of base cell models
        :param longer_idxs: A list of indices from which to pick a random index
        :return: A random cell model
        """

        rand_idx = longer_idxs[random.randrange(len(longer_idxs))]
        cell = list(base_models.values())[rand_idx]
        return cell

    def confine_dends(self, group_name, start_layer_name, end_layer_name, max_angle, height_start, height_end):
        """
        Confines the lateral dendrits of cells in a group to follow the curvature between two layers

         Height_start and height_end specify fractions that define the corridor between the layers to which
         the dendrites should be confined. E.g. To confine the dendrites in the halfway between the two layers,
         closer to the start layer, set height_start and height_end to 0 and 0.5. Set them to 0.5 and 1.0 to confine
         to the halfway that is closer to the end region. 1.0 corresponds to the local distance between the two layers.

         The two layers should be 'locally' parallel. E.g. two planes, or two concentric spheres. In OB model case,
         the OB layers are complex shapes, but concentric-like.

        :param group_name: The name of the group of cells to confine
        :param start_layer_name: The name of the first confinement layer
        :param end_layer_name: The name of the second confinement layer
        :param max_angle: The maximum angle that a dendritic branch can rotate to be confined between layers
        :param height_start: Fraction between 0-1
        :param height_end: Fraction between 0-1
        """

        group = self.node.groups[group_name]

        # Set the layers
        group.set_confiner_layers(start_layer_name, end_layer_name, max_angle, height_start, height_end)
        group.setup_confiner()
        group.confine_between_layers()

    def save_transform(self, group_name, instance_name):
        """
        Saves trnsformed cells in a cell group to a NEURON compatible file

        :param group_name: The name of the BlenderNEURON cell group to save
        :param instance_name: The name of the file to save. Can be the name of the cell soma segment.
        """


        group = self.node.groups[group_name]

        # Make instance name a valid python module name (eg: MC1[0].soma -> MC1_0.py)
        file_name = instance_name \
                        .replace("].soma", "") \
                        .replace("[", "_") + \
                        ".py"

        # Save to slice folder
        path = os.path.join(self.slice_dir, file_name)

        # Save cells part of the group as files
        group.to_file(path)

    def get_key_mctc_section_objects(self, base_model_dict, base_class, instance_name):
        """
        Gets the blender objects of soma, apical dendrite start (base), and apical dendrite end (tuft) sections

        :param base_model_dict: A dicttionary of base cell model metadata info
        :param base_class: The name of a base cell class
        :param instance_name: The name of blender object of cell soma section
        :return: Blender objects of soma, apical dendrite start, and apical dendrite end
        """

        cell_info = base_model_dict[base_class]
        apic_pattern = instance_name.replace(".soma", "") + '.apic[%s]'

        bpy_objects = bpy.data.objects
        soma = bpy_objects[instance_name]

        apic_start = bpy_objects.get(apic_pattern % cell_info["apical_dendrite_start"])
        apic_end = bpy_objects.get(apic_pattern % cell_info["apical_dendrite_end"])

        return soma, apic_start, apic_end

    def get_longest_apic_model(self, base_model_dict):
        """
        Returns the metadata of the base cell model that has that longest apical dendrite

        :param base_model_dict: The dict of base cell model metadatas
        :return: A metadata object
        """

        cell_names, apic_lengths = zip(*[(cell["class_name"], cell["apical_dendrite_reach"])
                                         for cell in base_model_dict.values()])

        max_apic_idx = np.argmax(apic_lengths)

        return base_model_dict[cell_names[max_apic_idx]]

    def dist_to_gloms(self, loc):
        """
        Returns an array of distances to glomeruli from a given location

        :param loc: XYZ coordinate
        :return: Distances to glomeruli
        """

        return self.dist_to(np.array([glom['loc'] for glom in self.glom_locs]), loc)

    @staticmethod
    def dist_to(targets_array, loc):
        """
        Returns an array of distances to the specified list of points from a given point

        :param targets_array: The array of xzy target coordinates
        :param loc: The target coordinate
        :return: An array of distances
        """

        return np.sqrt(np.sum(np.square(targets_array - loc), axis=1))

    def get_locs_within_slice(self, particle_obj_name, slice_obj_name, allowed_particles=None, limit=None):
        """
        Gets a list of particle locations that are contained by a slice object

        :param particle_obj_name: A blender particle object
        :param slice_obj_name: A blender mesh object that represents the virtual slice
        :param allowed_particles: A list of particle ids that are allowed to be included. If None, not restricted.
        :param limit: A list of dicts with ids and locs of matching points
        :return:
        """

        particles_obj = bpy.data.objects[particle_obj_name]
        particles = particles_obj.particle_systems[0].particles
        particles_wm = particles_obj.matrix_world
        slice_obj = bpy.data.objects[slice_obj_name]

        result = [{ 'id': pid, 'loc': np.array(particles_wm * ptc.location)}
                           for pid, ptc in enumerate(particles)
                           if (allowed_particles is None or pid in allowed_particles) and
                           self.is_inside(particles_wm * ptc.location, slice_obj)]

        if limit is not None:
            print('Selecting %s/%s %s locations inside slice'%(limit, len(result), particle_obj_name))
            result = result[:limit]

        return result

    @staticmethod
    def is_inside(target_pt_global, mesh_obj, tolerance=1):
        """
        Determines if a target point is inside a mesh

        :param target_pt_global: Target xyz point, in global coordinates
        :param mesh_obj: Target mesh object
        :param tolerance: A tolerance in degrees to account for rounding error in detecting points inside. <=1 is generally sufficient.
        :return: True or False
        """

        # Convert the point from global space to mesh local space
        target_pt_local = mesh_obj.matrix_world.inverted() * target_pt_global

        # Find the nearest point on the mesh and the nearest face normal
        _, pt_closest, face_normal, _ = mesh_obj.closest_point_on_mesh(target_pt_local)

        # Get the target-closest pt vector
        target_closest_pt_vec = (pt_closest - target_pt_local).normalized()

        # Compute the dot product = |a||b|*cos(angle)
        dot_prod = target_closest_pt_vec.dot(face_normal)

        # Get the angle between the normal and the target-closest-pt vector (from the dot prod)
        angle = acos(min(max(dot_prod, -1), 1)) * 180 / pi

        # Allow for some rounding error
        inside = angle < 90 - tolerance

        return inside

    def unparent(self, obj):
        """
        Removes the parent of a Blender object, and keeps the object in the same location

        :param obj: The blender object to unparent
        :return: The parent of the target object
        """

        prev_parent = obj.parent
        parented_wm = obj.matrix_world.copy()
        obj.parent = None
        obj.matrix_world = parented_wm
        return prev_parent

    def parent(self, obj, parent):
        """
        Make one blender object a parent of another

        :param obj: The child object
        :param parent: The child's parent object
        """

        obj.parent = parent
        obj.matrix_parent_inverse = parent.matrix_world.inverted()

    def position_orient_align_mctc(self, soma, apic_start, apic_end, loc, closest_glom_loc, apic_glom_loc):
        """
        TODO: this appears to have the same/similar function as self.position_orient_cell. This may be
        redundant code.

        :param soma:
        :param apic_start:
        :param apic_end:
        :param loc:
        :param closest_glom_loc:
        :param apic_glom_loc:
        :return:
        """

        # Position and 'point' cell towards closest glom
        self.position_orient_cell(soma, apic_end, loc, closest_glom_loc)

        # Temporarily unparent the apic start (location becomes global)
        apic_start_parent = self.unparent(apic_start)
        apic_end_parent = self.unparent(apic_end)

        # Compute the start and end alignment vectors (start apic->end apic TO start apic->glom)
        # Relative to apic_start
        apic_start_wmi = apic_start.matrix_world.inverted()
        apic_end_loc = apic_start_wmi * apic_end.location

        startVec = Vector(apic_end_loc)
        endVec = Vector(apic_start_wmi * Vector(apic_glom_loc))

        # Reparent the apic end (so it rotates with the start apic)
        self.parent(apic_end, apic_end_parent)

        # Compute rotation quaternion and rotate the start apic by it
        initMW = apic_start.matrix_world.copy()
        rotM = startVec.rotation_difference(endVec).to_matrix().to_4x4()
        apic_start.matrix_world = initMW * rotM

        # Reparent the apic end (so it rotates with the start apic)
        self.parent(apic_start, apic_start_parent)

        bpy.context.scene.update()

    def position_orient_cell(self, soma, apic_end, soma_loc, closest_glom_loc):
        """
        Position cell at a location, rotate it around its apical dendrite axis,
        and rotate it around the soma so it's apical dendrite 'points' towards
        the closest glomerulus location.

        :param soma: Blender object that holds the soma section
        :param apic_end: Blender object that holds the apical dendrite (furthest apical section on/near apical axis)
        :param soma_loc: The desired cell soma location
        :param closest_glom_loc: The location of the closest glomerulus
        """

        # Add random rotation around the apical axis
        soma.rotation_euler[2] = random.randrange(360) / 180.0 * pi

        # Position the soma
        soma.location = soma_loc

        # Update child matrices
        bpy.context.scene.update()

        # Align the soma to be orthogonal to the soma-closest glom vector
        soma_wmi = soma.matrix_world.inverted()
        apic_end_loc = Vector(soma_wmi * apic_end.matrix_world * apic_end.location)
        glom_loc = Vector(soma_wmi * Vector(closest_glom_loc))

        initMW = soma.matrix_world.copy()
        rotM = apic_end_loc.rotation_difference(glom_loc).to_matrix().to_4x4()
        soma.matrix_world = initMW * rotM

        # Update child matrices
        bpy.context.scene.update()

    def extend_apic(self, apic_start, apic_end, apic_glom_loc):
        """
        TODO: This is probably unused, leftover from an earlier version.

        :param apic_start:
        :param apic_end:
        :param apic_glom_loc:
        :return:
        """

        # Relative to apic_start
        glom_loc = apic_start.matrix_world.inverted() * Vector(apic_glom_loc)
        apic_end_loc = apic_start.matrix_world.inverted() * apic_end.matrix_world * apic_end.location

        apic_glom_diff = Vector(glom_loc - apic_end_loc)

        apic_start.location = apic_start.location.copy() + apic_glom_diff



# This makes it so the slice builder automatically runs when Blender loads
bpy.app.handlers.scene_update_post.append(auto_start)
