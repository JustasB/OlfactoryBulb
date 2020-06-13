import bpy
from mathutils import Vector

def print_apic_reach():
    '''
    Prints the distance between the selected apical dendrite and the cell soma. This assumes only one
    cell is in the scene and the cell is in 'Section' interaction level.

    :return: Nothing, prints to console
    '''
    for ob in bpy.data.objects:
        if "soma" in ob.name:
            start_ob = ob

    end_ob = bpy.context.object
    length = Vector(end_ob.location - start_ob.location).length

    print(end_ob.name, length)