import bpy
from mathutils import Vector
import numpy as np

mc2glom = {}

for ob in bpy.data.objects:
    if 'Sphere' in ob.name:
        bpy.data.objects.remove(ob)

    elif 'priden' in ob.name:
        last_pt_loc = ob.matrix_world * ob.data.splines[0].bezier_points[-1].co
        mc2glom[ob.name.replace('.priden', '').replace('Mitral[', '').replace(']', '')] = list(last_pt_loc)

particles = np.array([list(p.location) for pi, p in
                      enumerate(bpy.data.objects['0 GL Particles'].particle_systems['ParticleSystem'].particles)])

print(len(particles))

taken = set()

for mig_glom_id in mc2glom:
    mig_glom_loc = Vector(mc2glom[mig_glom_id])

    dists = np.array([(Vector(p) - mig_glom_loc).length for p in particles])

    closest_particle_dists = np.argsort(dists)

    for i_closest in closest_particle_dists:
        if i_closest in taken:
            continue

        taken.add(i_closest)
        break

    print(mig_glom_id, i_closest, dists[i_closest])

    bpy.ops.mesh.primitive_uv_sphere_add(location=particles[i_closest])

for ob in bpy.data.objects:
    if 'Sphere' in ob.name:
        ob.scale = [50] * 3