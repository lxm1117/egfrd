from _gfrd import (CuboidalRegion, SphericalSurface, CylindricalSurface, DiskSurface, PlanarSurface, Cylinder, Sphere, )
from single import (NonInteractionSingle)
from multi import (Multi)

structure_keywords = [(CuboidalRegion, 'CUBE'), (SphericalSurface, 'SPHERE'), (CylindricalSurface, 'CYLINDER'), (DiskSurface, 'DISK'), (PlanarSurface, 'PLANE')]


def format_shape(s):
   
    for object_type, key in structure_keywords:
        if isinstance(s, object_type):
            shapename = key

    if isinstance(s, CuboidalRegion):
        return '%s\tPosition=%e,%e,%e\tHalfExtent=%e,%e,%e' % (shapename, s.shape.position[0], s.shape.position[1], s.shape.position[2], s.shape.half_extent[0], s.shape.half_extent[1], s.shape.half_extent[2])
    elif isinstance(s, SphericalSurface):
        return '%s\tPosition=%e,%e,%e\tRadius=%e' % (shapename, s.shape.position[0], s.shape.position[1], s.shape.position[2], s.shape.radius)
    elif isinstance(s, CylindricalSurface):
        return '%s\tPosition=%e,%e,%e\tRadius=%e\tUnitZ=%e,%e,%e\tHalfLength=%e' % (shapename, s.shape.position[0], s.shape.position[1], s.shape.position[2], s.shape.radius, s.shape.unit_z[0], s.shape.unit_z[1], s.shape.unit_z[2], s.shape.half_length )
    elif isinstance(s, DiskSurface):
        return '%s\tPosition=%e,%e,%e\tRadius=%e\tUnitZ=%e,%e,%e' % (shapename, s.shape.position[0], s.shape.position[1], s.shape.position[2], s.shape.radius, s.shape.unit_z[0], s.shape.unit_z[1], s.shape.unit_z[2] )
    elif isinstance(s, PlanarSurface):
        return '%s\tPosition=%e,%e,%e\tHalfExtent=%e,%e\tUnitX=%e,%e,%e\tUnitY=%e,%e,%e\tOneSide=%d' % (shapename, s.shape.position[0], s.shape.position[1], s.shape.position[2], s.shape.half_extent[0], s.shape.half_extent[1], s.shape.unit_x[0], s.shape.unit_x[1], s.shape.unit_x[2],  s.shape.unit_y[0],s.shape.unit_y[1],s.shape.unit_y[2], s.shape.is_one_sided)
    return shapename

def format_shell(s):
   
    if isinstance(s, Sphere):
        return 'S'
    elif isinstance(s, Cylinder):
        return 'C\tUnitZ=%e,%e,%e\tHalfLength=%e' % (s.unit_z[0], s.unit_z[1], s.unit_z[2], s.half_length )
    return '';

def export(filename, s, a=None):
    w = s.world;
    oflag = 'w'
    if a != None and a : oflag = 'a'
    f = open(filename, oflag)
    
    f.write('sim_step=%d\tsim_time=%e\t' % (s.step_counter, s.t))
    f.write('world_size=%e\tmatrix_size=%d\n' % (w.world_size, w.matrix_size))

    for structure in s.get_structures():
        f.write('Structure=%d\tName=%s\tSid=%d\tShape=%s\n' % (structure.id.serial, structure.name, structure.sid.serial, format_shape(structure)))

    for pid, particle in w:
        f.write('Particle=%d\tRadius=%e\tPosition=%e,%e,%e\tSid=%d\n' % (pid.serial, particle.radius, particle.position[0], particle.position[1], particle.position[2], particle.sid.serial))

    for domain in s.domains.itervalues():
       f.write('Domain=%d\tMulti=%d\tCount=%d' % (domain.domain_id.serial, domain.multiplicity, domain.num_shells))
       for shell_id, shell in domain.shell_list:
          f.write('\tShell=%d\tCode=%d\tRadius=%e\tPosition=%e,%e,%e\t%s' % (shell_id.serial, 0 if isinstance(domain, NonInteractionSingle) and domain.is_reset() else 2 if isinstance(domain, Multi) else 1, shell.shape.radius, shell.shape.position[0], shell.shape.position[1], shell.shape.position[2], format_shell(shell.shape)))
       f.write('\n')
    f.write('End\n')
    f.close()
