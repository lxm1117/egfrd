#!/usr/bin/env python

import unittest
from egfrd import *
import _gfrd
import model
import myrandom

class CylindricalSurfaceTestCase(unittest.TestCase):
    '''Unittests creates cylinderical surfaces and tests the properties 
    position, projection-point, and average z-axis distance 
    of the cylinders.'''

    def setUp(self):
        '''Creates the world, particle model and structure id's.'''

        # Constants
        self.radius = 1
        self.L = 10

        # Create world with parent id
        self.world = _gfrd.World(self.L, 3)
        self.parent_id = self.world.get_def_structure_id()

        # Add cylinderical structure
        self.structure_type = _gfrd.StructureType()
        self.structure_type['name'] = 'Cylinder'

        # Connect to particle model
        self.particle_model = model.ParticleModel(self.L)
        self.particle_model.add_structure_type(self.structure_type)

        # Random radius values.
        self.r1 = 1.6
        self.r2 = 2.1
        self.r3 = 2.7
        self.r4 = 4.2


    def cylinder_at_position(self, h):
        '''Creates a cylinderical surface.
        Argurment - h: the closed point on the cylinder to the origin of the z-axis.'''
        return model.create_cylindrical_surface(self.structure_type.id, 'd', [5, 5, h], 
                                                self.radius, [0, 0, 1], self.L, self.parent_id)


    def test_random_positions(self):
        '''Tests whether the point closed to the origin is at the middle of the world
        on average, for a set of cylinders along the z-axis between 0 and L.'''
        for z in range(0, self.L):
            d = self.cylinder_at_position(z)

            positions = []
            for i in range(100):
                position = d.random_position(myrandom.rng)
                position = apply_boundary(position, self.world.world_size)
                positions.append(position)

            average_position = numpy.average(positions, 0)
            assert average_position[0] == 5
            assert average_position[1] == 5
            assert 4 < average_position[2] < 6


    def test_projected_point(self):
        '''Tests the projection-point value for a single cylinder.'''
        d = self.cylinder_at_position(self.r1)
        assert (d.project_point([5, 2, 2])[0] == [5, 5, 2]).all()
        assert d.project_point([5, 2, 2])[1][0] == 3


    def test_distance_to_cylinder(self):
        '''Tests whether the world distance function returns the correct value
        for a given a single cylinder.'''
        d = self.cylinder_at_position(self.r1)
        assert self.world.distance(d.shape, [7, 5, self.r2]) == 1
        assert self.world.distance(d.shape, [5, 5, self.r3]) == -1


# Runs the unittests
if __name__ == "__main__":
    unittest.main()

