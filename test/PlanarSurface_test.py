#!/usr/bin/env python

import unittest
from egfrd import *
import _gfrd
import myrandom
import model

class PlanarSurfaceTestCase(unittest.TestCase):
    '''Unittests creates planar surfaces and tests the properties 
    position, projection-point, and average distance to the plane.'''

    def setUp(self):
        '''Creates the world, particle model and structure id's.'''

        # Creates the world.
        self.radius = 1
        self.L = 10
        
        # Create world with parent id.
        self.world = _gfrd.World(self.L, 3)
        self.parent_id = self.world.get_def_structure_id()

        # Add planar structure
        self.structure_type = _gfrd.StructureType()
        self.structure_type['name'] = 'Planar'
        
        # Connect to particle model
        self.particle_model = model.ParticleModel(self.L)
        self.particle_model.add_structure_type(self.structure_type)

        # Some const distance values.
        self.r1 = 1.6
        self.r2 = 2.1
        self.r3 = 2.7
        self.r4 = 4.2


    def membrane_at_position(self, x, y):
        '''Creates a cylinderical surface.
        Argument - x: axis x origin of the plane.
                 - y: axis y origin of the plane.'''
        return model.create_planar_surface(self.structure_type.id, 'm', [x, y, 5],
                                           [1, 0, 0], [0, 1, 0], self.L, self.L, self.parent_id)


    def test_random_positions(self):
        '''Tests whether the point close to the origin is at the middle of the world
        on avarage, for a set of planes along the x-axis and y-axis between 0 and L.'''
        for x in range(0, self.L):
            for y in range(0, self.L):
                m = self.membrane_at_position(x, y)

                positions = []
                for i in range(100):
                    position = m.random_position(myrandom.rng)
                    position = apply_boundary(position, self.world.world_size)
                    positions.append(position)

                average_position = numpy.average(positions, 0)
                assert 4 < average_position[0] < 6
                assert 4 < average_position[1] < 6
                assert average_position[2] == 5


    def test_projected_point(self):
        '''Test the projection-point value for a single planar surface.'''
        m = self.membrane_at_position(self.r1, self.r2)
        assert (m.project_point([2, 2, 2])[0] == [2, 2, 5]).all()
        assert m.project_point([2, 2, 2])[1][0] == -3


    def test_distance_to_plane(self):
        '''Tests whether the world distance function returns the correct value
        for a given single planar surface.'''
        m = self.membrane_at_position(self.r1, self.r2)
        assert self.world.distance(m.shape, [self.r3, self.r4, 2]) == 3
        assert self.world.distance(m.shape, [self.r3, self.r4, 8]) == 3


if __name__ == "__main__":
    unittest.main()
