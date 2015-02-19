#!/usr/bin/env python

import unittest
import model
import _gfrd
import gfrdbase
import myrandom


class BDSimulatorTestCase(unittest.TestCase):
    '''Unittests creates Brownian Dynamics simulator instance and test
    the properties instance creation, particle mobility and simulation time.'''

    def setUp(self):
        '''Creates the world with particle species and simulator.'''
        self.m = model.ParticleModel(1e-5)

        self.A = model.Species('A', 0, 1e-8)
        self.B = model.Species('B', 2e-11, 5e-9)
        self.C = model.Species('C', 2e-11, 5e-8)
        self.m.add_species_type(self.A)
        self.m.add_species_type(self.B)
        self.m.add_species_type(self.C)

        self.m.set_all_repulsive()
        self.w = gfrdbase.create_world(self.m, 10)
        self.nrw = _gfrd.NetworkRulesWrapper(self.m.network_rules)

        self.s = _gfrd._BDSimulator(self.w, self.nrw, myrandom.rng)
        self.s.set_reaction_length_factor(0.05, 0.01)

    
    def test_instantiation(self):
        '''Tests whether an BDSimulator instance can be created.'''
        self.failIf(self.s == None)

    
    def test_one_particle(self):
        '''Tests whether the simulation time increases for a single particle simulation.'''
        gfrdbase.place_particle(self.w, self.C, [0.0,0.0,0.0])

        t = self.s.t
        for i in range(5):
            self.s.step()
        self.failIf(t == self.s.t)


    def test_two_particles(self):
        '''Tests whether the simulation time increases for a two particle simulation.'''
        gfrdbase.place_particle(self.w, self.C, [0.0,0.0,0.0])
        gfrdbase.place_particle(self.w, self.C, [5e-6,5e-6,5e-6])

        t = self.s.t
        for i in range(5):
            self.s.step()
        self.failIf(t == self.s.t)


    def test_three_particles(self):
        '''Tests whether the simulation time increases for a three particle simulation.'''
        gfrdbase.place_particle(self.w, self.C, [0.0,0.0,0.0])
        gfrdbase.place_particle(self.w, self.C, [5e-6,5e-6,5e-6])
        gfrdbase.place_particle(self.w, self.C, [1e-7,1e-7,1e-7])

        t = self.s.t
        for i in range(5):
            self.s.step()
        self.failIf(t == self.s.t)


    def test_inmobility(self):
        '''Tests whether the created particles are immobile.'''
        particle_a = gfrdbase.place_particle(self.w, self.A, [0.0,0.0,0.0])
        particle_b = gfrdbase.place_particle(self.w, self.B, [1.5000001e-8,0.0,0.0])
        initial_ypos_a = particle_a[1].position
        initial_ypos_b = particle_b[1].position

        for i in range(100): 
            self.s.step()
        
        new_ypos_a = particle_a[1].position
        dist_a = self.w.distance(initial_ypos_a, new_ypos_a)        
        self.failIf(dist_a != 0, 'Particle A (initial: %s, new: %s)' %(initial_ypos_a, new_ypos_a))

        new_ypos_b = particle_b[1].position
        dist_b = self.w.distance(initial_ypos_b, new_ypos_b)
        self.failIf(dist_b != 0, 'Particle B (initial: %s, new: %s)' %(initial_ypos_b, new_ypos_b))


if __name__ == "__main__":
    unittest.main()
