import logging
import unittest
import numpy
import _gfrd

from egfrd import *
from visualization import vtklogger

import model
import gfrdbase
import myrandom


class EGFRDSimulatorTestCase(unittest.TestCase):

    def setUp(self):
        self.m = model.ParticleModel(1e-5)

        self.S = model.Species('S', 2e-11, 5e-8)
        self.SS = model.Species('SS', 1e-12, 5e-9)
        self.A = model.Species('A', 0, 1e-8)
        self.B = model.Species('B', 2e-11, 5e-9)
        self.m.add_species_type(self.S)
        self.m.add_species_type(self.SS)
        self.m.add_species_type(self.A)
        self.m.add_species_type(self.B)
        self.m.set_all_repulsive()
        self.w = gfrdbase.create_world(self.m)
        self.nrw = _gfrd.NetworkRulesWrapper(self.m.network_rules)
        self.s = EGFRDSimulator(self.w, myrandom.rng, self.nrw)

    def tearDown(self):
        pass
    '''
    def test_instantiation(self):
        self.failIf(self.s == None)


    def test_one_particle(self):
        place_particle(self.s.world, self.S, [0.0,0.0,0.0])

        t = self.s.t
        for i in range(5):
            self.s.step()
        self.failIf(t == self.s.t)
    
    def test_two_particles(self):
        place_particle(self.s.world, self.S, [0.0,0.0,0.0])
        place_particle(self.s.world, self.S, [5e-6,5e-6,5e-6])

        t = self.s.t
        for i in range(5):
            self.s.step()
        self.failIf(t == self.s.t)

    def test_three_particles(self):
        place_particle(self.s.world, self.S, [0.0,0.0,0.0])
        place_particle(self.s.world, self.S, [5e-6,5e-6,5e-6])
        place_particle(self.s.world, self.S, [1e-7,1e-7,1e-7])

        t = self.s.t
        for i in range(5):
            self.s.step()
        self.failIf(t == self.s.t)


    def test_three_particles_in_contact(self):
        place_particle(self.s.world, self.S, [0.0,0.0,0.0])
        place_particle(self.s.world, self.S, [1e-7,0.0,0.0])

        # dummy
        place_particle(self.s.world, self.S, [2e-7,0.0,0.0])

        r = model.create_unimolecular_reaction_rule(self.S, self.A, 1e10)
        self.m.network_rules.add_reaction_rule(r)

        t = self.s.t
        for i in range(5):
            self.s.step()
            
            # Check if species ids are consistent after unimolecular 
            # multi reaction.
            for species in self.s.world.species:
                for pid in self.s.world.get_particle_ids(species.id):
                    particle = self.s.world.get_particle(pid)[1]
                    self.failIf(particle.sid != species.id)
        self.failIf(t == self.s.t)

    def test_four_particles_close(self):
        place_particle(self.s.world, self.SS, [2e-8,0.0,0.0])
        place_particle(self.s.world, self.SS, [3.003e-8,0.0,0.0])

        place_particle(self.s.world, self.SS, [0.994e-8-5e-10, 0.0, 0.0])

        t = self.s.t
        for i in range(10):
            self.s.step()
        self.failIf(t == self.s.t)

    def test_immobile_is_immobile(self):
        particleA = place_particle(self.s.world, self.A, [0.0,0.0,0.0])
        place_particle(self.s.world, self.B, [1.5000001e-8,0.0,0.0])

        initial_position = particleA[1].position

        for i in range(10):
            self.s.step()
        
        new_position = particleA[1].position
        dist = self.w.distance(initial_position, new_position)

        self.failIf(dist != 0, 'initial pos: %s,\tnew pos: %s' %
                    (initial_position, new_position))

    def test_pair_with_immobile(self):
        place_particle(self.s.world, self.A, [0.0, 0.0, 0.0])
        place_particle(self.s.world, self.B, [1.51e-8, 0.0, 0.0])

        for i in range(2):
            self.s.step()

    def test_pair_with_immobile_switched_order(self):
        place_particle(self.s.world, self.B, [1.51e-8, 0.0, 0.0])
        place_particle(self.s.world, self.A, [0.0, 0.0, 0.0])

        for i in range(2):
            self.s.step()
'''


class EGFRDSimulatorTestCaseBase(unittest.TestCase):
    """Base class for TestCases below.

    """
    def create_model(self):
        self.L = 1.e-6

        self.D = 1.e-12
        self.radius = 5.e-9

        self.m = model.ParticleModel(self.L)

        self.surface = _gfrd.StructureType()
        self.surface["name"]='surface'
        self.m.add_structure_type(self.surface)

        self.cylinder = _gfrd.StructureType()
        self.cylinder["name"]='cylinder'
        self.m.add_structure_type(self.cylinder)

        self.A = model.Species('A', self.D, self.radius)
        self.B = model.Species('B', self.D, self.radius, self.cylinder)
        self.C = model.Species('C', self.D, self.radius)

        self.kf_1 = 4000
        self.kf_2 = 5e-19
        self.kb_1 = 4000
        self.kb_2 = 4000


    def add_planar_surface(self):
        m1 = model.create_planar_surface(self.surface.id,'m1', [0, 0, 0], [1, 0, 0], [0, 1, 0], self.L, self.L, self.parent_id)
        self.w.add_structure(m1)


    def add_cylindrical_surface(self):
        d = model.create_cylindrical_surface(self.cylinder.id,'d', [self.L / 2, 0, self.L / 2], self.radius, [0, 1, 0], self.L, self.parent_id)
        self.w.add_structure(d)


    def add_species(self, domain_a=None, domain_b=None, domain_c=None):
        if domain_a is None : self.A = model.Species('A', self.D, self.radius)
        else: self.A = model.Species('A', self.D, self.radius, domain_a)

        if domain_b is None : self.B = model.Species('B', self.D, self.radius)
        else: self.B = model.Species('B', self.D, self.radius, domain_b)
        
        if domain_c is None : self.C = model.Species('C', self.D, self.radius)
        else: self.C = model.Species('C', self.D, self.radius, domain_c)

        self.m.add_species_type(self.A)
        self.m.add_species_type(self.B)
        self.m.add_species_type(self.C)


    def create_simulator(self):
        self.w = create_world(self.m)
        self.parent_id=self.w.get_def_structure_id()
        self.s = EGFRDSimulator(self.w)


    def add_reactions(self):
        A = self.A
        B = self.B
        C = self.C

        r = model.create_unimolecular_reaction_rule(A, B, self.kf_1)
        self.m.network_rules.add_reaction_rule(r)
        r = model.create_unimolecular_reaction_rule(B, A, self.kb_1)
        self.m.network_rules.add_reaction_rule(r)
        r = model.create_binding_reaction_rule(A, B, C, self.kf_2)
        self.m.network_rules.add_reaction_rule(r)
        r = model.create_unbinding_reaction_rule(C, A, B, self.kb_2)
        self.m.network_rules.add_reaction_rule(r)
        r = model.create_decay_reaction_rule(C, self.kb_1)
        self.m.network_rules.add_reaction_rule(r)

    def add_particles(self, n):
        throw_in_particles(self.w, self.A, n)
        throw_in_particles(self.w, self.B, n)

    def tearDown(self):
        pass


class CytosoleTestCase(EGFRDSimulatorTestCaseBase):
    """Events happening in the "world".

    """
    def setUp(self):
        self.create_model()

        self.add_species() 
        self.create_simulator() 

        self.add_reactions()
        self.add_particles(2)

    def test_run(self):
        for i in range(10):
            self.s.step()

    def test_vtklogger(self):
        vtk_logger = vtklogger.VTKLogger(self.s, 'vtk_temp_data')
        for i in range(10):
            vtk_logger.log()
            self.s.step()
        vtk_logger.stop()
        vtk_logger.cleanup()


class PlanarSurfaceTestCase(EGFRDSimulatorTestCaseBase):
    """Events happening *on* a planar surface.

    """
    def setUp(self):
        self.create_model()

        # All species on planar surface.
        self.A["structure"] = "m1"
        self.B["structure"] = "m1"
        self.C["structure"] = "m1"
        self.add_species(self.surface, self.surface, self.surface)
        self.create_simulator() 

        self.add_planar_surface()
        self.add_reactions()
        self.add_particles(2)


    def test_run(self):
        for i in range(10):
            self.s.step()


    def test_vtklogger(self):
        vtk_logger = vtklogger.VTKLogger(self.s, 'vtk_temp_data')
        for i in range(10):
            vtk_logger.log()
            self.s.step()
        vtk_logger.stop()
        vtk_logger.cleanup()


class CylindricalSurfaceTestCase(EGFRDSimulatorTestCaseBase):
    """Events happening *on* a cylindrical surface.

    """
    def setUp(self):
        self.create_model()
        
        # All species on cylindrical surface.
        self.A["structure"] = "d"
        self.B["structure"] = "d"
        self.C["structure"] = "d"
        self.add_species(self.cylinder, self.cylinder, self.cylinder)
        self.create_simulator() 

        self.add_cylindrical_surface()
        self.add_reactions()
        self.add_particles(2)


    def test_run(self):
        for i in range(10):
            self.s.step()


    def test_vtklogger(self):
        vtk_logger = vtklogger.VTKLogger(self.s, 'vtk_temp_data')
        for i in range(10):
            vtk_logger.log()
            self.s.step()
        vtk_logger.stop()
        vtk_logger.cleanup()


class PlanarSurfaceInteractionTestCase(EGFRDSimulatorTestCaseBase):
    """Events between the "world" and a planar surface.

    """
    def setUp(self):
        self.create_model()

        # Only species B is on the planar surface.
        self.B["structure"] = "m1"
        self.add_species(None, self.surface, None) 
        self.create_simulator() 
        
        self.add_planar_surface()
        self.add_particles(10)
        # Don't add all reactions.

    def test_interaction_single_is_formed(self):
        # Place a particle very close to the planar surface.
        z_position = float(self.A['radius']) * (1 + 1e-1)
        place_particle(self.w, self.A, [0.0,0.0, z_position])

        for i in range(10):
            self.s.step()


class CylindricalSurfaceInteractionTestCase(EGFRDSimulatorTestCaseBase):
    """Events between the "world" and a cylindrical surface.

    """
    def setUp(self):
        self.create_model()

        # Only species B is on the cylindrical surface.
        self.B["structure"] = "d"
        self.add_species(None, self.cylinder, None) 
        self.create_simulator() 
        # Don't add particles yet.
        # Don't add all reactions.
        self.add_cylindrical_surface()


    def test_simulation_does_not_come_to_a_halt(self):
        # Place a particle on the cylindrical surface.
        L = self.L
        place_particle(self.w, self.B, [-L / 2, -L / 2, -L / 2])

        # Place a particle very close to the cylindrical surface.
        offset = 5.0*self.radius
        place_particle(self.w, self.A,[L/2, L/2+offset, L/2+offset])

        vtk_logger = vtklogger.VTKLogger(self.s, 'vtk_temp_data')
        for i in range(10):
            self.s.step()
            vtk_logger.log()
        vtk_logger.stop()
        vtk_logger.cleanup()

        self.failIf(self.s.t < 1.e-10)
   
    
    def test_no_pair_between_species_on_different_surfaces(self):
        myrandom.seed(1)
        # Place a particle on the cylindrical surface.
        L = self.L
        place_particle(self.w, self.B, [L/2, L/2, L/2])

        # Place a particle very close to the cylindrical surface.
        offset = 5.0*self.radius
        place_particle(self.w, self.A,[L/2, L/2, L/2+offset])

        vtk_logger = vtklogger.VTKLogger(self.s, 'vtk_temp_data')

        # After 24 steps multi particle tries a move that would result 
        # in overlap with surface.
        #self.s.bd_dt_factor = 5
        for i in range(25):
            self.s.step()
            vtk_logger.log()
            self.s.check()
        vtk_logger.stop()
        vtk_logger.cleanup()
        #self.failIf(numpy.array(self.s.pair_steps.values()).sum() > 0)


    def test_no_overlap_interaction_and_cylindrical_single(self):
        myrandom.seed(1)
        # Place a particle on the cylindrical surface.
        L = self.L
        place_particle(self.w, self.B, [L/2, L/2, L/2])

        # Place a particle very close to the cylindrical surface.
        offset = 5.0*self.radius
        place_particle(self.w, self.A, [L/2, L/2, L/2+offset])

        vtk_logger = vtklogger.VTKLogger(self.s, 'vtk_temp_data')
        for i in range(10):
            self.s.step()
            vtk_logger.log()
            self.s.check()
        vtk_logger.stop()
        vtk_logger.cleanup()
    

    def test_longer_run(self):
        self.add_particles(20)

        for i in range(100):
            self.s.step()


if __name__ == "__main__":    
    unittest.main()

