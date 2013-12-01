import unittest

import ecell4.core

import world
import simulator


class WorldTestCase(unittest.TestCase):

    def setUp(self):
        pass

    def test1(self):
        w = world.EGFRDWorld(1.0, 3)

        self.assertEqual(w.num_particles(), 0)
        self.assertEqual(
            w.world.apply_boundary((1.5, 0.5, 0.5))[0], 0.5)

        sp1 = ecell4.core.Species("X")
        sp1.set_attribute("D", "1e-12")
        sp1.set_attribute("radius", "2.5e-9")
        species_list = [sp1]

        w.add_molecules(sp1, 60)
        self.assertEqual(len(w.world.species), 1)
        self.assertTrue(w.has_species(sp1))
        self.assertEqual(w.num_particles(), 60)

    def test2(self):
        w = world.EGFRDWorld(1.0, 3)

        self.assertEqual(len(w.world.species), 0)
        sp1 = ecell4.core.Species("X")
        sp1.set_attribute("D", "1e-12")
        sp1.set_attribute("radius", "2.5e-9")
        sid1 = w.add_species(sp1)
        self.assertEqual(len(w.world.species), 1)
        w.world.remove_species(sid1)
        self.assertEqual(len(w.world.species), 0)
        sid2 = w.add_species(sp1)
        self.assertEqual(len(w.world.species), 1)
        self.assertNotEqual(sid1, sid2)


if __name__ == "__main__":
    unittest.main()
