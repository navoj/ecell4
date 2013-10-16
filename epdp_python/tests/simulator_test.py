import unittest

import ecell4.core

import world
import simulator


class SimulatorTestCase(unittest.TestCase):

    def setUp(self):
        pass

    def test1(self):
        w = world.EGFRDWorld(1.0, 3)

        self.assertEqual(w.num_particles(), 0)
        self.assertEqual(
            w.world.apply_boundary((1.5, 0.5, 0.5))[0], 0.5)

        sp1 = ecell4.core.Species("X")
        sp1.set_attribute("D", "1.0")
        sp1.set_attribute("radius", "2.5e-3")
        species_list = [sp1]

        w.add_molecules(sp1, 60)
        self.assertEqual(len(w.world.world.species), 1)
        self.assertTrue(w.has_species(sp1))
        self.assertEqual(w.num_particles(), 60)

        sim = simulator.EGFRDSimulator(None, w)

        next_time, dt = 0.0, 0.02

        import sys
        import csv

        writer = csv.writer(sys.stdout, delimiter='\t')
        writer.writerow(['#t'] + [sp.name() for sp in species_list])
        writer.writerow(
            ['%.6e' % sim.t()]
            + ['%d' % w.num_molecules(sp) for sp in species_list])
        for i in xrange(5):
        # for i in xrange(10000):
            next_time += dt
            while sim.step(next_time):
                pass

            writer.writerow(
                ['%.6e' % sim.t()]
                + ['%d' % w.num_molecules(sp) for sp in species_list])
            # print list(w.world)[0][1].position


if __name__ == "__main__":
    unittest.main()
