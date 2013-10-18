import unittest

import ecell4.core
from ecell4.reaction_reader.decorator2 import create_species, create_reaction_rule

import world
import simulator

import math


class SimulatorTestCase(unittest.TestCase):

    def setUp(self):
        pass

    def test1(self):
        m = world.ModelWrapper()

        D, radius = "1.0", "2.5e-3"
        attr1 = create_species("_")
        attr1.set_attribute("D", D)
        attr1.set_attribute("radius", radius)
        m.add_species_attribute(attr1)

        kon, koff = 0.1, 3.0
        kD = (4 * math.pi * (float(radius) + float(radius))
            * (float(D) + float(D)))
        ka = kon * kD / (kD - kon)
        kd = ka * koff / kon
        m.add_reaction_rule(create_reaction_rule("X+Y>Z|%e" % ka))
        m.add_reaction_rule(create_reaction_rule("Z>X+Y|%e" % kd))

        w = world.EGFRDWorld(1.0, 3)
        w.bind_to(m)

        self.assertEqual(w.num_particles(), 0)
        self.assertEqual(
            w.world.apply_boundary((1.5, 0.5, 0.5))[0], 0.5)

        sp1 = ecell4.core.Species("X")
        sp2 = ecell4.core.Species("Y")
        w.add_molecules(sp1, 60)
        w.add_molecules(sp2, 60)
        self.assertEqual(len(w.world.world.species), 2)
        self.assertTrue(w.has_species(sp1))
        self.assertEqual(w.num_particles(), 120)

        sim = simulator.EGFRDSimulator(m, w)

        next_time, dt = 0.0, 0.02
        species_list = [sp1, sp2]

        import sys
        import csv

        with sys.stdout as fout:
        # with open('test.out', 'w') as fout:
            writer = csv.writer(fout, delimiter='\t')
            writer.writerow(['#t'] + [sp.name() for sp in species_list])
            writer.writerow(
                ['%.6e' % sim.t()]
                + ['%d' % w.num_molecules(sp) for sp in species_list]
                + ['%d' % w.num_molecules()])
            for i in xrange(500):
                next_time += dt
                while sim.step(next_time):
                    pass

                writer.writerow(
                    ['%.6e' % sim.t()]
                    + ['%d' % w.num_molecules(sp) for sp in species_list]
                    + ['%d' % w.num_molecules()])


if __name__ == "__main__":
    unittest.main()
