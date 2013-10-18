import unittest

import ecell4.core
from ecell4.reaction_reader.decorator2 import create_species, create_reaction_rule

import world
import simulator


class SimulatorTestCase(unittest.TestCase):

    def setUp(self):
        pass

    def test1(self):
        m = world.ModelWrapper()
        attr1 = create_species("_")
        attr1.set_attribute("D", "1.0")
        attr1.set_attribute("radius", "2.5e-3")
        m.add_species_attribute(attr1)

        m.add_reaction_rule(create_reaction_rule("X+X>Y|0.4896553344263186"))
        m.add_reaction_rule(create_reaction_rule("Y>X+X|6.0"))

        w = world.EGFRDWorld(1.0, 3)
        w.bind_to(m)

        self.assertEqual(w.num_particles(), 0)
        self.assertEqual(
            w.world.apply_boundary((1.5, 0.5, 0.5))[0], 0.5)

        sp1 = ecell4.core.Species("X")
        w.add_molecules(sp1, 60)
        self.assertEqual(len(w.world.world.species), 1)
        self.assertTrue(w.has_species(sp1))
        self.assertEqual(w.num_particles(), 60)

        sim = simulator.EGFRDSimulator(m, w)

        next_time, dt = 0.0, 0.02
        species_list = [sp1]

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
            # for i in xrange(100):
            for i in xrange(500):
                next_time += dt
                while sim.step(next_time):
                    pass

                writer.writerow(
                    ['%.6e' % sim.t()]
                    + ['%d' % w.num_molecules(sp) for sp in species_list]
                    + ['%d' % w.num_molecules()])
                # print list(w.world)[0][1].position


if __name__ == "__main__":
    unittest.main()
