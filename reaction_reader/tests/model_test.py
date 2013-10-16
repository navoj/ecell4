import unittest

import model
from decorator2 import create_species, create_reaction_rule


class ModelTestCase(unittest.TestCase):

    def setUp(self):
        pass

    def test1(self):
        m = model.Model()

        sp1 = create_species("_(bs)")
        sp1.set_attribute("D", "1")
        m.add_species_attribute(sp1)

        sp2 = create_species("_(bs^_)")
        sp2.set_attribute("D", "0")
        m.add_species_attribute(sp2)

        sp3 = create_species("A(bs^1).A(bs^1)")

        self.assertTrue(m.has_species_attribute(create_species("A(bs)")))
        self.assertTrue(m.has_species_attribute(sp3))
        self.assertFalse(m.has_species_attribute(create_species("A")))

        m.apply_species_attributes(sp3)
        self.assertTrue(sp3.has_attribute("D"))
        self.assertEqual(sp3.get_attribute("D"), "0")

    def test2(self):
        m = model.Model()

        m.add_reaction_rule(
            create_reaction_rule("_(bs)+_(bs)>_(bs^1)._(bs^1)|1.0"))
        m.add_reaction_rule(
            create_reaction_rule("_(bs^1)._(bs^1)>_(bs)+_(bs)|1.0"))

        sp1 = create_species("A(bs)")
        sp2 = create_species("A(bs^1).A(bs^1)")

        self.assertEqual(len(m.query_reaction_rules(sp1, sp2)), 0)

        retval = m.query_reaction_rules(sp1, sp1)
        self.assertEqual(len(retval), 2)
        self.assertEqual(len(retval[0]), 3)
        self.assertEqual(len(retval[0][1]), 1)
        self.assertEqual(retval[0][1][0], sp2)

        print retval[0]

        retval = m.query_reaction_rules(sp2)
        self.assertEqual(len(retval), 2)
        self.assertEqual(len(retval[0]), 3)
        self.assertEqual(len(retval[0][1]), 2)
        self.assertEqual(retval[0][1][0], sp1)
        self.assertEqual(retval[0][1][1], sp1)


if __name__ == '__main__':
    unittest.main()
