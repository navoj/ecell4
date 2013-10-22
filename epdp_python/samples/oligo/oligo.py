import numpy

import ecell4.core
from ecell4.reaction_reader.decorator2 import (
    species_attributes, reaction_rules, molecule_inits,
    create_species, create_reaction_rule)

import ecell4.epdp.world as world
import ecell4.epdp.simulator as simulator


N_A = 6.0221367e+23

@species_attributes
def attrs(D, radius):
    _ | {"D": D, "radius": radius}

@molecule_inits
def inits():
    A(l, r, bs=(l, r)) | 120

@reaction_rules
def rules(ka, kd):
    _(bs) + _(bs) == _(bs^1)._(bs^1) | (ka, kd)


if __name__ == "__main__":
    world_size, matrix_size = 1.0, 5
    D, radius = 1.0, 2.5e-3
    kD = 4 * numpy.pi * (radius + radius) * (D + D)
    ka, kd = kD, 0.0

    sim = simulator.create_simulator(
        world_size, matrix_size,
        attrs(D, radius), rules(ka, kd), inits())

    import sys

    # with open('oligo.dat', 'w') as fout:
    with sys.stdout as fout:
        species_set = set(sim.world.species)
        fout.write("%e, %3d, " % (sim.t(), sim.world.num_particles()))
        fout.write(", ".join([sp.name() for sp in species_set]))
        fout.write("\n")
        while sim.world.num_particles() > 1:
            sim.step()
            if len(species_set) != len(sim.world.species):
                fout.write("%e, %3d, " % (sim.t(), sim.world.num_particles()))
                fout.write(", ".join([sp.name() for sp in (set(sim.world.species) - species_set)]))
                fout.write("\n")
                species_set = set(sim.world.species)
