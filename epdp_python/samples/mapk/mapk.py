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
    mpk(bs, phos=YT) | 120
    kk(bs=on) | 30
    pp(bs=on) | 30

@reaction_rules
def rules(ka1, kd1, kf1, ka2, kd2, kf2, tau_rel=0.0):
    (mpk(bs, phos=YT) + kk(bs=on)
        == mpk(bs^1, phos=YT).kk(bs=on^1) | (ka1, kd1)
        > mpk(bs, phos=pYT) + kk(bs=off if tau_rel > 0 else on) | kf1)

    (mpk(bs, phos=pYT) + kk(bs=on)
        == mpk(bs^1, phos=pYT).kk(bs=on^1) | (ka2, kd2)
        > mpk(bs, phos=pYpT) + kk(bs=off if tau_rel > 0 else on) | kf2)

    (mpk(bs, phos=pYpT) + pp(bs=on)
        == mpk(bs^1, phos=pYpT).pp(bs=n^1) | (ka1, kd1)
        > mpk(bs, phos=pYT) + pp(bs=off if tau_rel > 0 else on) | kf1)

    (mpk(bs, phos=pYT) + pp(bs=on)
        == mpk(bs^1, phos=pYT).pp(bs=on^1) | (ka2, kd2)
        > mpk(bs, phos=YT) + pp(bs=off if tau_rel > 0 else on) | kf2)

    # for enz, before, after, ka, kd, kf in (
    #     (kk, YT, pYT, ka1, kd1, kf1), (kk, pYT, pYpT, ka2, kd2, kf2),
    #     (pp, pYpT, pYT, ka1, kd1, kf1), (pp, pYT, YT, ka2, kd2, kf2)):
    #     (mpk(bs, phos=before) + enz(bs=on)
    #         == mpk(bs^1, phos=before).enz(bs=on^1) | (ka, kd)
    #         > mpk(bs, phos=after) + enz(bs=off if tau_rel > 0 else on) | kf)

    if tau_rel > 0:
        krel = numpy.log(2.0) / tau_rel
        kk(bs=off) > kk(bs=on) | krel
        pp(bs=off) > pp(bs=on) | krel


if __name__ == "__main__":
    world_size, matrix_size = 1.0, 5
    D, radius = 1.0, 2.5e-3
    ka1, kd1, kf1 = 0.027 * 1e+24 / N_A, 1.35, 1.5
    ka2, kd2, kf2 = 0.056 * 1e+24 / N_A, 1.73, 15.0
    tau_rel = 1e-6

    m = world.create_model(
        attrs(D, radius),
        rules(ka1, kd1, kf1, ka2, kd2, kf2, tau_rel))
    w = world.create_world(world_size, matrix_size, m, inits())
    sim = simulator.EGFRDSimulator(m, w)

    import sys
    simulator.run_simulation(
        sim, 60.0, 0.1,
        observables=("mpk(phos=YT)", "mpk(phos=pYT)", "mpk(phos=pYpT)"),
        log=sys.stdout)
