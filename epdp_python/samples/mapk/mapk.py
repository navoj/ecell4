import numpy

import ecell4.core
from ecell4.reaction_reader.decorator2 import (
    species_attributes, reaction_rules, molecule_inits, observables,
    create_species, create_reaction_rule)

import ecell4.epdp.simulator as simulator


N_A = 6.0221367e+23

@species_attributes
def attrs(D, radius):
    _ | {"D": D, "radius": radius}

@molecule_inits
def inits():
    mpk(phos=YT) | 120
    kk(bs=on) | 30
    pp(bs=on) | 30

@reaction_rules
def rules(ka1, kd1, kf1, ka2, kd2, kf2, tau_rel=0.0):
    (mpk(phos=YT) + kk(bs=on)
        == mpk(phos=YT^1).kk(bs=on^1) | (ka1, kd1)
        > mpk(phos=pYT) + kk(bs=off if tau_rel > 0 else on) | kf1)

    (mpk(phos=pYT) + kk(bs=on)
        == mpk(phos=pYT^1).kk(bs=on^1) | (ka2, kd2)
        > mpk(phos=pYpT) + kk(bs=off if tau_rel > 0 else on) | kf2)

    (mpk(phos=pYpT) + pp(bs=on)
        == mpk(phos=pYpT^1).pp(bs=n^1) | (ka1, kd1)
        > mpk(phos=pYT) + pp(bs=off if tau_rel > 0 else on) | kf1)

    (mpk(phos=pYT) + pp(bs=on)
        == mpk(phos=pYT^1).pp(bs=on^1) | (ka2, kd2)
        > mpk(phos=YT) + pp(bs=off if tau_rel > 0 else on) | kf2)

    if tau_rel > 0:
        krel = numpy.log(2.0) / tau_rel
        kk(bs=off) > kk(bs=on) | krel
        pp(bs=off) > pp(bs=on) | krel

@observables
def observs():
    mpk(phos=YT) + mpk(phos=YT^_) | K
    mpk(phos=pYT) + mpk(phos=pYT^_) | Kp
    mpk(phos=pYpT) + mpk(phos=pYpT^_) | Kpp


if __name__ == "__main__":
    world_size, matrix_size = 1.0, 5
    D, radius = 1.0, 2.5e-3
    ka1, kd1, kf1 = 0.027 * 1e+24 / N_A, 1.35, 1.5
    ka2, kd2, kf2 = 0.056 * 1e+24 / N_A, 1.73, 15.0
    tau_rel = 1e-6

    sim = simulator.create_simulator(
        world_size, matrix_size,
        attrs(D, radius),
        rules(ka1, kd1, kf1, ka2, kd2, kf2, tau_rel),
        inits())

    import sys
    simulator.run_simulation(
        sim, 60.0, 0.02,
        observables=observs(),
        # log=sys.stdout)
        log="result.dat")
