import csv
import collections
import numpy

import ecell4.core

#XXX: epdp modules
import _gfrd
import egfrd
#XXX


def create_simulator(world_size, matrix_size, attrs, rules, inits):
    import world
    m = world.create_model(attrs, rules)
    w = world.create_world(world_size, matrix_size, m, inits)
    return EGFRDSimulator(m, w)

class TimecourseLogger:

    def __init__(self, func=None):
        self.__func = func
        self.__cache = []

    def __call__(self, *args, **kwargs):
        if self.__func is not None:
            self.__func(*args, **kwargs)

        self.log(*args)

    def log(self, *args):
        self.__cache.append(args)

    def get(self):
        return self.__cache

def run_simulation(sim, upto, dt, observables, log=None):
    if type(log) in (file, str):
        if type(log) is str:
            writer = csv.writer(open(log, 'w'), delimiter='\t')
        else:
            writer = csv.writer(log, delimiter='\t')
        def logfunc(t, *num_molecules):
            writer.writerow(['%.6e' % t] + ['%d' % n for n in num_molecules])
        writer.writerow(['#t'] + [rhs for lhs, rhs in observables])
        logobj = TimecourseLogger(logfunc)
    else:
        logobj= TimecourseLogger(log)

    next_time, retval = 0.0, []

    species_list = []
    for lhs, rhs in observables:
        if isinstance(lhs, collections.Iterable):
            for sp in lhs:
                sp.sort()
            species_list.append(
                tuple([ecell4.core.Species(str(sp)) for sp in lhs]))
        else:
            lhs.sort()
            species_list.append((ecell4.core.Species(str(lhs)), ))

    def num_molecules(sp):
        return sum([sim.world.num_molecules(elem) for elem in sp])

    logobj(sim.t(), *[num_molecules(sp) for sp in species_list])
    while sim.t() < upto:
        next_time += dt
        while sim.step(next_time):
            pass

        logobj(sim.t(), *[num_molecules(sp) for sp in species_list])
    return logobj.get()

class EGFRDSimulator:

    def __init__(self, m, w):
        self.model = m
        self.world = w

        nrw = m
        # nrw = _gfrd.NetworkRulesWrapper(m.network_rules)
        self.sim = egfrd.EGFRDSimulator(w.world, w.internal_rng, nrw)

    def t(self):
        return self.sim.t

    def dt(self):
        return self.sim.dt

    def num_steps(self):
        return self.sim.step_counter

    def step(self, upto=None):
        if upto is None:
            if self.world.num_particles() == 0:
                return

            self.sim.step()
            self.world.set_t(self.sim.t)
        else:
            if self.world.num_particles() == 0:
                return True
            elif upto <= self.sim.t:
                return False
            elif upto >= self.t() + self.dt():
                self.step()
                self.world.set_t(self.sim.t)
                return True
            else:
                self.sim.stop(upto)
                self.world.set_t(self.sim.t)
                return False

    def set_t(self, t):
        raise NotImplementedError

    def initialize(self):
        self.sim.initialize()
