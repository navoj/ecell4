#XXX: epdp modules
import _gfrd
import egfrd
#XXX


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
