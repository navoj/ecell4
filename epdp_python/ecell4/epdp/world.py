import numpy

#XXX: epdp modules
import _gfrd
import gfrdbase
#XXX

import ecell4.core
from ecell4.reaction_reader.decorator2 import create_species, create_reaction_rule
from ecell4.reaction_reader.model import Model


def create_model(attrs, rules):
    m = ModelWrapper()

    for attr in attrs:
        m.add_species_attribute(attr)

    for rr in rules:
        m.add_reaction_rule(rr)

    return m

def create_world(world_size, matrix_size, model=None, inits=[]):
    w = EGFRDWorld(world_size, matrix_size)

    if model is not None:
        w.bind_to(model)

    for sp, num in inits:
        sp.sort()
        newsp = ecell4.core.Species(str(sp))
        if model is not None:
            model.apply_species_attributes(newsp)
        w.add_molecules(newsp, num)

    return w

class ModelWrapper(Model):

    def __init__(self, *args, **kwargs):
        Model.__init__(self, *args, **kwargs)

        self.__sid_species_map = {}
        self.__sidgen = _gfrd.SpeciesIDGenerator(0)

    def add_species(self, sp):
        if self.has_species(sp):
            raise RuntimeError, "already exists [%s]." % (sp.serial())

        sid = self.__sidgen()
        self.__sid_species_map[sid] = sp
        return sid

    def has_species(self, sp):
        return (self.get_species_id(sp) is not None)

    def num_species(self):
        return len(self.__sid_species_map)

    def __get_species_id(self, serial):
        for sid, sp2 in self.__sid_species_map.items():
            if sp2.serial() == serial:
                return sid
        else:
            return None

    def get_species_id(self, sp):
        return self.__get_species_id(sp.serial())

    def get_species(self, sid):
        return self.__sid_species_map.get(sid)

    def query_reaction_rules(self, sid1, sid2=None):
        if sid2 is None:
            sp1 = create_species(self.get_species(sid1).serial())
            rrs = Model.query_reaction_rules(self, sp1)
            reactants = [sid1]
        else:
            sp1 = create_species(self.get_species(sid1).serial())
            sp2 = create_species(self.get_species(sid2).serial())
            rrs = Model.query_reaction_rules(self, sp1, sp2)
            reactants = [sid1, sid2]

            if len(rrs) == 0:
                rrs.append(((sp1, sp2), (), 0.0))
            # rr = _gfrd.ReactionRuleInfo(0, 0.0, [sid1, sid2], [])

        retval = []
        for rr in rrs:
            _, products, k = rr
            product_ids = []
            for sp in products:
                sp.sort()
                newsp = ecell4.core.Species(str(sp))
                sid = self.get_species_id(newsp)
                if sid is None:
                    self.apply_species_attributes(newsp)
                    sid = self.add_species(newsp)
                product_ids.append(sid)
            retval.append(
                _gfrd.ReactionRuleInfo(0, k, reactants, product_ids))
        return retval

    query_reaction_rule = query_reaction_rules

    def has_species_attribute(self, sp):
        return Model.has_species_attribute(self, create_species(sp.serial()))

    def apply_species_attributes(self, sp):
        sp1 = create_species(sp.serial())
        Model.apply_species_attributes(self, sp1)

        for key, value in sp1.attributes().items():
            sp.set_attribute(key, value)

class DummyModel:

    def __init__(self):
        pass

    def get_species_type_by_id(self, sid):
        return dict(name=str(sid))

class World:

    def __init__(self, world_size, matrix_size):
        self.world_size = world_size
        self.matrix_size = matrix_size

        # self.model.network_rules
        # self.structures

        self.__model = ModelWrapper()

        x = numpy.repeat(world_size * 0.5, 3)
        world_region = _gfrd.CuboidalRegion('world', _gfrd.Box(x, x))
        if not numpy.all(
            world_region.shape.half_extent == world_region.shape.half_extent[0]):
            raise NotImplementedError("non-cuboidal world is not supported")

        self.world = _gfrd.World(world_size, matrix_size)
        self.world.add_structure(world_region)

        self.world.model = DummyModel() # just for debugging
        self.model = self.world.model

    def get_num_particles(self):
        return self.world.num_particles

    num_particles = property(get_num_particles)

    def get_structures(self):
        return self.world.structures

    structures = property(get_structures)

    #XXX: HERE

    def has_species(self, sid):
        return sid in [sp.id for sp in self.world.species]

    def get_species(self, sid):
        if not self.has_species(sid):
            sp = self.__model.get_species(sid)
            if sp is None:
                raise RuntimeError, sid
            self.ecell4__add_species_to_world(sid, sp)

        return self.world.get_species(sid)

    def new_particle(self, sid, pos):
        return self.world.new_particle(sid, pos)

    def update_particle(self, pid_particle_pair):
        return self.world.update_particle(pid_particle_pair)

    #XXX: THERE

    def remove_particle(self, sid):
        self.world.remove_particle(sid)

    def check_overlap(self, region, ignore1, ignore2=None):
        pos, radius = region
        if ignore2 is None:
            return self.world.check_overlap((pos, radius), ignore1)
        else:
            return self.world.check_overlap((pos, radius), ignore1, ignore2)

    def __iter__(self):
        return self.world.__iter__()

    def get_structure(self, sid):
        return self.world.get_structure(sid)

    def apply_boundary(self, pos):
        return self.world.apply_boundary(pos)

    def cyclic_transpose(self, pos1, pos2):
        return self.world.cyclic_transpose(pos1, pos2)

    def calculate_pair_CoM(self, pos1, pos2, D1, D2):
        return self.world.calculate_pair_CoM(pos1, pos2, D1, D2)

    def distance(self, pos1, pos2):
        return self.world.distance(pos1, pos2)

    def list_species(self):
        return self.world.species

    species = property(list_species)

    #XXX: E-Cell4 special functions

    def ecell4__add_species_to_world(self, sid, sp):
        if sp.has_attribute("structure"):
            structure = sp.get_attribute("structure")
        else:
            structure = "world"

        self.world.add_species(
            _gfrd.SpeciesInfo(
                sid,  float(sp.get_attribute("D")),
                float(sp.get_attribute("radius")), structure, 0.0))

    def ecell4__bind_to(self, m):
        if self.__model.num_species() > 0:
            raise RuntimeError
        self.__model = m

    def ecell4__add_species(self, sp):
        # self.__model.apply_species_attributes(sp)
        sid = self.__model.add_species(sp)
        self.ecell4__add_species_to_world(sid, sp)
        return sid

    def ecell4__has_species(self, sp):
        return self.__model.has_species(sp)

    def ecell4__add_molecules(self, sp, num):
        if not self.ecell4__has_species(sp):
            self.ecell4__add_species(sp)
        sid = self.__model.get_species_id(sp)
        gfrdbase.throw_in_particles(self.world, sid, num)

    def ecell4__num_particles(self, sp=None):
        if sp is None:
            return self.world.num_particles
        elif self.ecell4__has_species(sp):
            sid = self.__model.get_species_id(sp)
            return len(self.world.get_particle_ids(sid))
        else:
            return 0

    def ecell4__num_molecules(self, sp=None):
        if sp is None:
            sp = ecell4.core.Species("_")
        newsp1 = create_species(sp.serial())
        num = 0
        for sp2 in self.world.species:
            retval = newsp1.match(
                create_species(self.__model.get_species(sp2.id).serial()))
            if len(retval) > 0:
                num += len(self.world.get_particle_ids(sp2.id)) * len(retval)
        return num

    def ecell4__get_species(self, sid):
        return self.__model.get_species(sid)

class EGFRDWorld:

    def __init__(self, world_size, matrix_size, rng=None):
        self.world = World(world_size, matrix_size)
        self.t = 0.0

        if rng is not None:
            raise RuntimeError
        else:
            import myrandom
            self.internal_rng = myrandom.rng

    def list_species(self):
        #XXX: use a generator instead
        retval = []
        for sp in self.world.species:
            retval.append(self.world.ecell4__get_species(sp.id))
        return retval

    species = property(list_species)

    def bind_to(self, m):
        self.world.ecell4__bind_to(m)

    def t(self):
        return self.t

    def set_t(self, t):
        self.t = t

    def add_species(self, sp):
        self.world.ecell4__add_species(sp)

    def has_species(self, sp):
        return self.world.ecell4__has_species(sp)

    def add_molecules(self, sp, num):
        return self.world.ecell4__add_molecules(sp, num)

    def num_particles(self, sp=None):
        return self.world.ecell4__num_particles(sp)

    def num_molecules(self, sp=None):
        # return self.num_particles(*args, **kwargs)
        return self.world.ecell4__num_molecules(sp)
