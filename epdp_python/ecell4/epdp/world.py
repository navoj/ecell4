import numpy

#XXX: epdp modules
import model
import _gfrd
import gfrdbase
#XXX

from ecell4.reaction_reader.decorator2 import create_species, create_reaction_rule
from ecell4.reaction_reader.model import Model


def create_world(world_size, matrix_size):
    x = numpy.repeat(world_size * 0.5, 3)
    world_region = _gfrd.CuboidalRegion('world', _gfrd.Box(x, x))
    if not numpy.all(
        world_region.shape.half_extent == world_region.shape.half_extent[0]):
        raise NotImplementedError("non-cuboidal world is not supported")

    world = _gfrd.World(world_size, matrix_size)
    world.add_structure(world_region)
    return world

class ModelWrapper(Model):

    def __init__(self, *args, **kwargs):
        Model.__init__(self, *args, **kwargs)

        self.__sid_species_map = {}
        self.__model = model.ParticleModel(1.0)

    def add_species(self, sp):
        if self.has_species(sp):
            raise RuntimeError, "already exists [%s]." % (sp.serial())

        st = model.Species(
            sp.serial(), sp.get_attribute("D"), sp.get_attribute("radius"))
        self.__model.add_species_type(st)
        self.__sid_species_map[st.id] = sp
        return st.id

    def has_species(self, sp):
        return (self.get_species_id(sp) is not None)

    def num_species(self):
        return len(self.__sid_species_map)

    def get_species_id(self, sp):
        for sid, sp2 in self.__sid_species_map.items():
            if sp2 == sp:
                return sid
        else:
            return None

    def get_species(self, sid):
        return self.__sid_species_map.get(sid)

    def query_reaction_rules(self, sp1, sp2=None):
        if sp2 is None:
            return None
        else:
            rr = _gfrd.ReactionRuleInfo(0, 0.0, [sp1, sp2], [])
            return [rr]

    query_reaction_rule = query_reaction_rules

    def has_species_attribute(self, sp):
        return Model.has_species_attribute(self, create_species(sp.serial()))

    def apply_species_attributes(self, sp):
        sp1 = create_species(sp.serial())
        Model.apply_species_attributes(self, sp1)

        for key, value in sp1.attributes().items():
            sp.set_attribute(key, value)

class World:

    def __init__(self, world_size, matrix_size):
        self.world_size = world_size
        self.matrix_size = matrix_size

        # self.model.network_rules
        # self.model.get_species_type_by_id(sid)
        # self.structures

        self.model = ModelWrapper()

        w = create_world(self.world_size, self.matrix_size)
        self.world = w

    def get_num_particles(self):
        return self.world.num_particles

    num_particles = property(get_num_particles)

    def get_structures(self):
        return self.world.structures

    structures = property(get_structures)

    #XXX: HERE

    def get_species(self, sid):
        return self.world.get_species(sid)

    def new_particle(self, sid, pos):
        return self.world.new_particle(sid, pos)

    def update_particle(self, pid_particle_pair):
        return self.world.update_particle(pid_particle_pair)

    #XXX: THERE

    def remove_particle(self, sid):
        self.world.remove_particle(sid)

    def check_overlap(self, region, ignore):
        pos, radius = region
        return self.world.check_overlap((pos, radius), ignore)

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

    def ecell4__add_species(self, sp):
        self.model.apply_species_attributes(sp)
        sid = self.model.add_species(sp)

        if sp.has_attribute("structure"):
            structure = st["structure"]
        else:
            structure = "world"

        self.world.add_species(
            _gfrd.SpeciesInfo(
                sid,  float(sp.get_attribute("D")),
                float(sp.get_attribute("radius")), structure, 0.0))
        return sid

    def ecell4__has_species(self, sp):
        return self.model.has_species(sp)

    def ecell4__add_molecules(self, sp, num):
        if not self.ecell4__has_species(sp):
            self.ecell4__add_species(sp)
        sid = self.model.get_species_id(sp)
        gfrdbase.throw_in_particles(self.world, sid, num)

    def ecell4__num_particles(self, sp=None):
        if sp is None:
            return self.world.num_particles
        elif self.ecell4__has_species(sp):
            sid = self.model.get_species_id(sp)
            return len(self.world.get_particle_ids(sid))
        else:
            return 0

class EGFRDWorld:

    def __init__(self, world_size, matrix_size, rng=None):
        self.world = World(world_size, matrix_size)
        self.t = 0.0

        if rng is not None:
            raise RuntimeError
        else:
            import myrandom
            self.internal_rng = myrandom.rng

    def bind_to(self, m):
        if self.world.model.num_species() > 0:
            raise RuntimeError
        self.world.model = m

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

    def num_molecules(self, *args, **kwargs):
        return self.num_particles(*args, **kwargs)
