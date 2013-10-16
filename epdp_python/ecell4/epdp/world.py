import numpy

#XXX: epdp modules
import model
import _gfrd
import gfrdbase
#XXX


def add_species_type(world, st):
    try:
        structure = st["structure"]
    except _gfrd.NotFound:
        structure = "world"

    world.add_species(
        _gfrd.SpeciesInfo(
            st.id,  float(st["D"]), float(st["radius"]),
            structure, float(st["v"])))

def create_world(m, world_size, matrix_size):
    x = numpy.repeat(world_size * 0.5, 3)
    world_region = _gfrd.CuboidalRegion('world', _gfrd.Box(x, x))
    # m.set_all_repulsive()
    # world_region = m.get_structure("world")
    # if not isinstance(world_region, _gfrd.CuboidalRegion):
    #     raise TypeError("the world should be a CuboidalRegion")

    if not numpy.all(
        world_region.shape.half_extent == world_region.shape.half_extent[0]):
        raise NotImplementedError("non-cuboidal world is not supported")

    # world_size = world_region.shape.half_extent[0] * 2

    world = _gfrd.World(world_size, matrix_size)

    for st in m.species_types:
        add_species_type(world, st)

    for r in m.structures.itervalues():
        world.add_structure(r)

    # world.model = m
    return world

class World:

    def __init__(self, world_size, matrix_size):
        self.world_size = world_size
        self.matrix_size = matrix_size

        # self.model.network_rules
        # self.model.get_species_type_by_id(sid)
        # self.structures

        m = model.ParticleModel(self.world_size)
        w = create_world(m, self.world_size, self.matrix_size)

        self.model = m
        self.world = w
        self.__species_id_map = {}

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
        if self.ecell4__has_species(sp):
            raise RuntimeError, "already exists [%s]." % (sp.serial())

        st = model.Species(
            sp.serial(), sp.get_attribute("D"), sp.get_attribute("radius"))
        self.model.add_species_type(st)
        add_species_type(self.world, st)
        self.__species_id_map[sp.serial()] = st.id

    def ecell4__has_species(self, sp):
        return (sp.serial() in self.__species_id_map.keys())

    def ecell4__add_molecules(self, sp, num):
        if not self.ecell4__has_species(sp):
            self.ecell4__add_species(sp)
        sid = self.__species_id_map[sp.serial()]
        gfrdbase.throw_in_particles(self.world, sid, num)

    def ecell4__num_particles(self, sp=None):
        if sp is None:
            return self.world.num_particles
        elif self.ecell4__has_species(sp):
            sid = self.__species_id_map[sp.serial()]
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
