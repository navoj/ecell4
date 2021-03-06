from libcpp.string cimport string
from libcpp cimport bool

from core cimport *


## Cpp_EGFRDWorld
#  ecell4::egfrd::EGFRDWorld
cdef extern from "ecell4/egfrd/egfrd.hpp" namespace "ecell4::egfrd":
    cdef cppclass Cpp_EGFRDWorld "ecell4::egfrd::EGFRDWorld":
        Cpp_EGFRDWorld() except +
        Cpp_EGFRDWorld(Cpp_Real3&) except +
        Cpp_EGFRDWorld(Cpp_Real3&, Cpp_Integer3&) except +
        Cpp_EGFRDWorld(
            Cpp_Real3&, Cpp_Integer3&,
            shared_ptr[Cpp_RandomNumberGenerator]&) except +
        #     shared_ptr[Cpp_GSLRandomNumberGenerator]&) except +
        Cpp_EGFRDWorld(string&) except +
        pair[pair[Cpp_ParticleID, Cpp_Particle], bool] new_particle(Cpp_Particle& p)
        pair[pair[Cpp_ParticleID, Cpp_Particle], bool] new_particle(Cpp_Species& sp, Cpp_Real3& pos)
        void set_t(Real t)
        Real t()
        Cpp_Real3 edge_lengths()
        Integer num_particles()
        Integer num_particles(Cpp_Species& sp)
        Integer num_particles_exact(Cpp_Species& sp)
        vector[pair[Cpp_ParticleID, Cpp_Particle]] list_particles()
        vector[pair[Cpp_ParticleID, Cpp_Particle]] list_particles(Cpp_Species& sp)
        vector[pair[Cpp_ParticleID, Cpp_Particle]] list_particles_exact(Cpp_Species& sp)
        bool has_particle(Cpp_ParticleID& pid)
        bool update_particle(Cpp_ParticleID& pid, Cpp_Particle& p)
        pair[Cpp_ParticleID, Cpp_Particle] get_particle(Cpp_ParticleID& pid)
        void remove_particle(Cpp_ParticleID& pid)
        vector[pair[pair[Cpp_ParticleID, Cpp_Particle], Real]] list_particles_within_radius(Cpp_Real3& pos, Real& radius)
        vector[pair[pair[Cpp_ParticleID, Cpp_Particle], Real]] list_particles_within_radius(Cpp_Real3& pos, Real& radius, Cpp_ParticleID& ignore)
        vector[pair[pair[Cpp_ParticleID, Cpp_Particle], Real]] list_particles_within_radius(Cpp_Real3& pos, Real& radius, Cpp_ParticleID& ignore1, Cpp_ParticleID& ignore2)
        # Cpp_Real3 periodic_transpose(Cpp_Real3& pos1, Cpp_Real3& pos2)
        Cpp_Real3 apply_boundary(Cpp_Real3& pos)
        # Real distance_sq(Cpp_Real3& pos1, Cpp_Real3& pos2)
        Real distance(Cpp_Real3& pos1, Cpp_Real3& pos2)
        Real volume()
        bool has_species(Cpp_Species& sp)
        Integer num_molecules(Cpp_Species& sp)
        Integer num_molecules_exact(Cpp_Species& sp)
        void add_molecules(Cpp_Species& sp, Integer num)
        void add_molecules(Cpp_Species& sp, Integer num, shared_ptr[Cpp_Shape])
        void remove_molecules(Cpp_Species& sp, Integer num)
        void save(string filename) except +
        void load(string filename)
        void bind_to(shared_ptr[Cpp_Model])
        shared_ptr[Cpp_RandomNumberGenerator] rng()

    cdef cppclass Cpp_EGFRDSimulator "ecell4::egfrd::EGFRDSimulator":
        #XXX: be carefull about the order of arguments
        Cpp_EGFRDSimulator(
            shared_ptr[Cpp_EGFRDWorld]&, shared_ptr[Cpp_Model]&) except +
        # Cpp_EGFRDSimulator(shared_ptr[Cpp_EGFRDWorld]&) except +
        Integer num_steps()
        void step()
        bool step(Real)
        Real t()
        void set_t(Real)
        void set_dt(Real)
        Real dt()
        Real next_time()
        vector[Cpp_ReactionRule] last_reactions()
        void initialize()
        # Cpp_GSLRandomNumberGenerator& rng()
        shared_ptr[Cpp_Model] model()
        shared_ptr[Cpp_EGFRDWorld] world()
        void run(Real)
        void run(Real, shared_ptr[Cpp_Observer])
        void run(Real, vector[shared_ptr[Cpp_Observer]])

    cdef cppclass Cpp_EGFRDFactory "ecell4::egfrd::EGFRDFactory":
        Cpp_EGFRDFactory() except +
        Cpp_EGFRDFactory(Integer) except +
        Cpp_EGFRDFactory(Integer, Real) except +
        Cpp_EGFRDFactory(Integer, Real, Real) except +
        Cpp_EGFRDFactory(Cpp_Integer3&) except +
        Cpp_EGFRDFactory(Cpp_Integer3&, Integer) except +
        Cpp_EGFRDFactory(Cpp_Integer3&, Integer, Real) except +
        Cpp_EGFRDFactory(Cpp_Integer3&, Integer, Real, Real) except +
        Cpp_EGFRDFactory(Cpp_Integer3&, shared_ptr[Cpp_RandomNumberGenerator]&) except +
        Cpp_EGFRDFactory(Cpp_Integer3&, shared_ptr[Cpp_RandomNumberGenerator]&, Integer) except +
        Cpp_EGFRDFactory(Cpp_Integer3&, shared_ptr[Cpp_RandomNumberGenerator]&, Integer, Real) except +
        Cpp_EGFRDFactory(Cpp_Integer3&, shared_ptr[Cpp_RandomNumberGenerator]&, Integer, Real, Real) except +
        Cpp_EGFRDWorld* create_world(string)
        Cpp_EGFRDWorld* create_world(Cpp_Real3&)
        Cpp_EGFRDWorld* create_world(shared_ptr[Cpp_Model])
        Cpp_EGFRDSimulator* create_simulator(shared_ptr[Cpp_Model], shared_ptr[Cpp_EGFRDWorld])
        Cpp_EGFRDSimulator* create_simulator(shared_ptr[Cpp_EGFRDWorld])

cdef class EGFRDWorld:
    cdef shared_ptr[Cpp_EGFRDWorld]* thisptr

cdef class EGFRDSimulator:
    cdef Cpp_EGFRDSimulator* thisptr

cdef class EGFRDFactory:
    cdef Cpp_EGFRDFactory* thisptr

cdef EGFRDWorld EGFRDWorld_from_Cpp_EGFRDWorld(
    shared_ptr[Cpp_EGFRDWorld] m)

cdef EGFRDSimulator EGFRDSimulator_from_Cpp_EGFRDSimulator(
    Cpp_EGFRDSimulator* s)
