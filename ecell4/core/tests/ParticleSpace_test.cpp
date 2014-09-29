#define BOOST_TEST_MODULE "ParticleSpace_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include <boost/test/test_case_template.hpp>

#include <ecell4/core/types.hpp>
#include <ecell4/core/SerialIDGenerator.hpp>
#include <ecell4/core/ParticleSpace.hpp>
#include <ecell4/core/ParticleSpaceCellListImpl.hpp>

using namespace ecell4;

template<typename Timpl_>
void ParticleSpace_test_edge_lengths_template()
{
    const Real L(1e-6);
    Timpl_ target(Position3(L, L, L));

    const Position3 edge_lengths(target.edge_lengths());
    for (Position3::size_type dim(0); dim < 3; ++dim)
    {
        BOOST_CHECK_CLOSE(edge_lengths[dim], L, 1e-6);
    }
}

BOOST_AUTO_TEST_CASE(ParticleSpace_test_edge_lengths)
{
    ParticleSpace_test_edge_lengths_template<ParticleSpaceVectorImpl>();
    ParticleSpace_test_edge_lengths_template<ParticleSpaceCellListImpl>();
}

template<typename Timpl_>
void ParticleSpace_test_manipulate_particle_template()
{
    const Real L(1e-6);
    const Real L_2(L * 0.5);
    Timpl_ target(Position3(L, L, L));
    SerialIDGenerator<ParticleID> pidgen;

    const Species sp1("A"), sp2("B");
    const Real radius(2.5e-9), D(1e-12);
    const ParticleID pid(pidgen());
    BOOST_CHECK(
        target.update_particle(pid, Particle(sp1, Position3(0, 0, 0), radius, D)));
    BOOST_CHECK(target.has_particle(pid));
    BOOST_CHECK_EQUAL(target.num_particles(), 1);
    BOOST_CHECK_EQUAL(target.num_particles(sp1), 1);
    BOOST_CHECK_EQUAL(target.num_particles_exact(sp1), 1);
    BOOST_CHECK_EQUAL(target.num_particles(sp2), 0);
    BOOST_CHECK_EQUAL(target.num_particles_exact(sp2), 0);

    const std::pair<ParticleID, Particle>
        pid_particle_pair1(target.get_particle(pid));
    BOOST_CHECK_EQUAL(pid_particle_pair1.first, pid);
    BOOST_CHECK(pid_particle_pair1.second.species() == sp1);
    BOOST_CHECK_CLOSE(pid_particle_pair1.second.radius(), radius, 1e-6);
    BOOST_CHECK_CLOSE(pid_particle_pair1.second.D(), D, 1e-6);

    for (Position3::size_type dim(0); dim < 3; ++dim)
    {
        BOOST_CHECK_CLOSE(pid_particle_pair1.second.position()[dim], 0, 1e-6);
    }

    BOOST_CHECK(
        !target.update_particle(pid, Particle(sp2, Position3(L_2, L_2, L_2), radius, D)));
    BOOST_CHECK(target.has_particle(pid));
    BOOST_CHECK_EQUAL(target.num_particles(), 1);
    BOOST_CHECK_EQUAL(target.num_particles(sp1), 0);
    BOOST_CHECK_EQUAL(target.num_particles_exact(sp1), 0);
    BOOST_CHECK_EQUAL(target.num_particles(sp2), 1);
    BOOST_CHECK_EQUAL(target.num_particles_exact(sp2), 1);

    const std::pair<ParticleID, Particle>
        pid_particle_pair2(target.get_particle(pid));
    BOOST_CHECK_EQUAL(pid_particle_pair2.first, pid);
    BOOST_CHECK(pid_particle_pair2.second.species() == sp2);
    for (Position3::size_type dim(0); dim < 3; ++dim)
    {
        BOOST_CHECK_CLOSE(pid_particle_pair2.second.position()[dim], L_2, 1e-6);
    }

    target.remove_particle(pid);
    BOOST_CHECK(!target.has_particle(pid));
    BOOST_CHECK_EQUAL(target.num_particles(), 0);
    BOOST_CHECK_EQUAL(target.num_particles(sp2), 0);
    BOOST_CHECK_EQUAL(target.num_particles_exact(sp2), 0);
}

BOOST_AUTO_TEST_CASE(ParticleSpace_test_manupulate_particle)
{
    ParticleSpace_test_manipulate_particle_template<ParticleSpaceVectorImpl>();
    ParticleSpace_test_manipulate_particle_template<ParticleSpaceCellListImpl>();
}
