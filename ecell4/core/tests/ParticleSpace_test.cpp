#define BOOST_TEST_MODULE "ParticleSpace_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include <boost/test/test_case_template.hpp>

#include <ecell4/core/types.hpp>
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
