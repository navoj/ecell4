#include <stdexcept>
#include <fstream>

#include "LatticeWorld.hpp"

namespace ecell4
{

namespace lattice
{

LatticeWorld* create_lattice_world_cell_list_impl(
    const Real3& edge_lengths, const Real& voxel_radius,
    const Integer3& matrix_sizes,
    const boost::shared_ptr<RandomNumberGenerator>& rng)
{
    return new LatticeWorld(
        new LatticeSpaceCellListImpl(edge_lengths, voxel_radius, matrix_sizes), rng);
}

LatticeWorld* create_lattice_world_vector_impl(
    const Real3& edge_lengths, const Real& voxel_radius,
    const boost::shared_ptr<RandomNumberGenerator>& rng)
{
    return new LatticeWorld(
        new LatticeSpaceVectorImpl(edge_lengths, voxel_radius), rng);
}

const Real& LatticeWorld::t() const
{
    return (*space_).t();
}

void LatticeWorld::set_t(const Real& t)
{
    (*space_).set_t(t);
}

const Real3& LatticeWorld::edge_lengths() const
{
    return (*space_).edge_lengths();
}

const Real LatticeWorld::volume() const
{
    return (*space_).volume();
}

Integer LatticeWorld::num_species() const
{
    return (*space_).num_species();
}

bool LatticeWorld::has_species(const Species &sp) const
{
    return (*space_).has_species(sp);
}

// bool LatticeWorld::has_species_exact(const Species &sp) const
// {
//     return (*space_).has_species_exact(sp);
// }

Integer LatticeWorld::num_molecules(const Species& sp) const
{
    return (*space_).num_molecules(sp);
}

Integer LatticeWorld::num_molecules_exact(const Species& sp) const
{
    return (*space_).num_molecules_exact(sp);
}

Integer LatticeWorld::num_molecules() const
{
    return (*space_).num_molecules();
}

Integer LatticeWorld::num_particles(const Species& sp) const
{
    return (*space_).num_particles(sp);
}

Integer LatticeWorld::num_particles_exact(const Species& sp) const
{
    return (*space_).num_particles_exact(sp);
}

Integer LatticeWorld::num_particles() const
{
    return (*space_).num_particles();
}

Integer LatticeWorld::num_voxels() const
{
    return (*space_).num_voxels();
}

Integer LatticeWorld::num_voxels(const Species& sp) const
{
    return (*space_).num_voxels(sp);
}

Integer LatticeWorld::num_voxels_exact(const Species& sp) const
{
    return (*space_).num_voxels_exact(sp);
}

bool LatticeWorld::has_particle(const ParticleID& pid) const
{
    return (*space_).has_particle(pid);
}

bool LatticeWorld::has_voxel(const ParticleID& pid) const
{
    return (*space_).has_voxel(pid);
}

std::vector<std::pair<ParticleID, Particle> >
LatticeWorld::list_particles() const
{
    return (*space_).list_particles();
}

std::vector<std::pair<ParticleID, Particle> >
LatticeWorld::list_particles(const Species& sp) const
{
    return (*space_).list_particles(sp);
}

std::vector<std::pair<ParticleID, Particle> >
LatticeWorld::list_particles_exact(const Species& sp) const
{
    return (*space_).list_particles_exact(sp);
}

bool LatticeWorld::update_particle(const ParticleID& pid, const Particle& p)
{
    return (*space_).update_particle(pid, p);
}

std::vector<Species> LatticeWorld::list_species() const
{
    return (*space_).list_species();
}

std::vector<std::pair<ParticleID, Voxel> > LatticeWorld::list_voxels() const
{
    return (*space_).list_voxels();
}

std::vector<std::pair<ParticleID, Voxel> >
    LatticeWorld::list_voxels(const Species& sp) const
{
    return (*space_).list_voxels(sp);
}

std::vector<std::pair<ParticleID, Voxel> >
    LatticeWorld::list_voxels_exact(const Species& sp) const
{
    return (*space_).list_voxels_exact(sp);
}

MolecularTypeBase* LatticeWorld::find_molecular_type(const Species& species)
{
    return (*space_).find_molecular_type(species);
}

MolecularTypeBase* LatticeWorld::get_molecular_type_private(
        const private_coordinate_type& coord)
{
    return (*space_).get_molecular_type(coord);
}

std::pair<std::pair<ParticleID, Voxel>, bool>
LatticeWorld::new_voxel(const Voxel& v)
{
    const private_coordinate_type private_coord(coord2private(v.coordinate()));
    return new_voxel_private(
        Voxel(v.species(), private_coord, v.radius(), v.D(), v.loc()));
}

std::pair<std::pair<ParticleID, Voxel>, bool>
LatticeWorld::new_voxel(const Species& sp, const coordinate_type& coord)
{
    const private_coordinate_type private_coord(coord2private(coord));
    const molecule_info_type minfo(get_molecule_info(sp));
    return new_voxel_private(
        Voxel(sp, private_coord, minfo.radius, minfo.D, minfo.loc));
}

std::pair<std::pair<ParticleID, Voxel>, bool>
LatticeWorld::new_voxel_private(const Voxel& v)
{
    ParticleID pid(sidgen_());
    const bool is_succeeded((*space_).update_voxel_private(pid, v));
    const coordinate_type coord(private2coord(v.coordinate()));
    return std::make_pair(
        std::make_pair(pid, Voxel(v.species(), coord, v.radius(), v.D(), v.loc())),
        is_succeeded);
}

std::pair<std::pair<ParticleID, Voxel>, bool>
LatticeWorld::new_voxel_private(
        const Species& sp, const private_coordinate_type& coord)
{
    const molecule_info_type minfo(get_molecule_info(sp));
    return new_voxel_private(
        Voxel(sp, coord, minfo.radius, minfo.D, minfo.loc));
}

std::pair<std::pair<ParticleID, Voxel>, bool>
LatticeWorld::new_voxel_structure(const Voxel& v)
{
    const bool is_succeeded((*space_).update_voxel_private(ParticleID(), v));
    const coordinate_type coord(private2coord(v.coordinate()));
    return std::make_pair(std::make_pair(ParticleID(),
                Voxel(v.species(), coord, v.radius(), v.D(), v.loc())),
        is_succeeded);
}

bool LatticeWorld::add_molecules(const Species& sp, const Integer& num)
{
    if (num < 0)
    {
        throw std::invalid_argument("The number of molecules must be positive.");
    }

    std::cerr << "[DEBUG]" << std::endl; // DEBUG
    const LatticeWorld::molecule_info_type info(get_molecule_info(sp));
    std::cerr << "  info.D = " << info.D << std::endl; // DEBUG

    Integer count(0);
    while (count < num)
    {
        const coordinate_type coord(rng()->uniform_int(0, (*space_).size() - 1));
        const Voxel v(sp, coord2private(coord), info.radius, info.D, info.loc);

        if ((*space_).on_structure(v))
        {
            continue;
        }
        else if (new_voxel_private(v).second)
        {
            ++count;
        }
    }
    return true;
}

bool LatticeWorld::add_molecules(
    const Species& sp, const Integer& num, const boost::shared_ptr<const Shape> shape)
{
    if (num < 0)
    {
        throw std::invalid_argument("The number of molecules must be positive.");
    }

    const LatticeWorld::molecule_info_type info(get_molecule_info(sp));

    Integer count(0);
    while (count < num)
    {
        const Real3 pos(shape->draw_position(rng_));
        const Voxel v(sp, (*space_).position2private(pos), info.radius, info.D, info.loc);

        if ((*space_).on_structure(v))
        {
            continue;
        }
        else if (new_voxel_private(v).second)
        {
            ++count;
        }
    }
    return true;
}

Integer LatticeWorld::add_structure(
    const Species& sp, const boost::shared_ptr<const Shape> shape)
{
    const LatticeWorld::molecule_info_type info(get_molecule_info(sp));
    (*space_).add_structure(sp, shape, info.loc);

    switch (shape->dimension())
    {
    case Shape::THREE:
        return add_structure3(sp, shape);
    case Shape::TWO:
        return add_structure2(sp, shape);
    }

    throw NotSupported("The dimension of a shape must be two or three.");
}

Integer LatticeWorld::add_structure3(const Species& sp, const boost::shared_ptr<const Shape> shape)
{
    // Real3 l, u;
    // shape->bounding_box(edge_lengths(), l, u);
    // const Real sigma(voxel_radius() * 2);
    // const unsigned int ndim(3);
    // for (unsigned int i(0); i != ndim; ++i)
    // {
    //     l[i] = std::max(0.0, l[i] - sigma);
    //     u[i] = std::min(edge_lengths()[i], u[i] + sigma);
    // }
    // const Integer3 lower(position2global(l)), upper(position2global(u));

    const LatticeWorld::molecule_info_type info(get_molecule_info(sp));
    Integer count(0);
    for (Integer col(0); col < col_size(); ++col)
    {
        for (Integer row(0); row < row_size(); ++row)
        {
            for (Integer layer(0); layer < layer_size(); ++layer)
            {
    // for (Integer col(lower[0]); col < upper[0]; ++col)
    // {
    //     for (Integer row(lower[1]); row < upper[1]; ++row)
    //     {
    //         for (Integer layer(lower[2]); layer < upper[2]; ++layer)
    //         {
                const Integer3 g(col, row, layer);
                const Real L(shape->is_inside(global2position(g)));
                if (L > 0)
                {
                    continue;
                }

                const Voxel v(sp, (*space_).global2private_coord(g),
                    info.radius, info.D, info.loc);
                if (new_voxel_structure(v).second)
                {
                    ++count;
                }
            }
        }
    }
    return count;
}

Integer LatticeWorld::add_structure2(const Species& sp, const boost::shared_ptr<const Shape> shape)
{
    // std::ofstream fout("shape.csv");
    // fout << "# " << sp.serial() << std::endl;
    // const unsigned int n(50);
    // const Real3 L(edge_lengths() / static_cast<Real>(n));
    // for (unsigned int i(0); i < n * n * n; ++i)
    // {
    //     const unsigned int x(i % n);
    //     const unsigned int y((i / n) % n);
    //     const unsigned int z(i / (n * n));
    //     if (shape->test_AABB(Real3(x * L[0], y * L[1], z * L[2]),
    //         Real3((x + 1) * L[0], (y + 1) * L[1], (z + 1) * L[2])))
    //     {
    //         fout << x << "," << y << "," << z << std::endl;
    //     }
    // }
    // fout.close();

    // Real3 l, u;
    // shape->bounding_box(edge_lengths(), l, u);
    // const Real sigma(voxel_radius() * 2);
    // const unsigned int ndim(3);
    // for (unsigned int i(0); i != ndim; ++i)
    // {
    //     l[i] = std::max(0.0, l[i] - sigma);
    //     u[i] = std::min(edge_lengths()[i], u[i] + sigma);
    // }
    // const Integer3 lower(position2global(l)), upper(position2global(u));

    const LatticeWorld::molecule_info_type info(get_molecule_info(sp));
    Integer count(0);
    for (Integer col(0); col < col_size(); ++col)
    {
        for (Integer row(0); row < row_size(); ++row)
        {
            for (Integer layer(0); layer < layer_size(); ++layer)
            {
    // for (Integer col(lower[0]); col < upper[0]; ++col)
    // {
    //     for (Integer row(lower[1]); row < upper[1]; ++row)
    //     {
    //         for (Integer layer(lower[2]); layer < upper[2]; ++layer)
    //         {
                const Integer3 g(col, row, layer);
                if (!is_surface_voxel(g, shape))
                {
                    continue;
                }

                const Voxel v(sp, (*space_).global2private_coord(g),
                    info.radius, info.D, info.loc);
                if (new_voxel_structure(v).second)
                {
                    ++count;
                }
            }
        }
    }
    return count;
}

bool LatticeWorld::is_surface_voxel(
    const Integer3& g, const boost::shared_ptr<const Shape> shape) const
{
    const Real L(shape->is_inside(global2position(g)));
    if (L > 0 || L < -2 * voxel_radius())
    {
        return false;
    }

    const LatticeWorld::private_coordinate_type
        private_coord((*space_).global2private_coord(g));
    for (Integer i(0); i < 12; ++i)
    {
        if (shape->is_inside(global2position((*space_).private_coord2global(
            (*space_).get_neighbor(private_coord, i)))) > 0)
        {
            return true;
        }
    }
    return false;
}

// TODO
Integer LatticeWorld::add_neighbors(const Species& sp,
    const LatticeWorld::private_coordinate_type center)
{
    Integer count(0);
    const LatticeWorld::molecule_info_type info(get_molecule_info(sp));
    for (Integer i(0); i < 12; ++i)
    {
        const private_coordinate_type n((*space_).get_neighbor(center, i));
        if (new_voxel_private(Voxel(sp, n, info.radius, info.D, info.loc)).second)
        {
            ++count;
        }
        else
        {
            throw "Error in add_neighbors()";
        }
    }
    return count;

    // Integer count(0);
    // const LatticeWorld::molecule_info_type info(get_molecule_info(sp));
    // std::vector<LatticeWorld::private_coordinate_type> neighbors(
    //         (*space_).get_neighbors(center));
    // for (std::vector<LatticeWorld::private_coordinate_type>::iterator itr(
    //             neighbors.begin()); itr != neighbors.end(); itr++)
    // {
    //     if (new_voxel_private(Voxel(sp, *itr, info.radius, info.D, info.loc)).second)
    //         ++count;
    //     else
    //         throw "Error in add_neighbors()";
    // }
    // return count;
}
// TODO

void LatticeWorld::remove_molecules(const Species& sp, const Integer& num)
{
    if (num < 0)
    {
        throw std::invalid_argument("The number of molecules must be positive.");
    }

    MolecularTypeBase* mtype(find_molecular_type(sp));
    if (!mtype->with_voxels())
    {
        throw NotSupported(
            "remove_molecuels for MolecularType with no voxel is not supported now");
    }
    else if (mtype->size() < num)
    {
        throw std::invalid_argument(
            "The number of molecules cannot be negative.");
    }

    Integer count(0);
    while (count < num)
    {
        const Integer idx(rng_->uniform_int(0, mtype->size() - 1));
        if (remove_voxel_private(mtype->at(idx).first))
        {
            ++count;
        }
    }
}

bool LatticeWorld::remove_voxel_private(const private_coordinate_type coord)
{
    return (*space_).remove_voxel_private(coord);
}

bool LatticeWorld::move(coordinate_type from, coordinate_type to)
{
    return (*space_).move(from, to);
}

std::pair<LatticeWorld::private_coordinate_type, bool>
LatticeWorld::move_to_neighbor(
    MolecularTypeBase* const& from_mt, MolecularTypeBase* const& loc,
    particle_info_type& info, const Integer nrand)
{
    return (*space_).move_to_neighbor(from_mt, loc, info, nrand);
}

std::pair<LatticeWorld::private_coordinate_type, bool>
LatticeWorld::check_neighbor_private(
    const private_coordinate_type coord, const std::string& loc)
{
    std::vector<private_coordinate_type> tmp;
    tmp.reserve(12);
    for (unsigned int rnd(0); rnd < 12; ++rnd)
    {
        const private_coordinate_type
            neighbor((*space_).get_neighbor(coord, rnd));
        const MolecularTypeBase* mt(get_molecular_type_private(neighbor));
        const std::string
            serial(mt->is_vacant() ? "" : mt->species().serial());
        if (serial == loc)
        {
            tmp.push_back(neighbor);
        }
    }

    if (tmp.size() == 0)
    {
        return std::make_pair(coord, false);
    }

    return std::make_pair(
        tmp[rng()->uniform_int(0, tmp.size() - 1)], true);

    // const Integer rnd(rng()->uniform_int(0, 11));
    // const private_coordinate_type neighbor((*space_).get_neighbor(coord, rnd));
    // bool flg = get_molecular_type_private(neighbor)->is_vacant(); //XXX: loc
    // return std::make_pair(neighbor, flg);
}

Shape::dimension_kind LatticeWorld::get_dimension_kind(const std::string& name) const
{
    if (name == "")
    {
        return Shape::THREE; // Default value for VACANT type.
    }
    return (*space_).get_structure_dimension(Species(name));
}

} // lattice

} // ecell4
