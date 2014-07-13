#ifndef ECELL4_EGFRD_WORLD_HPP
#define ECELL4_EGFRD_WORLD_HPP

#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/Species.hpp>
#include <ecell4/core/types.hpp>
#include <ecell4/core/Position3.hpp>
#include <ecell4/core/ParticleSpace.hpp>
#include <ecell4/core/functions.hpp>
//#include <ecell4/core/SerialIDGenerator.hpp>
#include "SerialIDGenerator.hpp"
#include "ParticleTraits.hpp"

#include <map>
#include <set>
#include <boost/lexical_cast.hpp>
#include "Structure.hpp"
#include "exceptions.hpp"

namespace ecell4 {

namespace egfrd {

// Following struct definition is to solve the 
// typedef definitions of World::ParticleContainerType, 
// but it is needed to adapt BDPropagator etc.
template <typename Ttraits_>
struct ParticleContainerCorrespondence 
{
    typedef Ttraits_ traits_type;
    typedef typename traits_type::length_type length_type;
    typedef typename traits_type::D_type D_type;
    typedef typename traits_type::molecule_info molecule_info;  //XXX
    typedef typename traits_type::position_type position_type;
    typedef typename traits_type::size_type size_type;
    typedef typename traits_type::species_type species_type;
    typedef typename traits_type::species_id_type species_id_type;
    typedef typename traits_type::structure_id_type structure_id_type;
    typedef typename traits_type::structure_type structure_type;
    typedef typename traits_type::particle_type particle_type;
    typedef typename traits_type::particle_id_type particle_id_type;
    typedef typename traits_type::particle_id_generator particle_id_generator;
    typedef typename traits_type::particle_shape_type particle_shape_type;
    typedef typename std::pair</*const*/ particle_id_type, particle_type> particle_id_pair;

    typedef sized_iterator_range<particle_id_pair> particle_id_pair_range;

    typedef std::pair<particle_id_pair, length_type> particle_id_pair_and_distance;
    typedef std::vector<particle_id_pair_and_distance> particle_id_pair_and_distance_list;
};

template <typename Ttraits_>
class EGFRDWorld
{
public:
    typedef Ttraits_ traits_type;
    typedef typename traits_type::length_type length_type;
    typedef typename traits_type::D_type D_type;
    typedef typename traits_type::molecule_info molecule_info;  //XXX
    typedef typename traits_type::position_type position_type;
    typedef typename traits_type::size_type size_type;
    typedef typename traits_type::species_type species_type;
    typedef typename traits_type::species_id_type species_id_type;
    typedef typename traits_type::structure_id_type structure_id_type;
    typedef typename traits_type::structure_type structure_type;
    typedef typename traits_type::particle_type particle_type;
    typedef typename traits_type::particle_id_type particle_id_type;
    typedef typename traits_type::particle_id_generator particle_id_generator;
    typedef typename traits_type::particle_shape_type particle_shape_type;
    typedef typename std::pair</*const*/particle_id_type, particle_type> particle_id_pair;
    typedef ParticleContainerCorrespondence<traits_type> particle_container_type;

protected:
    typedef std::map<species_id_type, molecule_info> species_map;
    typedef std::map<structure_id_type, boost::shared_ptr<structure_type> > structure_map;
    typedef std::set<particle_id_type> particle_id_set;
    typedef std::map<species_id_type, particle_id_set> per_species_particle_id_set;
    typedef select_second<typename species_map::value_type> species_second_selector_type;
    typedef select_second<typename structure_map::value_type> surface_second_selector_type;

public:
    typedef boost::transform_iterator<species_second_selector_type,
            typename species_map::const_iterator> species_iterator;
    typedef boost::transform_iterator<surface_second_selector_type,
            typename structure_map::const_iterator> surface_iterator;
    typedef sized_iterator_range<species_iterator> species_range;
    typedef sized_iterator_range<surface_iterator> structures_range;
    typedef sized_iterator_range<particle_id_pair> particle_id_pair_range;

    typedef std::pair<particle_id_pair, length_type> particle_id_pair_and_distance;
    typedef std::vector<particle_id_pair_and_distance> particle_id_pair_and_distance_list;

    struct particle_id_pair_and_distance_comparator {
        typedef particle_id_pair_and_distance arg_type;
        bool operator()(arg_type const &lhs, arg_type const &rhs) const
        {
            return lhs.second < rhs.second;
        }
    };
public:
    EGFRDWorld(length_type world_size = 1, size_type size = 1) : 
        world_size_(world_size), size_(size), 
        ps_(new ecell4::ParticleSpaceVectorImpl(
                    position_type(world_size, world_size, world_size)) )
    {}

    particle_id_pair new_particle(species_id_type const &sid, 
            position_type const &pos)
    {
        molecule_info const &info(get_species(sid));
        particle_id_pair retval(
                pidgen_(), 
                particle_type(sid, pos, info.radius(), info.D() ));
        update_particle(retval);
        return retval;
    }

    bool update_particle(particle_id_pair const &pi_pair)
    {
        return ps_->update_particle(pi_pair.first, pi_pair.second);
    }
    bool update_particle(const particle_id_type &pid, const particle_type &p)
    {
        return ps_->update_particle(pid, p);
    }

    bool remove_particle(particle_id_type const &id)
    {
        bool found(false);
        particle_id_pair pp(ps_->get_particle(id));
        if (!found)
        {
            return false;
        }
        this->particle_pool_[pp.second.sid()].erase(id);
        ps_->remove_particle(id);
        return true;
    }

    bool has_particle(particle_id_type const &pid) const
    {
        return ps_->has_particle(pid);
    }

    particle_id_set get_particle_ids(species_id_type const &sid) const
    {
        typename per_species_particle_id_set::const_iterator i(
                particle_pool_.find(sid));
        if (i == this->particle_pool_.end())
        {
            throw not_found("Unknown species");
        }
        return i->second;
    }

    particle_id_pair_range get_particles_range() const
    {
        return particle_id_pair_range(
                ps_->list_particles().begin(), 
                ps_->list_particles().end(), 
                ps_->list_particles().size());
    }

    size_type num_molecules(const ecell4::Species &sp) const
    {
        return ps_->num_particles(sp);
    }
    size_type num_molecules(const species_id_type &sid) const
    {
        return ps_->num_particles( ecell4::Species(sid) );
    }
    size_type num_particles() const
    {
        return ps_->num_particles();
    }
    size_type num_particles(const species_id_type &sid) const 
    {
        return ps_->num_particles(ecell4::Species(sid));
    }

    bool add_structure(boost::shared_ptr<structure_type> surface)
    {
        return structure_map_.insert(std::make_pair(surface->id(), surface)).second;
    }

    boost::shared_ptr<structure_type> get_structure(structure_id_type const &id) const
    {
        typename structure_map::const_iterator i(this->structure_map_.find(id));
        if (structure_map_.end() == i)
        {
            throw not_found(std::string("Unknown surface(id = ") + boost::lexical_cast<std::string>(id) + ")");
        }
        return i->second;
    }

    molecule_info get_molecule_info(const ecell4::Species &sp) const
    {
        const ecell4::Species::serial_type sid(sp.serial());
        const Real D(std::atof( sp.get_attribute("D").c_str() ));
        const Real radius(std::atof( sp.get_attribute("radius").c_str() ));
        const structure_id_type structure_id( sp.has_attribute("structure_id") ? 
                        sp.get_attribute("structure_id") : 
                        structure_id_type("world") );
        return molecule_info(sid, D, radius, structure_id);
    }

    void add_species(molecule_info const &info)
    {
        this->species_map_[info.id()] = info;
        this->particle_pool_[info.id()] = particle_id_set();
    }
    
    const molecule_info &get_species(species_id_type const &sid) const
    {
        typename species_map::const_iterator i(species_map_.find(sid));
        if (species_map_.end() == i)
        {
            throw not_found(std::string("Unknown Species"));
        }
        return i->second;
    }
    species_range get_species() const 
    {
        return species_range(
                species_iterator(species_map_.begin(), species_second_selector_type()),
                species_iterator(species_map_.end(), species_second_selector_type()),
                species_map_.size());
    }

    size_type world_size() const
    {
        return this->world_size_;
    }
    size_type cell_size() const
    {
        return world_size_ / size_;
    }
    size_type matrix_size() const
    {
        return 1;
    }

    particle_id_pair_and_distance_list list_particles_within_radius(
            const position_type &pos, const Real& radius) const
    {
        return ps_->list_particles_within_radius(pos, radius);
    }

    particle_id_pair_and_distance_list list_particles_within_radius(
            const position_type &pos, const Real& radius, 
            const particle_id_type &ignore) const
    {
        return ps_->list_particles_within_radius(pos, radius, ignore);
    }

    particle_id_pair_and_distance_list list_particles_within_radius(
            const position_type &pos, const Real& radius,
            const particle_id_type &ignore1, const particle_id_type &ignore2)
    {
        return ps_->list_particles_within_radius(pos, radius, ignore1, ignore2);
    }

    particle_id_pair_and_distance_list list_particles_within_radius(
            particle_shape_type const &s) const 
    {
        return ps_->list_particles_within_radius(s.position(), s.radius());
    }
    particle_id_pair_and_distance_list list_particles_within_radius(
            particle_shape_type const &s, particle_id_type const &ignore) const
    {
        return ps_->list_particles_within_radius(s.position(), s.radius(), ignore);
    }

    particle_id_pair_and_distance_list list_particles_within_radius(
            particle_shape_type const &s, 
            particle_id_type const &ignore1, particle_id_type const &ignore2) const
    {
        return ps_->list_particles_within_radius(s.position(), s.radius(), ignore1, ignore2);
    }

    template <typename T>
    length_type distance(const T &s1, const position_type &pos2) const
    {
        return ps_->distance(s1.position(), pos2);
    }

    length_type distance(const position_type &pos1, const position_type &pos2) const
    {
        return ps_->distance(pos1, pos2);
    }

    //
    position_type apply_boundary(position_type const &pos) const
    {
        // XXX  
        return ps_->apply_boundary(pos);
    }

    length_type apply_boundary(length_type const &v) const
    {
        return ecell4::modulo(v, world_size() );
    }

    position_type cyclic_transpose(position_type const &p0, position_type const &p1) const
    {   
        // XXX  
        return ps_->periodic_transpose(p0, p1);
    }
    length_type cyclic_transpose(length_type const &p0, length_type const &p1) const
    {
        const length_type diff(p1 - p0), half(world_size() / 2);
        if (diff > half)
        {
            return p0 + this->world_size();
        }
        else if (diff < -half)
        {
            return p0 - this->world_size();
        }
        else
        {
            return p0;
        }
    }

    position_type calculate_pair_CoM(
            position_type const& p1, position_type const& p2, 
            D_type const &D1, D_type const &D2)
    {
        position_type retval;
        const position_type p2t(this->cyclic_transpose(p2, p1));
        return ecell4::modulo(
                ecell4::divide(
                    ecell4::add( ecell4::multiply(p1, D2), ecell4::multiply(p2t, D1)),
                    ecell4::add(D1, D2)),
                world_size());
    }

    particle_id_pair_and_distance_list* check_overlap_adaptor(particle_id_pair_and_distance_list result) const
    {
        particle_id_pair_and_distance_list* retp = new particle_id_pair_and_distance_list(result);
        if (0 < retp->size()) 
        {
            particle_id_pair_and_distance_comparator c;
            std::sort(retp->begin(), retp->end(), c);
        }
        return retp;
    }

    particle_id_pair_and_distance_list* check_overlap(particle_shape_type const& s) const
    {
        return check_overlap_adaptor( list_particles_within_radius(s) );
    }
    particle_id_pair_and_distance_list* check_overlap(particle_shape_type const &s,
            particle_id_type const &ignore) const
    {
        return check_overlap_adaptor( list_particles_within_radius(s, ignore) );
    }

    particle_id_pair_and_distance_list* check_overlap(particle_shape_type const &s,
            particle_id_type const &ignore1, particle_id_type const &ignore2) const
    {
        return check_overlap_adaptor( list_particles_within_radius(s, ignore1, ignore2) );
    }
    
private:
    boost::scoped_ptr<ecell4::ParticleSpace> ps_;
    length_type world_size_;
    size_type size_;

    //ecell4::SerialIDGenerator<particle_id_type> pidgen_;
    particle_id_generator pidgen_;
    species_map species_map_;
    structure_map structure_map_;
    per_species_particle_id_set particle_pool_;
};

}   //egfrd

}   // ecell4

#endif
