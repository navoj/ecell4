#ifndef ECELL4_EGFRD_WORLD_HPP
#define ECELL4_EGFRD_WORLD_HPP

#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/Species.hpp>
#include <ecell4/core/types.hpp>
#include <ecell4/core/Position3.hpp>
#include <ecell4/core/ParticleSpace.hpp>
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

template <typename Ttraits_>
class EGFRDWorld
{
public:
    typedef Ttraits_ traits_type;
    typedef typename traits_type::length_type length_type;
    typedef typename traits_type::molecule_info molecule_info;  //XXX
    typedef typename traits_type::position_type position_type;
    typedef typename traits_type::size_type size_type;
    typedef typename traits_type::species_id_type species_id_type;
    typedef typename traits_type::structure_id_type structure_id_type;
    typedef typename traits_type::structure_type structure_type;
    typedef typename traits_type::particle_type particle_type;
    typedef typename traits_type::particle_id_type particle_id_type;
    typedef typename traits_type::particle_id_generator particle_id_generator;
    typedef typename traits_type::particle_shape_type particle_shape_type;
    typedef typename std::pair<const particle_id_type, particle_type> particle_id_pair;

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
        this->partile_pool_[pp.second.sid()].erase(id);
        ps_->remove_particle(id);
        return true;
    }

    particle_id_set get_particle_ids(species_id_type const &sid) const
    {
        typename per_species_particle_id_set::const_iterator i(
                particle_pool_find(sid));
        if (i == this->particle_pool_.end())
        {
            throw not_found("Unknown species");
        }
        return i->second;
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


    length_type distance(const position_type &pos1, const position_type &pos2) const
    {
        return ps_->distance(pos1, pos2);
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
    
    molecule_info get_species(species_id_type const &sid) const
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
