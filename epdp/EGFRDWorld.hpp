#ifndef EGFRD_WORLD_HPP
#define EGFRD_WORLD_HPP

#include <map>
#include <boost/lexical_cast.hpp>
#include <boost/array.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_same.hpp>
#include "exceptions.hpp"
#include "generator.hpp"
#include "filters.hpp"
#include "Particle.hpp"
#include "ParticleID.hpp"
#include "SpeciesTypeID.hpp"
#include "SpeciesInfo.hpp"
#include "SerialIDGenerator.hpp"
#include "Transaction.hpp"
#include "Structure.hpp"
#include "Surface.hpp"
#include "Region.hpp"
#include "geometry.hpp"
//#include "GSLRandomNumberGenerator.hpp"
#include "Point.hpp" // XXX: workaround. should be removed later.
#include "utils/pair.hpp"

#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/Species.hpp>

#include <ecell4/core/types.hpp>
#include <ecell4/core/SerialIDGenerator.hpp>
#include <ecell4/core/Position3.hpp>
#include <ecell4/core/Particle.hpp>
#include <ecell4/core/ParticleSpace.hpp>

#include "World.hpp"

namespace ecell4 {

namespace egfrd {
// set of the species informations which are reqruired by EGFRDSimulator
struct MoleculeInfo
{
    const ecell4::Species::serial_type sid;
    const Real radius;
    const Real D;
};
}   //egfrd
}   // ecell4

template<typename Ttraits_>
class EGFRDWorld
{
    /* species id type should be fixed */
public:
    typedef Ttraits_ traits_type;
    typedef typename traits_type::length_type length_type;
    typedef typename traits_type::species_type species_type;
    typedef typename traits_type::position_type position_type;
    typedef typename traits_type::size_type size_type;
    typedef typename traits_type::species_id_type species_id_type;
    typedef typename traits_type::structure_id_type structure_id_type;
    typedef typename traits_type::structure_type structure_type;

    typedef ecell4::ParticleID particle_id_type;
    typedef ecell4::Particle particle_type;
    typedef ecell4::ParticleSpace particle_space_type;

    typedef std::pair<const particle_id_type, particle_type> particle_id_pair;
    typedef std::pair<particle_id_pair, length_type> particle_id_pair_and_distance;
    typedef sized_iterator_range<
        typename particle_space_type::particle_container_type::const_iterator> particle_id_pair_range;

protected:
    typedef std::map<species_id_type, species_type> species_map;
    typedef std::map<structure_id_type, boost::shared_ptr<structure_type> > structure_map;
    typedef std::set<particle_id_type> particle_id_set;
    typedef std::map<species_id_type, particle_id_set> per_species_particle_id_set;
    typedef select_second<typename structure_map::value_type> surface_second_selector_type;
    
public:
    typedef boost::transform_iterator<surface_second_selector_type, typename structure_map::const_iterator> surface_iterator;
    typedef sized_iterator_range<surface_iterator> structures_range;

public:
    EGFRDWorld(length_type world_size = 1, size_type size = 1) :
        world_size_(world_size),
        size_(size),
        ps_(new ecell4::ParticleSpaceVectorImpl(
                    position_type(world_size,world_size,world_size) ))
    {}

    size_type num_molecules(const ecell4::Species &sp) const
    {   return num_particles(sp);   }

    size_type num_molecules(const species_id_type &sid) const
    {   return num_particles(ecell4::Species &sp);  }

    size_type num_particles() const
    {   return ps_->num_particles();    }

    size_type world_size() const 
    {   return this->world_size();  }

    inline length_type distance(const position_type &pos1, const position_type &pos2) const 
    {   return ps_->distance(pos1, pos2);   }

    // epdp
    particle_id_pair new_particle(species_id_type const &sid, position_type &pos)
    {
        
    }

    // particlecontainer
    position_type apply_boundary(position_type const &v) const
    {   return traits_type::apply_boundary(v, world_size());    }

    length_type apply_boundary(length_type const &v) const
    {   return traits_type::apply_boundary(v, world_size());   }

    position_type cyclic_transpose(position_type const &p0, position_type const &p1) const
    {   return traits_type::cyclic_transpose(p0, p1, world_size());  }

    length_type cyclic_transpose(length_type const &p0, length_type const &p1) const
    {   return traits_type::cyclic_transpose(p0, p1, world_size()); }

    particle_id_pair get_particle(particle_id_type const &pid, bool &found) const
    {
        if (ps_->has_particle(pid) == false) {
            found = false;
            return particle_id_pair();
        }
        found = true;
        return ps_->get_particle(pid);
    }

    particle_id_pair get_particle(particle_id_type const &pid) const
    {
        return ps_->get_particle(pid);
    }

    /*
    particle_id_pair_generator *get_particles() const
    {
    }
    */

    particle_id_pair_range get_particles_range() const
    {
        return particle_id_pair_range(
                ps_->list_particles().begin(), 
                ps_->list_particles().end(), 
                ps_->num_particles() );
    }

    bool has_particle(particle_id_type const &pid) const
    {
        return ps_->has_particle(pid);
    }
    
    // bd
    bool update_particle(const /*ecell4::ParticleID */particle_id_type &pid,
            const /* ecell4::Particle */ particle_type &p)
    {   return ps_->update_particle(pid, p);    }

    // epdp
    bool update_particle(particle_id_pair const &pi_pair) 
    {   return ps_->update_particle(pi_pair.first, pi_pair.second); }

    // epdp bd
    void remove_particle(const particle_id_type &pid)
    {
        ps_->remove_particle(pid);
    }

    ecell4::egfrd::MoleculeInfo get_molecule_info(const ecell4::Species &sp) const 
    {
        const ecell4::Species::serial_type serial(sp.name());
        const Real radius(std::atof(sp.get_attribute("radius").c_str()));
        const Real D(std::atof(sp.get_attribute("D").c_str()));
        ecell4::egfrd::MoleculeInfo info = {serial, sp.serial(), radius, D};
        return info;
    }

    // species_type should be fixed
    void add_species(species_id_type const &sid, MoleculeInfo const &info, structure_id_type structure_id = structure_id_type("world") )
    {
        //XXX
        species_type sp(sid, info.D, info.radius, structure_id);
        species_map_[sid] = sp;
        particle_pool_[sid] = particle_id_set();
    }

    void add_species(species_type const &species)
    {
        species_map_[species.id()] = species;
        particle_pool_[species.id()] = particle_id_set();
    }

    species_type const &get_species(species_id_type const &id) const
    {
        // XXX
        typename species_map::const_iterator i(species_map_.find(id));
        if (species_map_.end() == i)
        {
            throw not_found("unknown species");
        }
        return i->second;
    }

    species_range get_species() const
    {
    }



    //bd
    std::vector<std::pair<particle_id_pair, ecell4::Real> >
    list_particles_within_radius(
            const position_type &pos, const length_type &radius) const
    {   return ps_->list_particles_within_radius(pos, radius);  }

    // bd
    std::vector<std::pair<particle_id_pair, ecell4::Real> >
    list_particles_within_radius(
            const position_type &pos, const length_type &radius,
            const particle_id_type &ignore) const
    {   return ps_->list_particles_within_radius(pos, radius, ignore);  }
    
    // bd
    std::vector<std::pair<particle_id_pair, ecell4::Real> >
    list_particles_within_radius(
            const position_type &pos, const length_type &radius,
            const particle_id_type &ignore1, const particle_id_type &ignore2)
    {   return ps_->list_particles_within_radius(pos, radius, ignore1, ignore2); }

    // CompartmentSpace Traits
    // bd
    void add_molecules(const ecell4::Species &sp, const Integer &num)
    {

    }

    bool add_structure(boost::shared_ptr<structure_type> surface)
    {
        return structure_map_.insert(std::make_pair(surface->id(), surface)).second;
    }

    boost::shared_ptr<structure_type> get_structure(structure_id_type const &id) const
    {
        typename structure_map::const_iterator i(structure_map_.find(id));
        if (structure_map_.end() == i)
        {
            throw not_found("Unknown surface");
        }
        return i->second;
    }

    structures_range get_structures() const 
    {
        return structures_range(
                surface_iterator(structure_map_.begin(), surface_second_selector_type()),
                surface_iterator(structure_map_.end(), surface_second_selector_type() ), 
                structure_map_.size());
    }
    particle_id_set get_particle_ids(species_id_type const &sid) const
    {
        typename per_species_particle_id_set::const_iterator i(
                particle_pool_.find(sid));
        if (i == particle_pool_.end())
        {
            throw(not_found("Unknown species"));
        }
        return i->second;
    }

private:
    boost::scoped_ptr<ecell4::ParticleSpace> ps_;
    length_type world_size_;
    size_type size_;

    ecell4::SerialIDGenerator<particle_id_type> pidgen_;
    species_map species_map_;
    structure_map structure_map_;
    per_species_particle_id_set particle_pool_;
};

#endif 
