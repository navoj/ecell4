#ifndef __ECELL4_NETWORK_MODEL_HPP
#define __ECELL4_NETWORK_MODEL_HPP

// #include "get_mapper_mf.hpp"

#include <map>
#include <set>
#include <algorithm>
#include <iterator>

#include "types.hpp"
#include "Species.hpp"
#include "ReactionRule.hpp"
#include "Model.hpp"

#include "Context.hpp"


namespace ecell4
{

class NetworkModel
    : public Model
{
public:

    typedef Model base_type;
    typedef base_type::species_container_type species_container_type;
    typedef base_type::reaction_rule_container_type reaction_rule_container_type;
    typedef base_type::parameter_container_type parameter_container_type;

protected:

    typedef std::map<Species::serial_type,
                     std::vector<reaction_rule_container_type::size_type> >
        first_order_reaction_rules_map_type;
    typedef std::map<std::pair<Species::serial_type, Species::serial_type>,
                     std::vector<reaction_rule_container_type::size_type> >
        second_order_reaction_rules_map_type;

public:

    NetworkModel()
        : base_type(), species_attributes_(), reaction_rules_(),
        first_order_reaction_rules_map_(), second_order_reaction_rules_map_()
    {
        ;
    }

    virtual ~NetworkModel()
    {
        ;
    }

    // ModelTraits

    std::vector<ReactionRule> query_reaction_rules(const Species& sp) const;
    std::vector<ReactionRule> query_reaction_rules(
        const Species& sp1, const Species& sp2) const;

    Integer apply(const Species& pttrn, const Species& sp) const;
    std::vector<ReactionRule> apply(
        const ReactionRule& rr,
        const ReactionRule::reactant_container_type& reactants) const;

    Species apply_species_attributes(const Species& sp) const
    {
        for (species_container_type::const_iterator
            i(species_attributes_.begin()); i != species_attributes_.end(); ++i)
        {
            if (sp == *i)
            {
                Species retval(sp);
                retval.set_attributes(*i);
                return retval;
            }
        }
        return sp;
    }

    boost::shared_ptr<Model> expand(
        const std::vector<Species>& sp, const Integer max_itr,
        const std::map<Species, Integer>& max_stoich) const
    {
        return boost::shared_ptr<Model>(new NetworkModel(*this));
    }

    boost::shared_ptr<Model> expand(
        const std::vector<Species>& sp, const Integer max_itr) const
    {
        return boost::shared_ptr<Model>(new NetworkModel(*this));
    }

    boost::shared_ptr<Model> expand(const std::vector<Species>& sp) const
    {
        return boost::shared_ptr<Model>(new NetworkModel(*this));
    }

    // NetworkModelTraits

    void add_species_attribute(const Species& sp);
    bool has_species_attribute(const Species& sp) const;
    void remove_species_attribute(const Species& sp);

    void add_reaction_rule(const ReactionRule& rr);
    void remove_reaction_rule(const ReactionRule& rr);
    bool has_reaction_rule(const ReactionRule& rr) const;

    bool has_parameter(const Species::serial_type& name) const
    {
        parameter_container_type::const_iterator i(
            std::find(parameters_.begin(), parameters_.end(), Species(name)));
        return (i != parameters_.end());
    }

    const Species& get_parameter(const Species::serial_type& name) const
    {
        parameter_container_type::const_iterator i(
            std::find(parameters_.begin(), parameters_.end(), Species(name)));
        if (i != parameters_.end())
        {
            return (*i);
        }

        throw NotFound("Parameter not found.");
    }

    void add_parameter(const Species& sp)
    {
        parameter_container_type::iterator i(
            std::find(parameters_.begin(), parameters_.end(), sp));
        if (i != parameters_.end())
        {
            (*i).overwrite_attributes(sp);
        }
        else
        {
            parameters_.push_back(sp);
        }
    }

    const parameter_container_type& parameters() const
    {
        return parameters_;
    }

    // Optional functions

    const reaction_rule_container_type& reaction_rules() const
    {
        return reaction_rules_;
    }

    const species_container_type& species_attributes() const
    {
        return species_attributes_;
    }

protected:

    species_container_type species_attributes_;
    reaction_rule_container_type reaction_rules_;
    parameter_container_type parameters_;

    first_order_reaction_rules_map_type first_order_reaction_rules_map_;
    second_order_reaction_rules_map_type second_order_reaction_rules_map_;
};

} // ecell4

#endif /* __ECELL4_NETWORK_MODEL_HPP */
