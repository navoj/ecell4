#ifndef __ECELL4_NETFREE_MODEL_HPP
#define __ECELL4_NETFREE_MODEL_HPP

// #include "get_mapper_mf.hpp"

#include <map>
#include <set>
#include <algorithm>
#include <iterator>
#include <boost/shared_ptr.hpp>

#include "types.hpp"
#include "Species.hpp"
#include "ReactionRule.hpp"
#include "Model.hpp"

#include "Context.hpp"
#include "NetworkModel.hpp"


namespace ecell4
{

class NetfreeModel
    : public Model
{
public:

    typedef Model base_type;
    typedef base_type::species_container_type species_container_type;
    typedef base_type::reaction_rule_container_type reaction_rule_container_type;
    typedef base_type::parameter_container_type parameter_container_type;

public:

    NetfreeModel()
        : base_type(), species_attributes_(), reaction_rules_()
    {
        ;
    }

    virtual ~NetfreeModel()
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
            if (spmatch(*i, sp))
            {
                Species retval(sp);
                retval.set_attributes(*i);
                return retval;
            }
        }
        return sp;
    }

    // NetfreeModelTraits

    void add_species_attribute(const Species& sp);
    bool has_species_attribute(const Species& sp) const;
    bool has_species_attribute_exact(const Species& sp) const;
    void remove_species_attribute(const Species& sp);

    void add_reaction_rule(const ReactionRule& rr);
    void remove_reaction_rule(const ReactionRule& rr);
    bool has_reaction_rule(const ReactionRule& rr) const;

    // Optional functions

    const reaction_rule_container_type& reaction_rules() const
    {
        return reaction_rules_;
    }

    const species_container_type& species_attributes() const
    {
        return species_attributes_;
    }

    boost::shared_ptr<Model> expand(
        const std::vector<Species>& sp, const Integer max_itr,
        const std::map<Species, Integer>& max_stoich) const;
    boost::shared_ptr<Model> expand(
        const std::vector<Species>& sp, const Integer max_itr) const;
    boost::shared_ptr<Model> expand(const std::vector<Species>& sp) const;

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

protected:

    species_container_type species_attributes_;
    reaction_rule_container_type reaction_rules_;
    parameter_container_type parameters_;
};

namespace extras
{

std::pair<boost::shared_ptr<NetworkModel>, bool> generate_network_from_netfree_model(
    const NetfreeModel& nfm, const std::vector<Species>& seeds, const Integer max_itr,
    const std::map<Species, Integer>& max_stoich);

inline std::pair<boost::shared_ptr<NetworkModel>, bool> generate_network_from_netfree_model(
    const NetfreeModel& nfm, const std::vector<Species>& seeds, const Integer max_itr)
{
    const std::map<Species, Integer> max_stoich;
    return generate_network_from_netfree_model(
        nfm, seeds, max_itr, max_stoich);
}

} // extras

} // ecell4

#endif /* __ECELL4_NETFREE_MODEL_HPP */
