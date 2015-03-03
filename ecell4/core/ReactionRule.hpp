#ifndef __ECELL4_REACTION_RULE_HPP
#define __ECELL4_REACTION_RULE_HPP

// #include <set>
#include <stdexcept>

#include "types.hpp"
#include "Species.hpp"
#include "Ratelaw.hpp"


namespace ecell4
{

class ReactionRule
{
public:

    /**
     * a type of the container of reactants
     * std::multiset allows multiple keys with equal values,
     * but looses the original order at the registration.
     * when changing this type into the ordered one,
     * please modify NetworkModel too.
     */
    typedef std::vector<Species> reactant_container_type;
    typedef std::vector<Species> product_container_type;

    // XXX unite {reactant|product}_coeff_container_type 
    // with above {reactant|product}_container_type??
    typedef std::vector<Real> reactant_coeff_container_type;
    typedef std::vector<Real> product_coeff_container_type;

public:

    ReactionRule()
        : k_(0), reactants_(), products_(), 
        reactant_coeffs_(reactants_.size(), 1.0), product_coeffs_(products_.size(), 1.0)
    {
        ;
    }

    ReactionRule(
        const reactant_container_type& reactants,
        const product_container_type& products)
        : k_(0), reactants_(reactants), products_(products), 
        reactant_coeffs_(reactants_.size(), 1.0), product_coeffs_(products_.size(), 1.0)
    {
        ;
    }

    ReactionRule(
        const reactant_container_type& reactants,
        const product_container_type& products,
        const Real& k)
        : k_(k), reactants_(reactants), products_(products),
        reactant_coeffs_(reactants_.size(), 1.0), product_coeffs_(products_.size(), 1.0)
    {
        ;
    }

    Real k() const
    {
        return k_;
    }

    const reactant_container_type& reactants() const
    {
        return reactants_;
    }

    const product_container_type& products() const
    {
        return products_;
    }

    void set_k(const Real& k)
    {
        if (k < 0)
        {
            throw std::invalid_argument("a kinetic rate must be positive.");
        }
        k_ = k;
    }

    void add_reactant(const Species& sp, const Real coefficient = 1.0)
    {
        reactants_.push_back(sp);
        reactant_coeffs_.push_back(coefficient);
    }

    void add_product(const Species& sp, const Real coefficient = 1.0)
    {
        products_.push_back(sp);
        product_coeffs_.push_back(coefficient);
    }

    const std::string as_string() const;
    Integer count(const reactant_container_type& reactants) const;
    std::vector<ReactionRule> generate(const reactant_container_type& reactants) const;

    /** Ratelaw related functions.
      */

    void set_ratelaw(const boost::shared_ptr<Ratelaw> ratelaw)
    {
        this->ratelaw_ = ratelaw;
    }

    boost::shared_ptr<Ratelaw> get_ratelaw() const
    {
        return this->ratelaw_.lock();
    }

    bool has_ratelaw() const
    {
        return !(this->ratelaw_.expired());
    }

    /**  Coefficient specific functions.
     */
    void set_reactant_coefficient(Integer const reactant_index, Integer const coefficient)
    {
        if (reactants_.size() < reactant_index) {
            throw std::invalid_argument("reactant_index is out of range");
        }
        reactant_coeffs_[reactant_index] = coefficient;
    }
    void set_product_coefficient(Integer product_index, Integer coefficient)
    {
        if (products_.size() < product_index) {
            throw std::invalid_argument("product_index is out of range");
        }
        product_coeffs_[product_index] = coefficient;
    }

protected:

    Real k_;
    reactant_container_type reactants_;
    product_container_type products_;

    reactant_coeff_container_type reactant_coeffs_;
    product_coeff_container_type product_coeffs_;

    boost::weak_ptr<Ratelaw> ratelaw_;
};

inline bool operator<(const ReactionRule& lhs, const ReactionRule& rhs)
{
    if (lhs.reactants() < rhs.reactants())
    {
        return true;
    }
    else if (lhs.reactants() > rhs.reactants())
    {
        return false;
    }
    return (lhs.products() < rhs.products());
}

inline bool operator==(const ReactionRule& lhs, const ReactionRule& rhs)
{
    return ((lhs.reactants() == rhs.reactants())
            && (lhs.products() == rhs.products()));
}

inline bool operator!=(const ReactionRule& lhs, const ReactionRule& rhs)
{
    return !(lhs == rhs);
}

} // ecell4

#endif /* __ECELL4_REACTION_RULE_HPP */
