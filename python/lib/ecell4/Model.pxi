from cython.operator cimport dereference as deref, preincrement as inc
from cython cimport address

from libcpp.vector cimport vector


cdef class Model:

    def __cinit__(self):
        # self.thisptr = new shared_ptr[Cpp_Model](
        #     new Cpp_Model()) # cannot allocate an object of abstract type
        self.thisptr = new shared_ptr[Cpp_Model](
            <Cpp_Model*>(new Cpp_NetworkModel()))

    def __dealloc__(self):
        del self.thisptr

    def add_species_attribute(self, Species sp):
        self.thisptr.get().add_species_attribute(deref(sp.thisptr))

    def has_species_attribute(self, Species sp):
        return self.thisptr.get().has_species_attribute(deref(sp.thisptr))

    def remove_species_attribute(self, Species sp):
        self.thisptr.get().remove_species_attribute(deref(sp.thisptr))

    def add_reaction_rule(self, ReactionRule rr):
        self.thisptr.get().add_reaction_rule(deref(rr.thisptr))

    def remove_reaction_rule(self, ReactionRule rr):
        self.thisptr.get().remove_reaction_rule(deref(rr.thisptr))

    def has_reaction_rule(self, ReactionRule rr):
        return self.thisptr.get().has_reaction_rule(deref(rr.thisptr))

    def apply_species_attributes(self, Species sp):
        cdef Cpp_Species retval = self.thisptr.get().apply_species_attributes(
            deref(sp.thisptr))
        return Species_from_Cpp_Species(address(retval))

    def num_reaction_rules(self):
        return self.thisptr.get().num_reaction_rules()

    def reaction_rules(self):
        cdef vector[Cpp_ReactionRule] c_rr_vector = self.thisptr.get().reaction_rules()
        retval = []
        cdef vector[Cpp_ReactionRule].iterator it = c_rr_vector.begin()
        while it != c_rr_vector.end():
            retval.append(ReactionRule_from_Cpp_ReactionRule(
                <Cpp_ReactionRule*>(address(deref(it)))))
            inc(it)
        return retval

    def species_attributes(self):
        cdef vector[Cpp_Species] species = self.thisptr.get().species_attributes()
        retval = []
        cdef vector[Cpp_Species].iterator it = species.begin()
        while it != species.end():
            retval.append(Species_from_Cpp_Species(
                <Cpp_Species*>(address(deref(it)))))
            inc(it)
        return retval

    def list_species(self):
        cdef vector[Cpp_Species] species = self.thisptr.get().list_species()
        retval = []
        cdef vector[Cpp_Species].iterator it = species.begin()
        while it != species.end():
            retval.append(Species_from_Cpp_Species(
                <Cpp_Species*>(address(deref(it)))))
            inc(it)
        return retval

    def query_reaction_rules(self, Species sp1, Species sp2 = None):
        cdef vector[Cpp_ReactionRule] rules
        if sp2 is None:
            rules = self.thisptr.get().query_reaction_rules(
                deref(sp1.thisptr))
        else:
            rules = self.thisptr.get().query_reaction_rules(
                deref(sp1.thisptr), deref(sp2.thisptr))
        retval = []
        cdef vector[Cpp_ReactionRule].iterator it = rules.begin()
        while it != rules.end():
            retval.append(ReactionRule_from_Cpp_ReactionRule(
                <Cpp_ReactionRule*>(address(deref(it)))))
            inc(it)
        return retval

    def add_species_attributes(self, attrs):
        cdef vector[Cpp_Species] species
        for sp in attrs:
            species.push_back(deref((<Species>sp).thisptr))
        self.thisptr.get().add_species_attributes(species)

    def add_reaction_rules(self, rrs):
        cdef vector[Cpp_ReactionRule] reaction_rules
        for rr in rrs:
            reaction_rules.push_back(deref((<ReactionRule>rr).thisptr))
        self.thisptr.get().add_reaction_rules(reaction_rules)

    def expand(self, seeds, max_itr=None, max_stoich=None):
        cdef vector[Cpp_Species] _seeds
        cdef map[Cpp_Species, Integer] _max_stoich
        for sp in seeds:
            if not isinstance(sp, Species):
                raise ValueError(
                    'seeds must be given as a list of Species.'
                    + ' {0} given.'.format(repr(sp)))
            _seeds.push_back(deref((<Species>sp).thisptr))

        if max_stoich is not None:
            for sp, n in max_stoich.items():
                _max_stoich[deref((<Species>sp).thisptr)] = <Integer>n
            return Model_from_Cpp_Model(
                self.thisptr.get().expand(_seeds, <Integer>max_itr, _max_stoich))
        elif max_itr is not None:
            return Model_from_Cpp_Model(
                self.thisptr.get().expand(_seeds, <Integer>max_itr))
        else:
            return Model_from_Cpp_Model(
                self.thisptr.get().expand(_seeds))

cdef Model Model_from_Cpp_Model(shared_ptr[Cpp_Model] m):
    r = Model()
    r.thisptr.swap(m)
    return r
