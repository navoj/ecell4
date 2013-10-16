import copy


class Model(object):

    def __init__(self):
        self.__species_attributes = []
        self.__reaction_rules = [[], []]

    def query_reaction_rules(self, sp1, sp2=None):
        retval = []
        if sp2 is None:
            for rr in self.__reaction_rules[0]:
                tmp = rr.generate([sp1])
                if tmp is not None:
                    retval.extend(tmp)
        else:
            for rr in self.__reaction_rules[1]:
                tmp = rr.generate([sp1, sp2])
                if tmp is not None:
                    retval.extend(tmp)

                tmp = rr.generate([sp2, sp1])
                if tmp is not None:
                    retval.extend(tmp)
        return retval

    def add_species_attribute(self, sp):
        self.__species_attributes.append(copy.deepcopy(sp))

    def has_species_attribute(self, sp):
        for sp2 in self.__species_attributes:
            contexts = sp2.match(sp)
            if contexts is not None and len(contexts) != 0:
                return True
        else:
            return False

    def remove_species_attribute(self, sp):
        raise NotImplementedError

    def apply_species_attributes(self, sp):
        for sp2 in self.__species_attributes:
            contexts = sp2.match(sp)
            if contexts is not None and len(contexts) != 0:
                for key, value in sp2.attributes().items():
                    sp.set_attribute(key, value)
                return

    def add_reaction_rule(self, rr):
        if len(rr.reactants()) == 1:
            self.__reaction_rules[0].append(copy.deepcopy(rr))
        elif len(rr.reactants()) == 2:
            self.__reaction_rules[1].append(copy.deepcopy(rr))
        else:
            raise NotImplementedError

    def remove_reaction_rule(self, rr):
        raise NotImplementedError

    def has_reaction_rule(self, rr):
        raise NotImplementedError
