class Gases(object):
    def __init__(self, name, mol_weight, tcritical, pcritical, afactor):
        self.afactor = afactor
        self.pcritical = pcritical
        self.tcritical = tcritical
        self.mol_weight = mol_weight
        self.name = name
