class Packing(object):
    def __init__(self, name, cl, cv, void_fraction, specific_area, packing_diameter):
        self.packing_diameter = packing_diameter
        self.specific_area = specific_area
        self.cl = cl
        self.cv = cv
        self.void_fraction = void_fraction
        self.name = name
