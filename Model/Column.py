class Column(object):
    def __init__(self, diameter, height):
        self.diameter = diameter
        self.height = height

    def cross_section_area(self):
        return 3.1415926/4*self.diameter**2
