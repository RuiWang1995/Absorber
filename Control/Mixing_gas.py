import math
import numpy as np
from Control.SolvePrEos import SolvePrEos
from Model.HeatCapacity import molar_heat_capacity
from Model.Viscosity import viscosity
from Model.binary_diffusivity import binary_diffusivity


class Mixing_gas(object):
    def __init__(self, Gb, properties, gas1, mole_fraction1, gas2, mole_fraction2, gas3, mole_fraction3, k12, k13, k23):
        self.Gb = Gb
        self.properties = properties
        self.k23 = k23
        self.k13 = k13
        self.k12 = k12
        self.mole_fraction3 = mole_fraction3
        self.gas3 = gas3
        self.mole_fraction2 = mole_fraction2
        self.gas2 = gas2
        self.mole_fraction1 = mole_fraction1
        self.gas1 = gas1

    def matrix(self):
        return np.array([[self.gas1, self.mole_fraction1], [self.gas2, self.mole_fraction2],
                         [self.gas3, self.mole_fraction3]])

    def k_matrix(self):
        return np.array([[0, self.k12, self.k13], [self.k12, 0, self.k23], [self.k13, self.k23, 0]])

    def a12(self, comp1, comp2, k12):
        preos1 = SolvePrEos(comp1, self.properties)
        preos2 = SolvePrEos(comp2, self.properties)
        a_comp1 = preos1.aValue()
        a_comp2 = preos2.aValue()
        return (1 - k12) * math.sqrt(a_comp1 * a_comp2)

    def new_avalue(self):
        summed_value = 0
        for i in range(0, 3):
            for j in range(0, 3):
                xixjaij = self.matrix()[i, 1] * self.matrix()[j, 1] * self.a12(self.matrix()[i, 0], self.matrix()[j, 0],
                                                                               self.k_matrix()[i, j])
                summed_value += xixjaij
        return summed_value

    def new_bvalue(self):
        summed_value = 0
        for i in range(0, 3):
            preos = SolvePrEos(self.matrix()[i, 0], self.properties)
            xibi = self.matrix()[i, 1] * preos.bValue()
            summed_value += xibi
        return summed_value

    def mix_molar_volume(self):
        from sympy import solveset, S
        from sympy.abc import x
        R = 8.314
        roots = solveset(self.properties.pressure - R * self.properties.temperature / (
                x - self.new_bvalue()) + self.new_avalue() / (
                                 x ** 2 + 2 * self.new_bvalue() * x - self.new_bvalue() ** 2), x, domain=S.Reals)
        if len(roots) == 1:
            return roots.args[0]
        else:
            return roots.args[2]

    def mix_molecular_weight(self):
        mw = 0
        for i in range(0, 3):
            mw += self.matrix()[i, 0].mol_weight * self.matrix()[i, 1]
        return mw

    # kg/m3
    def mix_density(self):
        return self.mix_molecular_weight() * 0.001 / self.mix_molar_volume()

    def mix_molar_heat_capacity(self):
        mix_molar_cp = 0
        for i in range(0, 3):
            mix_molar_cp += molar_heat_capacity(self.matrix()[i, 0], self.properties.temperature) * self.matrix()[i, 1]
        return mix_molar_cp

    def mix_molar_viscosity(self):
        mix_miu_part1 = 0
        mix_miu_part2 = 0
        for i in range(0, 3):
            mix_miu_part1 += viscosity(self.matrix()[i, 0], self.properties.temperature) * self.matrix()[i, 1] \
                             * self.matrix()[i, 0].mol_weight ** 0.5
        for i in range(0, 3):
            mix_miu_part2 += self.matrix()[i, 1] * self.matrix()[i, 0].mol_weight ** 0.5
        return mix_miu_part1 / mix_miu_part2

    def mix_diffusivity_matrix(self):
        matrix = np.array([[self.gas1, self.mole_fraction1, 0], [self.gas2, self.mole_fraction2, 0],
                           [self.gas3, self.mole_fraction3, 0]])
        for i in range(0, 3):
            sum_of_rest_diff = 0
            for j in range(0, 3):
                if j == i:
                    sum_of_rest_diff += 0
                else:
                    sum_of_rest_diff += matrix[j, 1] / binary_diffusivity(matrix[i, 0], matrix[j, 0],
                                                                          self.properties.temperature,
                                                                          self.properties.pressure)
            diffusivity_of_i = (1 - matrix[i, 1]) / sum_of_rest_diff
            matrix[i, 2] = diffusivity_of_i
        return matrix

    def average_Diff(self):
        average_D = 0
        for i in range(0, 3):
            average_D += self.mix_diffusivity_matrix()[i, 1] * self.mix_diffusivity_matrix()[i, 2]
        return average_D

    def gas_composition(self):
        matrix = self.matrix()
        new_matrix = np.array([[self.gas1, 0, 0], [self.gas2, 0, 0], [self.gas3, 0, 0]])
        for i in range(0, 3):
            if matrix[i, 0].name == 'AIR':
                yb = matrix[i, 1]
        for i in range(0, 3):
            new_matrix[i, 1] = matrix[i, 1]
            new_matrix[i, 2] = matrix[i, 1] / yb
        return new_matrix
