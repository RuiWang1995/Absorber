import numpy as np
import math


# J/mol CO2
def delta_Hr():
    return 88910


class Mixing_liquid(object):
    def __init__(self, velocity, properties, liquid1, mole_fraction1, concentration1, liquid2, mole_fraction2,
                 concentration2,
                 liquid3, mole_fraction3, concentration3):
        self.velocity = velocity
        self.concentration3 = concentration3
        self.concentration2 = concentration2
        self.concentration1 = concentration1
        self.mole_fraction3 = mole_fraction3
        self.liquid3 = liquid3
        self.mole_fraction2 = mole_fraction2
        self.liquid2 = liquid2
        self.mole_fraction1 = mole_fraction1
        self.liquid1 = liquid1
        self.properties = properties

    def matrix(self):
        return np.array([[self.liquid1, self.mole_fraction1, self.concentration1],
                         [self.liquid2, self.mole_fraction2, self.concentration2],
                         [self.liquid3, self.mole_fraction3, self.concentration3]])

    # mol/m2/s
    def molar_velocity(self):
        return self.velocity * (self.concentration1 + self.concentration2 + self.concentration3)

    def concentration_of_MEA(self):
        for i in range(0, 3):
            if self.matrix()[i, 0].name == 'MEA':
                return self.matrix()[i, 2]

    # kg/m3
    def mix_density(self):
        for i in range(0, 3):
            if self.matrix()[i, 0].name == 'MEA':
                C_MEA = self.matrix()[i, 2]
        TL = self.properties.temperature
        rhol = (-0.00311*TL**2+1.564*TL+807.466)+(1.358*10**(-8)*TL**3-1.296*10**(-5)*TL**2+4.1034*10**(-3)*TL-0.42931)*C_MEA+\
               (-4.19753*10**(-12)*TL**3+4.1*10**(-9)*TL**2-1.33528*10**(-6)*TL+1.45078*10**(-4))*C_MEA**2
        return rhol

    # N/m
    def mix_surface_tension(self):
        for i in range(0, 3):
            if self.matrix()[i, 0].name == 'MEA':
                x2 = self.matrix()[i, 1]
        zigma = 0.001 * math.exp((1.172791 - 0.023527 * math.log(x2)) *
                                 (3.396525 - 0.0022515 * (self.properties.temperature - 273.15)))
        return zigma

    # Pa.s
    def mix_viscosity(self):
        for i in range(0, 3):
            if self.matrix()[i, 0].name == 'MEA':
                C_MEA = self.matrix()[i, 2]
        miul = 0.001 * ((
                                0.0001222 * self.properties.temperature ** 2 - 0.08882 * self.properties.temperature + 16.48706) + (
                                1.9753 * 10 ** -10 * self.properties.temperature ** 3 - 1.84 * 10 ** -7 * self.properties.temperature ** 2 +
                                5.5521 * 10 ** -5 * self.properties.temperature - 5.3279 * 10 ** -3) * C_MEA
                        + (
                                -9.87654 * 10 ** -14 * self.properties.temperature ** 3 + 1.07556 * 10 ** -10 * self.properties.temperature ** 2
                                - 3.90871 * 10 ** -8 * self.properties.temperature + 4.7463 * 10 ** -6) * C_MEA ** 2)
        return miul

    # m2/s
    def mix_kinematic_viscosity(self):
        return self.mix_viscosity() / self.mix_density()

    # m2/s
    def diffusivity_CO2_aqueous_MEA(self):
        for i in range(0, 3):
            if self.matrix()[i, 0].name == 'MEA':
                C_MEA = self.matrix()[i, 2]
        Dl = 2.35 * 10 ** -6 * math.exp(
            -2119 / self.properties.temperature - C_MEA / 7796.94 + (C_MEA / 9614.1847) ** 2)
        return Dl

    # Pa.m3/mol
    def Henry_constant_CO2_aqueous_MEA(self):
        for i in range(0, 3):
            if self.matrix()[i, 0].name == 'MEA':
                C_MEA = self.matrix()[i, 2]
        HCO2 = 2.8249 * 10 ** 6 * math.exp(
            -2044 / self.properties.temperature + C_MEA / 413305 - (C_MEA / 34849.34) ** 2)
        return HCO2

    # m2/s
    def diffusivity_MEA_H2O(self):
        for i in range(0, 3):
            if self.matrix()[i, 0].name == 'MEA':
                mol_weight_MEA = self.matrix()[i, 0].mol_weight
        return 7.4 * 10 ** -12 * (2.6 * 18.01528) ** 0.5 * self.properties.temperature / (
                self.mix_viscosity() * 1000) / (mol_weight_MEA / 0.814051) ** 0.6

    # J/mol/K
    def molar_heat_capacity_aqueous_MEA(self):
        for i in range(0, 3):
            if self.matrix()[i, 0].name == 'MEA':
                x2 = self.matrix()[i, 1]
        A1 = -146.65 + 0.4898 * self.properties.temperature
        A2 = 22.6 - 0.079 * self.properties.temperature
        A3 = 10.876 - 0.0532 * self.properties.temperature
        X1 = x2
        X2 = 1 - x2
        CpE = X1 * X2 * (A1 + A2 * (X1 - X2) + A3 * (X1 - X2) ** 2)
        Cpl_2 = x2 * (78.578 + 0.2927 * self.properties.temperature) + (1 - x2) * (
                1000.182 / (
                1 - self.properties.temperature / 647.1081) + 53867.96 + 282.323 * self.properties.temperature -
                1.173532 * self.properties.temperature ** 2 + 0.001503196 * self.properties.temperature ** 3) / 1000 + CpE
        return Cpl_2

    def reaction_rate_constant(self):
        return 4.4 * 10 ** 8 * math.exp(-5400 / self.properties.temperature)
