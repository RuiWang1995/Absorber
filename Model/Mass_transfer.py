from Model.Gases import Gases
from Model.HeatCapacity import mass_heat_capacity
from Model.heat_transfer_coeff import heat_transfer_coefficient


class mass_transfer(object):
    def __init__(self, Mixing_liquid, Mixing_gas, Packing):
        self.Mixing_liquid = Mixing_liquid
        self.Mixing_gas = Mixing_gas
        self.Packing = Packing

    # kg/m2/s
    def superficial_liquid_mass_velocity(self):
        return self.Mixing_liquid.velocity * self.Mixing_liquid.mix_density()

    def hydraulic_diameter(self):
        return 4 * self.Packing.void_fraction / self.Packing.specific_area

    def reynolds_number(self):
        return self.hydraulic_diameter() * self.superficial_liquid_mass_velocity() * self.Mixing_liquid.mix_viscosity()

    def froude_number(self):
        g = 9.81
        return self.superficial_liquid_mass_velocity() ** 2 * self.Packing.specific_area / (
                self.Mixing_liquid.mix_density() ** 2 * g)

    def weber_number(self):
        return self.superficial_liquid_mass_velocity() ** 2 / (
                self.Mixing_liquid.mix_density() * self.Packing.specific_area
                * self.Mixing_liquid.mix_surface_tension())

    def renewal_surface(self):
        return 0.74 * self.reynolds_number() ** 0.37

    # 1/m
    def specific_wet_area(self):
        return 0.25 * self.reynolds_number() ** 0.61 * 100

    def hl(self):
        return (12 / 9.81 * self.Mixing_liquid.mix_viscosity() / self.Mixing_liquid.mix_density() *
                self.Mixing_liquid.velocity * self.Packing.specific_area ** 2) ** (1 / 3)

    def mass_transfer_liquid(self):
        return (self.Mixing_liquid.diffusivity_CO2_aqueous_MEA() * self.renewal_surface()) ** 0.5

    # mol/m3
    def vapor_concentration(self):
        R = 8.314
        return self.Mixing_gas.properties.pressure / R / self.Mixing_gas.properties.temperature

    def mass_transfer_gas(self):
        matrix = self.Mixing_gas.mix_diffusivity_matrix()
        for i in range(0, 3):
            kg = self.Packing.cv / (self.Packing.void_fraction - self.hl()) ** 0.5 \
                 * (self.Packing.specific_area / self.hydraulic_diameter()) ** 0.5 * matrix[i, 2] \
                 * (self.Mixing_gas.Gb / self.vapor_concentration() / self.Packing.specific_area
                    / (self.Mixing_gas.mix_molar_viscosity() / self.Mixing_gas.mix_density())) ** 0.75 \
                 * (self.Mixing_gas.mix_molar_viscosity() / self.Mixing_gas.mix_density() / matrix[i, 2]) ** (1 / 3) \
                 / 8.314 / self.Mixing_gas.properties.temperature
            matrix[i, 2] = kg
        return matrix

    def average_kg(self):
        average_k = 0
        for i in range(0, 3):
            average_k += self.mass_transfer_gas()[i, 1] * self.mass_transfer_gas()[i, 2]
        return average_k * 8.314 * self.Mixing_gas.properties.temperature

    def heat_transfer_gas(self):
        temperature = self.Mixing_gas.properties.temperature
        for i in range(0, 3):
            if self.Mixing_gas.matrix()[i, 0].name == 'CO2':
                ya = self.Mixing_gas.matrix()[i, 1]
            if self.Mixing_gas.matrix()[i, 0].name == 'H2O':
                ys = self.Mixing_gas.matrix()[i, 1]
        AIR = self.Mixing_gas.gas1
        CO2 = self.Mixing_gas.gas2
        TB_co2 = 185.35
        TB_air = 78.8
        S1 = 1.5 * TB_co2
        S2 = 1.5 * TB_air
        S12 = (S1 * S2) ** 0.5
        Cp1 = mass_heat_capacity(CO2, temperature) * 0.2388458966
        Cp2 = mass_heat_capacity(AIR, temperature) * 0.2388458966
        M1 = CO2.mol_weight
        M2 = AIR.mol_weight
        lamda_co2 = heat_transfer_coefficient(CO2, temperature)
        lamda_air = heat_transfer_coefficient(AIR, temperature)
        miu1by2 = lamda_co2 / lamda_air * (Cp2 + 1.25 * 1.98 / M2) / (Cp1 + 1.25 * 1.98 / M1)
        miu2by1 = 1 / miu1by2
        A12 = 0.25 * (
                1 + (miu1by2 * (M2 / M1) ** 0.75 * (1 + S1 / temperature) / (1 + S2 / temperature)) ** 0.5) ** 2 * (
                      1 + S12 / temperature) / (1 + S1 / temperature)
        A21 = 0.25 * (
                1 + (miu2by1 * (M1 / M2) ** 0.75 * (1 + S2 / temperature) / (1 + S1 / temperature)) ** 0.5) ** 2 * (
                      1 + S12 / temperature) / (1 + S2 / temperature)
        A11 = 1
        A22 = 1
        lamda_ = ya * lamda_co2 / (A11 * ya + A12 * ys) + ys * lamda_air / (A21 * ya + A22 * ys)
        rho = 1/self.Mixing_gas.mix_molar_volume()
        return self.average_kg() * rho * self.Mixing_gas.mix_molar_heat_capacity() * \
               (lamda_ / rho / self.Mixing_gas.mix_molar_heat_capacity() /
                self.Mixing_gas.average_Diff()) ** (2 / 3)
