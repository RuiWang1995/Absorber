import math

from Control.Mixing_gas import Mixing_gas
from Control.Mixing_liquid import delta_Hr, Mixing_liquid
from Model.HeatCapacity import molar_heat_capacity
from Model.Mass_transfer import mass_transfer
from Model.Properties import Properties
from Model.Vapor_pressure import water_vapor_pressure
from Model.stage import stage


class iteration(object):
    def __init__(self, dz, column, packing, initial_stage):
        self.initial_stage = initial_stage
        self.packing = packing
        self.column = column
        self.stages = [initial_stage]
        self.dz = dz

    def flow_rate(self):
        return self.initial_stage.Mixing_liquid.velocity * self.column.cross_section_area()

    def contact_area(self, specific_wet_area):
        return self.dz * 3.1415926 / 4 * self.column.diameter ** 2 * specific_wet_area

    def number_of_stages(self):
        return int(self.column.height / self.dz)

    def number_of_iterations(self):
        return self.number_of_stages() + 1

    # GLI refers to gas_liquid_interface
    def iterate_for_GLI(self, previous_stage):
        ya = previous_stage.Mixing_gas.gas_composition()[1, 1]
        ys = previous_stage.Mixing_gas.gas_composition()[2, 1]
        yaj = previous_stage.Mixing_gas.gas_composition()[1, 1]
        yai = 0.01
        pressure = previous_stage.Mixing_gas.properties.pressure
        p_H2O = water_vapor_pressure(previous_stage.Mixing_gas.properties.temperature)
        HCO2 = previous_stage.Mixing_liquid.Henry_constant_CO2_aqueous_MEA()
        Dl = previous_stage.Mixing_liquid.diffusivity_CO2_aqueous_MEA()
        D_MEA = previous_stage.Mixing_liquid.diffusivity_MEA_H2O()
        k2 = previous_stage.Mixing_liquid.reaction_rate_constant()
        C_MEA = previous_stage.Mixing_liquid.concentration_of_MEA()
        mtc = mass_transfer(previous_stage.Mixing_liquid, previous_stage.Mixing_gas, self.packing)
        aph = mtc.specific_wet_area()
        a_contact = self.contact_area(aph)
        kl = mtc.mass_transfer_liquid()
        kg_CO2 = mtc.mass_transfer_gas()[1, 2]
        kg_H2O = mtc.mass_transfer_gas()[2, 2]
        pai = yai * pressure
        paj = yaj * pressure
        while abs(pai - paj) > 0.00001:
            yai = yaj
            pai = yai * pressure
            C_co2i = pai / HCO2
            b = 2
            M = Dl * k2 * C_MEA / (kl ** 2)
            Ei = 1 + (C_MEA * D_MEA / (b * Dl * C_co2i))
            E1 = (M ** 0.5) / math.tanh(M ** 0.5)
            E = 1 + 1 / (((1 / (Ei - 1)) ** 1.35 + (1 / (E1 - 1)) ** 1.35) ** (1 / 1.35))
            Nj = E * kl * (C_co2i - 0)
            paj = pressure * ya / (1 + kl * E / kg_CO2 / HCO2)
            yaj = paj / pressure
        deltaYa = -kg_CO2 * aph * pressure * (ya - yai) / previous_stage.Mixing_gas.Gb * self.dz
        deltaYs = kg_H2O * (p_H2O - pressure * ys) * aph / previous_stage.Mixing_gas.Gb * self.dz
        deltaC_MEA = 2 * Nj * a_contact / self.flow_rate()
        deltaC_H2O = previous_stage.Mixing_gas.Gb * deltaYs * self.column.cross_section_area() / self.flow_rate()
        deltaC_MEACOO_MEAH = -2 * Nj * a_contact / self.flow_rate()

        TGin = previous_stage.Mixing_gas.properties.temperature
        TLout = previous_stage.Mixing_liquid.properties.temperature
        TG_out = TGin + 0.0002
        TGout = 0
        cpmxg = previous_stage.Mixing_gas.mix_molar_heat_capacity()
        Cpl_2 = previous_stage.Mixing_liquid.molar_heat_capacity_aqueous_MEA()
        CpCO2_2 = molar_heat_capacity(previous_stage.Mixing_gas.gas2, previous_stage.Mixing_gas.properties.temperature)
        CpH2O_2 = molar_heat_capacity(previous_stage.Mixing_gas.gas3, previous_stage.Mixing_gas.properties.temperature)
        LM = previous_stage.Mixing_liquid.molar_velocity()
        hG = mtc.heat_transfer_gas()
        while abs(TG_out - TGout) > 0.01:
            latent = -3.068e-4 * TLout ** 3 + 0.27778 * TLout ** 2 - 126.7 * TLout + 65197
            TGout = TG_out
            dTG = TGout - TGin
            dYs = deltaYs
            HR = delta_Hr()
            dYa = deltaYa
            TLin = TLout + previous_stage.Mixing_gas.Gb / LM / Cpl_2 * (cpmxg * dTG + latent * dYs + HR * dYa)
            hG_ = -1 * previous_stage.Mixing_gas.Gb * (CpCO2_2 * dYa / self.dz + CpH2O_2 * dYs / self.dz) / (
                    1 - math.exp(
                previous_stage.Mixing_gas.Gb * (CpCO2_2 * dYa / self.dz + CpH2O_2 * dYs / self.dz) / (hG * aph))) / aph
            Q = hG_ * aph * (TLin - TGin)
            TG_out = Q / (previous_stage.Mixing_gas.Gb * cpmxg) * self.dz + TGin
        return [Nj, deltaYa, deltaYs, deltaC_MEA, deltaC_H2O, deltaC_MEACOO_MEAH, TLin, TG_out]

    def next_stage(self, previous_stage):
        Ya = previous_stage.Mixing_gas.gas_composition()[1, 1] + self.iterate_for_GLI(previous_stage)[1]
        Ys = previous_stage.Mixing_gas.gas_composition()[2, 1] + self.iterate_for_GLI(previous_stage)[2]
        yb = 1 / (Ya + Ys + 1)
        ya = Ya / (Ya + Ys + 1)
        ys = Ys / (Ya + Ys + 1)
        C_MEA = previous_stage.Mixing_liquid.matrix()[0, 2] + self.iterate_for_GLI(previous_stage)[3]
        C_H2O = previous_stage.Mixing_liquid.matrix()[1, 2] + self.iterate_for_GLI(previous_stage)[4]
        C_MEACOO_MEAH = previous_stage.Mixing_liquid.matrix()[2, 2] + self.iterate_for_GLI(previous_stage)[5]
        y_MEA = C_MEA / (C_MEA + C_H2O + C_MEACOO_MEAH)
        y_H2O = C_H2O / (C_MEA + C_H2O + C_MEACOO_MEAH)
        y_MEACOO_MEAH = C_MEACOO_MEAH / (C_MEA + C_H2O + C_MEACOO_MEAH)
        TL = self.iterate_for_GLI(previous_stage)[6]
        TG = self.iterate_for_GLI(previous_stage)[7]
        properties_gas = Properties(TG, previous_stage.Mixing_gas.properties.pressure)
        properties_liquid = Properties(TL, previous_stage.Mixing_liquid.properties.pressure)
        gas1 = previous_stage.Mixing_gas.gas1
        gas2 = previous_stage.Mixing_gas.gas2
        gas3 = previous_stage.Mixing_gas.gas3
        Gb = previous_stage.Mixing_gas.Gb
        current_mixing_gas = Mixing_gas(Gb, properties_gas, gas1, yb, gas2, ya, gas3, ys,
                                        previous_stage.Mixing_gas.k12, previous_stage.Mixing_gas.k13,
                                        previous_stage.Mixing_gas.k23)
        velocity = previous_stage.Mixing_liquid.velocity
        liquid1 = previous_stage.Mixing_liquid.liquid1
        liquid2 = previous_stage.Mixing_liquid.liquid2
        liquid3 = previous_stage.Mixing_liquid.liquid3
        current_mixing_liquid = Mixing_liquid(velocity, properties_liquid, liquid1, y_MEA, C_MEA,
                                              liquid2, y_H2O, C_H2O, liquid3, y_MEACOO_MEAH, C_MEACOO_MEAH)
        current_packing = previous_stage.Mass_transfer.Packing
        current_mass_transfer = mass_transfer(current_mixing_liquid, current_mixing_gas, current_packing)
        current_stage = stage(current_mixing_gas, current_mixing_liquid, current_mass_transfer)
        return current_stage

    def going_up(self):
        self.stages[0] = self.initial_stage
        for i in range(0, 2):
            self.stages.append(self.next_stage(self.stages[i]))
        return self.stages
