from Control.Iterations import iteration
from Control.Mixing_gas import Mixing_gas
from Control.Mixing_liquid import Mixing_liquid
from Model.Column import Column
from Model.Gases import Gases
from Model.HeatCapacity import mass_heat_capacity, molar_heat_capacity
from Model.Liquid import Liquid
from Model.Mass_transfer import mass_transfer
from Model.Packing import Packing
from Model.Properties import Properties
from Control.SolvePrEos import SolvePrEos
from Model.binary_diffusivity import binary_diffusivity
import math

from Model.heat_transfer_coeff import heat_transfer_coefficient
from Model.stage import stage

AIR = Gases('AIR', 28.9647, 132.63, 3785800, 0.036)
CO2 = Gases('CO2', 44, 304.1, 7380000, 0.225)
H2O = Gases('H2O', 18, 647.15, 22050000, 0.344)
MEA = Liquid('MEA', 61.08)
H2O_L = Liquid('H2O_L', 18)
MEACOO = Liquid('MEACOO', 104)
PT_gas = Properties(292.15, 101325)
PT_liquid = Properties(312.47, 101325)
mixture_gas = Mixing_gas(14.8, PT_gas, AIR, 0.8050, CO2, 0.1950, H2O, 0, 0, 0, 0)
mixture_liquid = Mixing_liquid(0.0037, PT_liquid, MEA, 0.0023, 130, H2O, 0.9647, 54714, MEACOO, 0.033, 1870)
berl_saddle = Packing('berl_saddle', 1.364, 0.232, 0.65, 545, 0.0127)
column = Column(0.1, 6.55)
mass_transfer1 = mass_transfer(mixture_liquid, mixture_gas, berl_saddle)
stage1 = stage(mixture_gas, mixture_liquid, mass_transfer1)
it_test = iteration(0.01, column, berl_saddle, stage1)

print(mass_transfer1.heat_transfer_gas())
print(mixture_gas.mix_molar_heat_capacity())
print(mass_transfer1.average_kg())
print(mixture_gas.mix_density())
print(mixture_gas.average_Diff())
print(mass_heat_capacity(CO2, 298.15) * 0.238845)
print(mass_heat_capacity(AIR, 298.15) * 0.238845)
print(heat_transfer_coefficient(CO2, 298.15))
print(heat_transfer_coefficient(AIR, 298.15))
print(it_test.going_up()[0].Mixing_gas.gas_composition())




