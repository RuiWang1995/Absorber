import math


# return J/mol/K
def diffusion_volume(gas):
    if gas.name == 'CO2':
        return 26.9
    if gas.name == 'AIR':
        return 20.1
    if gas.name == 'H2O':
        return 12.7


# m2/s
def binary_diffusivity(gas1, gas2, temperature, pressure):
    return 10 ** -7 * temperature ** 1.75 * math.sqrt(1 / gas1.mol_weight + 1 / gas2.mol_weight) / (
                pressure * 0.986923 / 10 ** 5) \
           / (diffusion_volume(gas1) ** (1 / 3) + diffusion_volume(gas2) ** (1 / 3)) ** 2
