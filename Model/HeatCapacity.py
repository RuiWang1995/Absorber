import math


# return J/mol/K
def molar_heat_capacity(gas, temperature):
    if gas.name == 'CO2':
        return 4.184 * (7.54056 + 7.51625 * (1442.7 / temperature / math.sinh(1442.7 / temperature)) ** 2 +
                        5.38023 * (647.502 / temperature / math.cosh(647.502 / temperature)) ** 2)
    if gas.name == 'AIR':
        return 7 * 10 ** (-12) * temperature ** 4 - 3 * 10 ** (-8) * temperature ** 3 + \
               4 * 10 ** (-5) * temperature ** 2 - 0.0144 * temperature + 31.648
    if gas.name == 'H2O':
        return 4.184 * (7.97183 + 6.27078 * (2572.63 / temperature / math.sinh(2572.63 / temperature)) ** 2 +
                        2.051 * (1156.72 / temperature / math.cosh(1156.72 / temperature)) ** 2)


# return J/g/K
def mass_heat_capacity(gas, temperature):
    if gas.name == 'CO2':
        return 4.184 * (7.54056 + 7.51625 * (1442.7 / temperature / math.sinh(1442.7 / temperature)) ** 2 +
                        5.38023 * (647.502 / temperature / math.cosh(647.502 / temperature)) ** 2)/44
    if gas.name == 'AIR':
        return (7 * 10 ** (-12) * temperature ** 4 - 3 * 10 ** (-8) * temperature ** 3 + \
               4 * 10 ** (-5) * temperature ** 2 - 0.0144 * temperature + 31.648)/28.9647
    if gas.name == 'H2O':
        return 4.184 * (7.97183 + 6.27078 * (2572.63 / temperature / math.sinh(2572.63 / temperature)) ** 2 +
                        2.051 * (1156.72 / temperature / math.cosh(1156.72 / temperature)) ** 2)/18
