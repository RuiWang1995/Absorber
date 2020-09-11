import math


# return J/mol/K
def viscosity(gas, temperature):
    if gas.name == 'CO2':
        a_0 = -8.84999 * 10 ** -7
        a_1 = 5.84534 * 10 ** -8
        a_2 = -1.88244 * 10 ** -11
        a_3 = 3.32298 * 10 ** -15
        return a_3 * temperature ** 3 + a_2 * temperature ** 2 + a_1 * temperature + a_0

    if gas.name == 'AIR':
        T0air = 524.07
        miu0air = 0.01827
        Cair = 120
        Tair = 1.8 * temperature
        aair = 0.555 * T0air + Cair
        bair = 0.555 * Tair + Cair
        return miu0air * aair / bair * (Tair / T0air) ** (3 / 2) / 1000

    if gas.name == 'H2O':
        a_0 = -6.607122 * 10 ** -6
        a_1 = 6.283711 * 10 ** -8
        a_2 = -3.955387 * 10 ** -11
        a_3 = 1.946393 * 10 ** -14
        return a_3 * temperature ** 3 + a_2 * temperature ** 2 + a_1 * temperature + a_0

