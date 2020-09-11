# return W/m/K
def heat_transfer_coefficient(gas, temperature):
    if gas.name == 'CO2':
        return 0.001530777 + 1.364015 * 10 ** (-5) * temperature + 1.527046 * 10 ** (-7) * temperature ** 2 \
               - 1.075683 * 10 ** (-10) * temperature ** 3

    if gas.name == 'AIR':
        return -3.9333 * 10 ** (-4) + 1.0184 * 10 ** (-4) * temperature - 4.8574 * 10 ** (-8) * temperature ** 2 \
               + 1.5207 * 10 ** (-11) * temperature ** 3

    if gas.name == 'H2O':
        return 0
