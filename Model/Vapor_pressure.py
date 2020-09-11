import math


# return Pa
def water_vapor_pressure(temperature):
    return 61.094 * math.exp(17.625 * (temperature - 273.15) / (temperature - 30.11))
