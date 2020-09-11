import math

'''
@Author: Rui
@Date: 2020-03-26 23:45:30
@LastEditTime: 2020-03-30 04:19:09
@LastEditors: Please set LastEditors
@Description: P-R EoS
@FilePath: /ThermoProject/CH4.py
'''
R = 8.314
Na = 6.022 * 10 ** 23
kb = 1.38064852 * 10 ** -23
planck = 6.626 * 10 ** -34
naturalE = 2.71828183


class SolvePrEos(object):

    def __init__(self, gas, properties):
        self.gas = gas
        self.properties = properties

    def Tr(self):
        return self.properties.temperature / self.gas.tcritical

    def ratioB(self):
        return 0.0778 * self.properties.pressure / self.gas.pcritical / self.Tr()

    def kappa(self):
        return 0.37464 + 1.5422 * self.gas.afactor - 0.26992 * self.gas.afactor ** 2

    def alpha(self):
        return (1 + (1 - math.sqrt(self.Tr())) * self.kappa()) ** 2

    def aValue(self):
        return 0.45724 * (R * self.gas.tcritical) ** 2 / self.gas.pcritical * self.alpha()

    def bValue(self):
        return 0.0778 * R * self.gas.tcritical / self.gas.pcritical

    def solve(self):
        from sympy import solveset, S
        from sympy.abc import x
        roots = solveset(self.properties.pressure - R * self.properties.temperature / (
                x - self.bValue()) + self.aValue()  / (
                                 x ** 2 + 2 * self.bValue() * x - self.bValue() ** 2), x, domain=S.Reals)
        if len(roots) == 1:
            return roots.args[0]
        else:
            return roots.args[2]

    def compress(self):
        # return (R * self.properties.temperature / self.properties.pressure)/self.solve()
        return self.solve() / (self.solve() - self.bValue()) - self.alpha() * self.aValue() / (
            R * self.properties.temperature * self.solve() + 2 * R * self.properties.temperature * self.bValue() -
            R * self.properties.temperature * (self.bValue()) ** 2 / self.solve())
