from ece230 import ECE230
from mosfet import MOSFET
from constants import *
import sympy as sp
from sympy.physics.units import *
import numpy as np


def parallel(Rs):
    return 1 / np.sum(1 / Rs)

#Simple MOSFET Amplifier
class MOSAMP(MOSFET):

    def __init__(self, T = 300 * kelvin):
        super().__init__(T)
        self.calculation_functions.extend([])

        mosfet_quantities = {
            'Q_V_DS' : None, #Q point value of V_DS
            'Q_I_DS' : None, #Q point value of I_DS
            'V_DD' : None, #Supply Voltage NOTE: this is set by the user
            'RD' : None, #Drain Resistance NOTE: this is set by the user
            'RS' : None, #Source Resistance NOTE: this is set by the user
            'R1' : None, #Gate to VDD Resistance NOTE: this is set by the user
            'R2' : None, #Gate to GND Resistance NOTE: this is set by the user
            'Rin' : None, #Input Resistance NOTE: this is set by the user
        }

        self.known_quantities.update(mosfet_quantities)

        print('Run calc_Q() to find the Q point. This will modify V_DS to find the Q point.')

    def calc_Q(self):

        print('Note: V_DS is modified in calc_Q() to find the Q point.')

        if not self.should_calculate(needed=['V_DD', 'RD', 'RS', 'K_n', 'V_DSat']):
            return
        V_DD = remove_units(self.known_quantities['V_DD'])
        RD = remove_units(self.known_quantities['RD'])
        RS = remove_units(self.known_quantities['RS'])

        MAX_V_DS = 15 #Assume V_DS is less than 15V for the Q point
        V_DS = np.linspace(0, MAX_V_DS, 5000)
        I_DS = V_DD / (RD + RS) - 1 / (RD + RS) * V_DS

        I_D = np.array([])
        for V in V_DS:
            self.add_known_quantity('V_DS', V * volt)
            I_D = np.append(I_D, remove_units(self.get_known_quantity('I_D')))

        V_DS_Q = V_DS[np.argmin(np.abs(I_DS - I_D))]
        I_DS_Q = I_DS[np.argmin(np.abs(I_DS - I_D))]

        self.known_quantities['Q_V_DS'] = V_DS_Q * volt
        self.known_quantities['Q_I_DS'] = I_DS_Q * ampere