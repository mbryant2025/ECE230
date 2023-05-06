from ece230 import ECE230
from mosfet import MOSFET
from constants import *
import sympy as sp
from sympy.physics.units import *
import numpy as np

#Simple MOSFET Amplifier
class MOSINVERTER(MOSFET):

    def __init__(self, T = 300 * kelvin):
        super().__init__(T)
        self.calculation_functions.extend([self.calc_V_It, self.calc_V_OL, self.calc_I_D_max, self.calc_P_max])

        mosfet_quantities = {
            'V_DD' : None, #Supply Voltage NOTE: this is set by the user
            'RD' : None, #Drain Resistance NOTE: this is set by the user
            'V_It' : None, #Transition voltage
            'V_OL' : None, #Output Voltage Low
            'I_D_max' : None, #Maximum Drain Current
            'P_max' : None, #Maximum Power Dissipation
            }

        self.known_quantities.update(mosfet_quantities)

    def calc_V_It(self):
        if not self.should_calculate(needed=['K_n', 'RD', 'V_TB', 'V_DD']):
            return
        K_n = remove_units(self.known_quantities['K_n'])
        RD = remove_units(self.known_quantities['RD'])
        V_TB = remove_units(self.known_quantities['V_TB'])
        V_DD = remove_units(self.known_quantities['V_DD'])

        V_It = sp.Symbol('V_It')
        V_It = sp.solve(K_n * RD * (V_It - V_TB)**2 + V_It - V_TB - V_DD, V_It)
        #take positive solution (higher voltage)
        V_It = V_It[-1]
        self.known_quantities['V_It'] = V_It * volt

        #need V_GS
        if not self.should_calculate(needed=['V_GS', 'V_TB']):
            return
        V_GS = remove_units(self.known_quantities['V_GS'])
        V_TB = remove_units(self.known_quantities['V_TB'])
        if V_GS < V_TB:
            #cut off
            self.known_quantities['Bias Range'] = 'Cut Off'
        elif V_GS > V_It:
            #linear range
            self.known_quantities['Bias Range'] = 'Linear'
        else:
            #saturation range
            self.known_quantities['Bias Range'] = 'Saturation'

    def calc_V_OL(self):
        #need V_DD K_N RD and V_TB
        if not self.should_calculate(needed=['K_N', 'RD', 'V_TB', 'V_DD']):
            return
        K_N = remove_units(self.known_quantities['K_N'])
        RD = remove_units(self.known_quantities['RD'])
        V_TB = remove_units(self.known_quantities['V_TB'])
        V_DD = remove_units(self.known_quantities['V_DD'])
        VOL = V_DD / (1 + K_N * RD * (V_DD - V_TB)) * volt
        self.known_quantities['V_OL'] = VOL

    def calc_I_D_max(self):
        if not self.should_calculate(needed=['V_DD', 'V_OL', 'RD']):
            return
        V_DD = self.get_known_quantity('V_DD')
        V_OL = self.get_known_quantity('V_OL')
        RD = self.get_known_quantity('RD')
        I_D_max = (V_DD - V_OL) / RD * ampere / volt * ohm
        self.known_quantities['I_D_max'] = I_D_max

    def calc_P_max(self):
        if not self.should_calculate(needed=['I_D_max', 'V_DD']):
            return
        I_D_max = self.get_known_quantity('I_D_max')
        V_DD = self.get_known_quantity('V_DD')
        P_max = I_D_max * V_DD * watt / volt / ampere
        self.known_quantities['P_max'] = P_max
    
    