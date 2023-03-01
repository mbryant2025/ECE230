from ece230 import ECE230
from constants import *
import sympy as sp
from sympy.physics.units import *

class MOSFET(ECE230):

    def __init__(self, T = 300 * kelvin):
        super().__init__(T)
        self.calculation_functions.extend([])

        mosfet_quantities = {
            'V_GS': None, #Gate to Source Voltage
            'V_DS': None, #Drain to Source Voltage
            'V_BS': None, #Substrate (Bulk) to Source Voltage
            'L_Channel': None, #Channel Length
            'W_Channel': None, #Channel Width
            'X_ox': None, #Oxide Thickness
            'C_ox': None, #Oxide Capacitance per unit area
            'mu_nch': None, #Low-field electron Mobility in channel
            'N_aB' : None, #Substrate dopant concentration
            'phi_PMOS': None, #Contact potential difference

        }

        self.known_quantities.update(mosfet_quantities)