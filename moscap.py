from ece230 import ECE230
from constants import *
import sympy as sp
from sympy.physics.units import *

class MOSCAP(ECE230):

    def __init__(self, T = 300 * kelvin):
        super().__init__(T)
        self.calculation_functions.extend([])

        moscap_quantities = {
            'N_aB' : None, #Substrate acceptor concentration
            'N_dB' : None, #Substrate donor concentration
            'X_ox' : None, #Gate oxide thickness
            'X_w' : None, #Wafer thickness
            'Q_f' : None, #Fixed oxide charge density
            'phi_PM' : None, #Metal-p-side contact potential difference
        }

        self.known_quantities.update(moscap_quantities)