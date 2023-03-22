from ece230 import ECE230
from moscap import MOSCAP
from constants import *
import sympy as sp
from sympy.physics.units import *

class MOSFET(MOSCAP):

    def __init__(self, T = 300 * kelvin):
        super().__init__(T)
        self.calculation_functions.extend([self.calc_mu_pP, self.calc_mu_nP, self.calc_V_DSat, self.calc_V_TN_at_V_BS])

        mosfet_quantities = {
            'V_GS': None, #Gate to Source Voltage
            'V_DS': None, #Drain to Source Voltage
            'V_BS': None, #Substrate (Bulk) to Source Voltage
            'L_ch': None, #Channel Length
            'W_ch': None, #Channel Width
            'X_ox': None, #Oxide Thickness
            'C_ox': None, #Oxide Capacitance per unit area
            'mu_nch': None, #Low-field electron Mobility in channel
            'N_aB' : None, #Substrate dopant concentration
            'phi_PMOS': None, #Contact potential difference
            'mu_nP' : None, #Electron mobility in p-type side
            'mu_pP' : None, #Hole mobility in p-type side
            'mu_nch' : None, #Electron mobility in channel #NOTE: this is set by the user
            'V_GS' : None, #Gate to Source Voltage
            'V_TN@V_BS' : None, #Threshold voltage at V_BS
            'V_DSat' : None, #Saturation Source-Drain Voltage #NOTE: this is set by the user

        }

        self.known_quantities.update(mosfet_quantities)

    def calc_mu_pP(self):
        if not self.should_calculate(needed=['N_aB']):
            return
        N_aB = remove_units(self.known_quantities['N_aB'])
        mu_pP = (49.7 + 418.3 / (1+(N_aB/(1.6e17))**(0.70))) * cm**2 / (volt * second)
        self.known_quantities['mu_pP'] = mu_pP

    def calc_mu_nP(self):
        if not self.should_calculate(needed=['N_aB']):
            return
        N_aB = remove_units(self.known_quantities['N_aB'])
        mu_nP = (232 + 1180 / (1+(N_aB/(8.0e16))**(0.90))) * cm**2 / (volt * second)
        self.known_quantities['mu_nP'] = mu_nP

    def calc_V_DSat(self):
        if not self.should_calculate(needed=['V_GS', 'V_TN@V_BS']):
            return
        V_GS = self.get_known_quantity('V_GS')
        V_TN = self.get_known_quantity('V_TN@V_BS')
        V_DSat = V_GS - V_TN
        self.known_quantities['V_DSat'] = V_DSat

    def calc_V_TN_at_V_BS(self):
        if not self.should_calculate(needed=['V_FBN', 'phi_FB', 'C_ox', 'N_aB', 'V_BS']):
            return
        V_FB = remove_units(self.known_quantities['V_FBN'])
        phi_FB = remove_units(self.known_quantities['phi_FB'])
        C_ox = remove_units(self.known_quantities['C_ox'])
        N_aB = remove_units(self.known_quantities['N_aB'])
        V_BS = remove_units(self.known_quantities['V_BS'])
        V_TN = V_FB * volt + 2 * phi_FB  * volt + 1/C_ox * sp.sqrt(2 * q_ * epsilon_Si_ * N_aB * (2 * phi_FB - V_BS)) * volt
        self.known_quantities['V_TN@V_BS'] = V_TN

