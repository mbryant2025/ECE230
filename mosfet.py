from ece230 import ECE230
from moscap import MOSCAP
from constants import *
import sympy as sp
from sympy.physics.units import *
import numpy as np

class MOSFET(MOSCAP):

    def __init__(self, T = 300 * kelvin):
        super().__init__(T)
        self.calculation_functions.extend([self.calc_V_TB, self.calc_mu_pP, self.calc_mu_nP, self.calc_V_DSat, self.calc_V_TN_at_V_BS, self.calc_V_DSat,
                                           self.calc_Bias_Range, self.calc_K_n, self.calc_K_N, self.calc_I_D, self.calc_V_GS, self.calc_k_prime_n,
                                           self.calc_width_to_length, self.calc_cutoff_freq])

        mosfet_quantities = {
            'Bias Range': None, #Bias Range
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
            'V_TN@V_BS' : None, #Threshold voltage at V_BS
            'V_DSat' : None, #Saturation Source-Drain Voltage
            'I_D' : None, #Drain Current
            'K_n' : None, #Electron mobility factor
            'K_N' : None, #Electron mobility factor
            'k_prime_n' : None, #Electron mobility factor
            'width_to_length' : None, #Width to Length Ratio
            'cutoff_freq' : None, #Cutoff Frequency
        }

        self.known_quantities.update(mosfet_quantities)

    def calc_V_TB(self):
        if not self.should_calculate(needed=['V_FBN', 'phi_FB', 'C_ox', 'N_aB', 'C_ox']):
            return
        V_FB = remove_units(self.known_quantities['V_FBN'])
        phi_FB = remove_units(self.known_quantities['phi_FB'])
        C_ox = remove_units(self.known_quantities['C_ox'])
        N_aB = remove_units(self.known_quantities['N_aB'])
        V_TB = (V_FB + 2 * phi_FB + 1/C_ox * sp.sqrt(4 * q_ * epsilon_Si_ * N_aB * phi_FB)) * volt
        self.known_quantities['V_TB'] = V_TB


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
        V_TN = V_FB * volt + 2 * phi_FB * volt + 1/C_ox * sp.sqrt(2 * q_ * epsilon_Si_ * N_aB * (2 * phi_FB - V_BS)) * volt
        self.known_quantities['V_TN@V_BS'] = V_TN

    def calc_V_DSat(self):
        if not self.should_calculate(needed=['V_GS', 'V_TB']):
            return
        V_GS = self.get_known_quantity('V_GS')
        V_TB = self.get_known_quantity('V_TB')
        V_DSat = V_GS - V_TB
        self.known_quantities['V_DSat'] = V_DSat

    def calc_V_GS(self):
        if not self.should_calculate(needed=['V_DS', 'V_TB']):
            return
        V_DSat = self.get_known_quantity('V_DSat')
        V_TB = self.get_known_quantity('V_TB')
        V_GS = V_DSat + V_TB
        self.known_quantities['V_GS'] = V_GS

    def calc_Bias_Range(self):
        if not self.should_calculate(needed=['V_DSat', 'V_DS', 'V_TB', 'V_GS']):
            return
        V_DSat = self.get_known_quantity('V_DSat')
        V_DS = self.get_known_quantity('V_DS')
        V_TB = self.get_known_quantity('V_TB')
        V_GS = self.get_known_quantity('V_GS')
        if V_DS > V_DSat and V_GS > V_TB:
            self.known_quantities['Bias Range'] = 'Saturation'
        else:
            self.known_quantities['Bias Range'] = 'Linear'
        if V_GS < V_TB:
            self.known_quantities['Bias Range'] = 'Cut-Off'


    def calc_K_n(self):
        if self.should_calculate(needed=['mu_nch', 'C_ox', 'W_ch', 'L_ch']):
            mu_nch = self.get_known_quantity('mu_nch')
            C_ox = self.get_known_quantity('C_ox')
            W_ch = self.get_known_quantity('W_ch')
            L_ch = self.get_known_quantity('L_ch')
            K_n = mu_nch * C_ox * W_ch / (2*L_ch)
            self.known_quantities['K_n'] = K_n
        #if have K_N, then K_n = 0.5*K_N
        if self.should_calculate(needed=['K_N']):
            K_N = self.get_known_quantity('K_N')
            self.known_quantities['K_n'] = 0.5*K_N


    def calc_K_N(self):
        if self.should_calculate(needed=['mu_nch', 'C_ox', 'W_ch', 'L_ch']):
            mu_nch = self.get_known_quantity('mu_nch')
            C_ox = self.get_known_quantity('C_ox')
            W_ch = self.get_known_quantity('W_ch')
            L_ch = self.get_known_quantity('L_ch')
            K_N = mu_nch * C_ox * W_ch / (L_ch)
            self.known_quantities['K_N'] = K_N

        # If we have I_D, Bias Range, V_DS, V_GS, V_TB, then we can calculate K_N
        if self.should_calculate(needed=['I_D', 'Bias Range', 'V_DS', 'V_GS', 'V_TB']):
            I_D = self.get_known_quantity('I_D')
            V_DS = self.get_known_quantity('V_DS')
            V_GS = self.get_known_quantity('V_GS')
            V_TB = self.get_known_quantity('V_TB')
            bias_range = self.get_known_quantity('Bias Range')
            if bias_range == 'Linear':
                K_N = I_D / ((V_GS - V_TB) * V_DS - 0.5 * V_DS**2)
                self.known_quantities['K_N'] = K_N

        #if have K_n, then K_N = 2*K_n
        if self.should_calculate(needed=['K_n']):
            K_n = self.get_known_quantity('K_n')
            self.known_quantities['K_N'] = 2*K_n

    def calc_k_prime_n(self):
        if self.should_calculate(needed=['mu_nch', 'C_ox']):
            mu_nch = self.get_known_quantity('mu_nch')
            C_ox = self.get_known_quantity('C_ox')
            k_prime_n = mu_nch * C_ox
            self.known_quantities['k_prime_n'] = k_prime_n

    def calc_I_D(self):
        if not self.should_calculate(needed=['K_n', 'V_DS', 'Bias Range', 'V_DSat']):
            return
        K_n = self.get_known_quantity('K_n')
        V_DS = self.get_known_quantity('V_DS')
        Bias_Range = self.get_known_quantity('Bias Range')
        V_DSat = self.get_known_quantity('V_DSat')
        if Bias_Range == 'Saturation':
            I_D = K_n * V_DSat**2
        elif Bias_Range == 'Linear':
            I_D = K_n*2 * (V_DSat * V_DS - 0.5 * V_DS**2)
        else:
            I_D = 0
        self.known_quantities['I_D'] = I_D

    def calc_width_to_length(self):
        if self.should_calculate(needed=['W_ch', 'L_ch']):
            W_ch = self.get_known_quantity('W_ch')
            L_ch = self.get_known_quantity('L_ch')
            width_to_length = W_ch / L_ch
            self.known_quantities['width_to_length'] = width_to_length
        #If we have K_N and k_prime_n, then we can calculate width_to_length
        if self.should_calculate(needed=['K_N', 'k_prime_n']):
            K_N = self.get_known_quantity('K_N')
            k_prime_n = self.get_known_quantity('k_prime_n')
            width_to_length = K_N / k_prime_n
            self.known_quantities['width_to_length'] = width_to_length


    def calc_cutoff_freq(self):
        #Need V_GS, V_TB, L_Ch
        if self.should_calculate(needed=['V_GS', 'V_TB', 'L_ch', 'mu_nch']):
            V_GS = self.get_known_quantity('V_GS')
            V_TB = self.get_known_quantity('V_TB')
            L_ch = self.get_known_quantity('L_ch')
            mu_nch = self.get_known_quantity('mu_nch')
            cutoff_freq = 3 * mu_nch * (V_GS - V_TB) / (4 * math.pi * L_ch**2)
            self.known_quantities['cutoff_freq'] = cutoff_freq
