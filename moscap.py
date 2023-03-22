from ece230 import ECE230
from constants import *
import sympy as sp
from sympy.physics.units import *

class MOSCAP(ECE230):

    def __init__(self, T = 300 * kelvin):
        super().__init__(T)
        self.calculation_functions.extend([self.calc_phi_PM, self.calc_phi_NM, self.calc_phi_FB, self.calc_V_SCB_at_V_TN,
                                        self.calc_Vdep_SCB, self.calc_W_dB_at_V_TN, self.calc_L_DB, self.calc_Q_SCB_at_V_SCB,
                                        self.calc_Q_G_at_V_GB, self.calc_epsilon_xOX_at_V_GB, self.calc_V_ox, self.calc_epsilon_xB,
                                        self.calc_V_TB, self.calc_C_ox, self.calc_V_FBN, self.calc_C_substrate_at_V_FBN,
                                        self.calc_C_GB_at_V_FBN, self.calc_C_substrate_at_V_TB, self.calc_C_GB_at_V_TB,
                                        self.calc_W_dB_at_V_GB, self.calc_C_substrate_at_V_GB, self.calc_C_GBHF_DIVIDED_C_ox])

        moscap_quantities = {
            'X_ox' : None, #Gate-oxide thickness
            'N_aB' : None, #Substrate acceptor concentration
            'N_dB' : None, #Substrate donor concentration
            'Q_f' : None, #Fixed oxide charge density
            'Q_it' : None, #The interface-trap charge density
            'phi_PM' : None, #Metal-p-side contact potential difference
            'phi_NM' : None, #Metal-n-side contact potential difference
            'phi_FB' : None, #Fermi potential difference
            'V_GB' : 0 * volt, #Gate bulk potential difference
            'Vdep_SCB' : None, #Voltage across the substrate space-charge region
            'C_ox' : None, #Gate oxide capacitance
            'V_SCB@V_TN' : None, #Voltage across the substrate space-charge region at the threshold voltage
            'W_dB@V_TN' : None, #Depletion-region width at the threshold voltage
            'W_dB@V_GB' : None, #Depletion-region width at the gate bulk voltage
            'Q_SCB@V_SCB' : None, #Substrate space charge density at the substrate space-charge region voltage
            'L_DB' : None, #Extrinsic Debye length
            'Q_G@V_GB' : None, #Gate charge density at the gate bulk voltage
            'epsilon_xOX@V_GB' : None, #Oxide electric field
            'V_ox' : None, #Oxide potential difference
            'epsilon_xB' : None, #Transverse electric-field distribution
            'V_TB' : None, #Threshold voltage
            'V_FBN' : None, #Flat-band voltage
            'C_substrate@V_FBN' : None, #Small signal substrate capacitance at the flat-band voltage
            'C_GB@V_FBN' : None, #Small signal gate to bulk capacitance at the flat-band voltage
            'C_substrate@V_TB' : None, #Small signal substrate capacitance at the threshold voltage
            'C_GB@V_TB' : None, #Small signal gate to bulk capacitance at the threshold voltage
            'C_GBHF_DIVIDED_C_ox' : None, #Normalized high-frequency differential gate-to-bulk capacitance
            'C_substrate@V_GB' : None, #Small signal substrate capacitance at the gate bulk voltage

            'X_w' : None, #Wafer thickness
        }

        self.known_quantities.update(moscap_quantities)

    def calc_phi_PM(self):
        #assumes aluminum metal contact
        if not self.should_calculate(needed=['N_aB']):
            return
        N_aB = self.known_quantities['N_aB']
        phi_PM = -0.51165 * volt - k_BT_300K_DIVIDED_q * sp.log(N_aB / n_i_300K)
        self.known_quantities['phi_PM'] = phi_PM

    def calc_phi_NM(self):
        #assumes aluminum metal contact
        if not self.should_calculate(needed=['N_dB']):
            return
        N_dB = self.known_quantities['N_dB']
        phi_NM = -0.51165 * volt + k_BT_300K_DIVIDED_q * sp.log(N_dB / n_i_300K)
        self.known_quantities['phi_NM'] = phi_NM

    def calc_phi_FB(self):
        if not self.should_calculate(needed=['N_aB']):
            return
        N_aB = self.known_quantities['N_aB']
        phi_FB = k_BT_300K_DIVIDED_q * sp.log(N_aB / n_i_300K)
        self.known_quantities['phi_FB'] = phi_FB

    def calc_V_SCB_at_V_TN(self):
        if not self.should_calculate(needed=['phi_FB']):
            return
        phi_FB = self.known_quantities['phi_FB']
        V_SCB_at_V_TN = phi_FB * 2
        self.known_quantities['V_SCB@V_TN'] = V_SCB_at_V_TN

    def calc_Vdep_SCB(self):
        if not self.should_calculate(needed=['V_GB', 'V_FBN', 'N_aB', 'C_ox']):
            return
        V_GB = remove_units(self.known_quantities['V_GB'])
        V_FBN = remove_units(self.known_quantities['V_FBN'])
        N_aB = remove_units(self.known_quantities['N_aB'])
        C_ox = remove_units(self.known_quantities['C_ox'])
        Vdep_SCB = V_GB * volt - V_FBN * volt + (q_ * epsilon_Si_ * N_aB / C_ox**2) * (1 - sp.sqrt(1 + 2 * C_ox**2 * (V_GB - V_FBN) / (q_ * epsilon_Si_ * N_aB))) * volt
        self.known_quantities['Vdep_SCB'] = Vdep_SCB

    def calc_W_dB_at_V_TN(self):
        if not self.should_calculate(needed=['N_aB', 'V_SCB@V_TN']):
            return
        N_aB = remove_units(self.known_quantities['N_aB'])
        Vdep_SCB_at_V_TN = remove_units(self.known_quantities['V_SCB@V_TN'])
        W_dep_dB_at_V_TN = sp.sqrt(2 * epsilon_Si_ * Vdep_SCB_at_V_TN / (q_ * N_aB)) * cm
        self.known_quantities['W_dB@V_TN'] = W_dep_dB_at_V_TN

    def calc_L_DB(self):
        if not self.should_calculate(needed=['N_aB']):
            return
        N_aB = remove_units(self.known_quantities['N_aB'])
        L_DB = sp.sqrt(epsilon_Si_ * k_BT_300K_ / (q_**2 * N_aB)) * cm
        self.known_quantities['L_DB'] = L_DB
        

    def calc_Q_SCB_at_V_SCB(self):
        if not self.should_calculate(needed=['L_DB', 'N_aB', 'V_SCB@V_TN']):
            return
        L_DB = remove_units(self.known_quantities['L_DB'])
        N_aB = remove_units(self.known_quantities['N_aB'])
        Vdep_SCB = remove_units(self.known_quantities['V_SCB@V_TN'])
        beta = 1/(k_BT_300K_DIVIDED_q_)
        no_B = n_i_300K_**2 / N_aB
        po_B = N_aB
        inner1 = sp.exp(-beta * Vdep_SCB) + beta * Vdep_SCB - 1
        inner2 = sp.exp(beta * Vdep_SCB) - beta * Vdep_SCB - 1
        Q_SCB_at_V_SCB = -math.sqrt(2) * epsilon_Si_ * k_BT_300K_ / (q_ * L_DB) * sp.sqrt(inner1 + (no_B / po_B) * inner2) * C / cm**2
        self.known_quantities['Q_SCB@V_SCB'] = Q_SCB_at_V_SCB

    def calc_Q_G_at_V_GB(self):
        if not self.should_calculate(needed=['Q_f', 'Q_it', 'Q_SCB@V_SCB']):
            return
        Q_f = self.known_quantities['Q_f']
        Q_it = self.known_quantities['Q_it']
        Q_SCB_at_V_SCB = self.known_quantities['Q_SCB@V_SCB']
        Q_G_at_V_GB = -(Q_f + Q_it + Q_SCB_at_V_SCB)
        self.known_quantities['Q_G@V_GB'] = Q_G_at_V_GB

    def calc_epsilon_xOX_at_V_GB(self):
        if not self.should_calculate(needed=['Q_G@V_GB']):
            return
        Q_G_at_V_GB = remove_units(self.known_quantities['Q_G@V_GB'])
        epsilon_xOX = Q_G_at_V_GB / epsilon_ox_ * volt / cm
        self.known_quantities['epsilon_xOX@V_GB'] = epsilon_xOX

    def calc_V_ox(self):
        if not self.should_calculate(needed=['epsilon_xOX@V_GB', 'X_ox']):
            return
        epsilon_xOX = self.known_quantities['epsilon_xOX@V_GB']
        X_ox = self.known_quantities['X_ox']
        V_ox = epsilon_xOX * X_ox
        self.known_quantities['V_ox'] = V_ox

    def calc_epsilon_xB(self):
        if not self.should_calculate(needed=['Q_SCB@V_SCB', 'X_ox']):
            return
        Q_SCB_at_V_SCB = remove_units(self.known_quantities['Q_SCB@V_SCB'])
        epsilon_xB = -Q_SCB_at_V_SCB / epsilon_Si_ * volt / cm
        self.known_quantities['epsilon_xB'] = epsilon_xB

    def calc_V_TB(self):
        if not self.should_calculate(needed=['V_ox', 'V_SCB@V_TN', 'phi_PM']):
            return
        V_ox = self.known_quantities['V_ox']
        V_SCB_at_V_TN = self.known_quantities['V_SCB@V_TN']
        phi_PM = self.known_quantities['phi_PM']
        V_TN = V_SCB_at_V_TN + phi_PM + V_ox
        self.known_quantities['V_TB'] = V_TN

    def calc_C_ox(self):
        if not self.should_calculate(needed=['X_ox']):
            return
        X_ox = self.known_quantities['X_ox']
        C_ox = epsilon_ox / X_ox
        self.known_quantities['C_ox'] = C_ox

    def calc_V_FBN(self):
        if not self.should_calculate(needed=['phi_PM', 'Q_f', 'Q_it', 'C_ox']):
            return
        phi_PM = self.known_quantities['phi_PM']
        Q_f = self.known_quantities['Q_f']
        Q_it = self.known_quantities['Q_it']
        C_ox = self.known_quantities['C_ox']
        V_FBN = phi_PM - ((Q_f + Q_it) / C_ox) * COULOUMB_PER_FARAD_TO_VOLT
        self.known_quantities['V_FBN'] = V_FBN

    def calc_W_dB_at_V_GB(self):
        if not self.should_calculate(needed=['N_aB', 'Vdep_SCB']):
            return
        N_aB = remove_units(self.known_quantities['N_aB'])
        Vdep_SCB = remove_units(self.known_quantities['Vdep_SCB'])
        W_dB_at_V_GB = sp.sqrt(2 * epsilon_Si_ * Vdep_SCB / (q_ * N_aB)) * cm
        self.known_quantities['W_dB@V_GB'] = W_dB_at_V_GB 

    def calc_C_substrate_at_V_FBN(self):
        if not self.should_calculate(needed=['V_FBN', 'N_aB']):
            return
        self.known_quantities['C_substrate@V_FBN'] = float('inf') * farad / cm**2

    def calc_C_GB_at_V_FBN(self):
        if not self.should_calculate(needed=['C_ox']):
            return
        C_ox = self.known_quantities['C_ox']
        C_GB_at_V_FBN = C_ox
        self.known_quantities['C_GB@V_FBN'] = C_GB_at_V_FBN

    def calc_C_substrate_at_V_TB(self):
        if not self.should_calculate(needed=['W_dB@V_TN']):
            return
        W_dB_at_V_TN = self.known_quantities['W_dB@V_TN']
        C_substrate_at_V_TB = epsilon_Si / W_dB_at_V_TN
        self.known_quantities['C_substrate@V_TB'] = C_substrate_at_V_TB

    def calc_C_substrate_at_V_GB(self):
        if not self.should_calculate(needed=['W_dB@V_GB']):
            return
        W_dB_at_V_GB = self.known_quantities['W_dB@V_GB']
        C_substrate_at_V_TB = epsilon_Si / W_dB_at_V_GB
        self.known_quantities['C_substrate@V_GB'] = C_substrate_at_V_TB

    def calc_C_GB_at_V_TB(self):
        if not self.should_calculate(needed=['C_substrate@V_TB', 'C_ox']):
            return
        C_substrate_at_V_TB = self.known_quantities['C_substrate@V_TB']
        C_ox = self.known_quantities['C_ox']
        C_GB_at_V_TB = (C_ox * C_substrate_at_V_TB) / (C_ox + C_substrate_at_V_TB)
        self.known_quantities['C_GB@V_TB'] = C_GB_at_V_TB

    def calc_C_GBHF_DIVIDED_C_ox(self):
        if not self.should_calculate(needed=['C_GB@V_TB', 'C_ox']):
            return
        C_GB_at_V_TB = self.known_quantities['C_GB@V_TB']
        C_ox = self.known_quantities['C_ox']
        C_GBHF_DIVIDED_C_ox = C_GB_at_V_TB / C_ox
        self.known_quantities['C_GBHF_DIVIDED_C_ox'] = C_GBHF_DIVIDED_C_ox

